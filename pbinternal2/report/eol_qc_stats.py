#!/usr/bin/env python

import time
import sys
import os
import logging
import numpy as np
import glob
from pbcommand.models.report import Report, Attribute, PlotGroup, Plot
from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from multiprocessing import Pool

from pbcore.io import openDataFile, SubreadSet, AlignmentSet
from pbcore.io.dataset.utils import sampleHolesUniformly

from pbinternal2 import get_version

log = logging.getLogger(__name__)
__version__ = get_version()

# ZMW CSV:

def Get(key, altName=None, transform=(lambda x: x)):
    def getter(x):
        return transform(getattr(x, key))
    if altName is None:
        altName = key
    getter.__name__ = altName
    return getter

def getPkmid(read):
    try:
        return np.array(read.peer.opt('pm'), dtype=np.int)
    except KeyError:
        raise IOError("Input Subreads BAM file must be PacBio Internal Bam")

def getPulseLabels(read):
    try:
        return np.array(list(read.peer.opt('pc')), dtype='S1')
    except KeyError:
        raise IOError("Input Subreads BAM file must be PacBio Internal Bam")

def pkmid_mean(read):
    return np.mean(read.get('pms', getPkmid))

def pkmid_channel_mean(channel):
    def midmean(read):
        pms = read.get('pms', getPkmid)
        pls = read.get('pls', getPulseLabels)
        return np.mean(pms[pls == channel])
    midmean.__name__ = 'pkmid_{c}_mean'.format(c=channel)
    return midmean

class ReadShare(object):

    def __init__(self, read):
        self.sharedstate = {}
        self.sharedstate['read'] = read

    def get(self, key, get_func=None):
        res = self.sharedstate.get(key)
        if res is None and not get_func is None:
            self.sharedstate[key] = get_func(self.sharedstate['read'])
            res = self.sharedstate.get(key)
        return res

    def __getattr__(self, key):
        return getattr(self.sharedstate['read'], key)

def eol_qc_zmw_stats(aset, outcsv, nproc=1):
    start = time.clock()
    acc = [Get('movieName', 'moviename'), Get('holeNumber', 'holenumber'),
           Get('qStart'), Get('qEnd'), Get('identity', 'concordance'),
           Get('hqRegionSnr', 'snrA', lambda x: x[0]),
           Get('hqRegionSnr', 'snrC', lambda x: x[1]),
           Get('hqRegionSnr', 'snrG', lambda x: x[2]),
           Get('hqRegionSnr', 'snrT', lambda x: x[3]),
           pkmid_mean, pkmid_channel_mean('A'), pkmid_channel_mean('C'),
           pkmid_channel_mean('G'), pkmid_channel_mean('T'),
          ]
    csv = []
    for read in aset:
        row = []
        sread = ReadShare(read)
        for fun in acc:
            row.append(fun(sread))
        csv.append(row)
    log.info("ZMW info processing time: {:}".format(time.clock() - start))
    write_csv(outcsv, [a.__name__ for a in acc], csv)
    return 0

# MOVIE CSV:

def eol_qc_movie_stats(sset, aset, outcsv, nproc=1):
    csv = []
    start = time.clock()
    header = ['moviename',
              'concordance',
              'nreads',
              'nreads_mapped',
              'nsubreads',
              'nsubreads_mapped',
              'polrl_mean',
              'polrl_std',
              'substrate_id',
              'substrate_barcode',
              'substrate_lot_number',
              'coupler_laser_power',
              'templateprepkit_barcode',
              'templateprepkit_lot_number',
              'bindingkit_barcode',
              'bindingkit_lot_number',
              'sequencingkitplate_barcode',
              'sequencingkitplate_lot_number',
              'movie_length',
              'insert_len_mean',
              'insert_len_std',
              'movie_index_in_cell',
              'pd_Empty', 'pd_Productive', 'pd_Other', 'pd_Undefined',
              'BaselineLevelMean_A',
              'BaselineLevelMean_C',
              'BaselineLevelMean_G',
              'BaselineLevelMean_T',
              'BaselineLevelStdMean_A',
              'BaselineLevelStdMean_T',
              'BaselineLevelStdMean_G',
              'BaselineLevelStdMean_C',
              'HqBasPkMidMean_A',
              'HqBasPkMidMean_T',
              'HqBasPkMidMean_G',
              'HqBasPkMidMean_C',
              ]
    # TODO (mdsmith)(7-14-2016): Clean this up, use per external-resouce
    # sts.xml accessor
    for movieName, movie in sset.movieIds.items():

        # Super sketchy: movie number in cell isn't in the metadata. You have
        # Collection Number and CellIndex, but not how many movies there are in
        # a cell...
        try:
            cellpath = os.path.dirname(sset.toExternalFiles()[0])
            cellcode = cellpath.split('/')[-1]
            runpath = '/'.join(cellpath.split('/')[:-1])
            cellname = cellcode.split('_')[-1]
            pattern = os.path.join(runpath, '_'.join(['*', cellname]))
            moviesincell = glob.glob(pattern)
            movieincell = moviesincell.index(cellpath)
        except ValueError:
            # VAlueError: the path isn't as expected, fname not in list
            movieincell = -1

        row = []
        row.append(movieName)
        # concordance
        row.append(1.0 - np.mean(np.true_divide(
            aset.index.nMM + aset.index.nIns + aset.index.nDel,
            aset.index.aEnd - aset.index.aStart)))
        # nreads
        row.append(len(set(sset.index[sset.index.qId == movie].holeNumber)))
        # nreads mapped
        row.append(len(set(aset.index[aset.index.qId == movie].holeNumber)))
        # nsubreads
        row.append(len(sset.index.qId == movie))
        # nsubreads mapped
        row.append(len(aset.index.qId == movie))
        # polrl mean
        row.append(sset.metadata.summaryStats.readLenDist.sampleMean)
        # polrl stdev
        row.append(sset.metadata.summaryStats.readLenDist.sampleStd)
        # substrate id
        row.append(sset.metadata.collections[0].cellPac.barcode[:8])
        # substrate barcode
        row.append(sset.metadata.collections[0].cellPac.barcode)
        # substrate barcode
        row.append(sset.metadata.collections[0].cellPac.lotNumber)
        # coupler laser power
        row.append(sset.metadata.collections[
            0].automation.automationParameters['CouplerLaserPower'].value)
        # templateprepkit barcode
        row.append(sset.metadata.collections[0].templatePrepKit.barcode)
        # templateprepkit barcode
        row.append(sset.metadata.collections[0].templatePrepKit.lotNumber)
        # bindingkit barcode
        row.append(sset.metadata.collections[0].bindingKit.barcode)
        # bindingkit barcode
        row.append(sset.metadata.collections[0].bindingKit.lotNumber)
        row.append(sset.metadata.collections[0].sequencingKitPlate.barcode)
        row.append(sset.metadata.collections[0].sequencingKitPlate.lotNumber)
        # movie length (minutes)
        row.append(sset.metadata.collections[
            0].automation.automationParameters['MovieLength'].value)
        # insert len
        row.append(sset.metadata.summaryStats.insertReadLenDist.sampleMean)
        row.append(sset.metadata.summaryStats.insertReadLenDist.sampleStd)
        row.append(movieincell)
        # totalloading
        row.extend(sset.metadata.summaryStats.prodDist.bins)
        # baselineleveldist
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineLevelDist']['A'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineLevelDist']['C'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineLevelDist']['G'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineLevelDist']['T'][0].sampleMean)
        # baselinestddev
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineStdDist']['A'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineStdDist']['C'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineStdDist']['G'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'BaselineStdDist']['T'][0].sampleMean)
        # pkmid
        row.append(sset.metadata.summaryStats.channelDists[
            'HqBasPkMidDist']['A'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'HqBasPkMidDist']['C'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'HqBasPkMidDist']['G'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
            'HqBasPkMidDist']['T'][0].sampleMean)
        csv.append(row)
    log.info("Movie info processing time: {:}".format(time.clock() - start))
    write_csv(outcsv, header, csv)
    return 0

def write_csv(fname, header, csv):
    print fname
    log.info("Writing {:}".format(fname))
    with open(fname, 'w') as fh:
        fh.write(','.join(header))
        fh.write('\n')
        for row in csv:
            fh.write(','.join(map(str, row)))
            fh.write('\n')

def sample_first(aset, num):
    hn = np.unique(aset.index.holeNumber)[num]
    aset.filters.addRequirement(zm=[('<=', hn)])

def sample_random(aset, num):
    hns = np.random.choice(aset.index.holeNumber, num, False)
    aset.filters.addRequirement(zm=[('in', hns)])

def sample_uniform(aset, num, faillimit=25):
    possible = set(aset.index.holeNumber)
    hns = sampleHolesUniformly(num, possible,
                               faillimit=faillimit)
    olen = len(hns)
    hns = np.unique(hns)
    if len(hns) < olen:
        log.warn("Overlapping samples found: {:}".format(olen - len(hns)))
    aset.filters.addRequirement(zm=[('in', hns)])

def eol_qc_stats(sset, aset, zmwscsv, moviescsv, nholes):
    sset = SubreadSet(sset)
    aset = AlignmentSet(aset)
    try:
        sset[0].peer.opt('pm')
    except KeyError:
        raise IOError("Input Subreads BAM file must be PacBio Internal Bam")
    eol_qc_movie_stats(sset, aset, moviescsv)
    sample_uniform(aset, nholes)
    eol_qc_zmw_stats(aset, zmwscsv)
    return 0

def run(args):
    # open files so that isn't counted towards the movie analysis
    args.subreadset.index
    args.alignmentset.index
    log.info("Starting movie analysis")
    eol_qc_movie_stats(args.subreadset, args.alignmentset,
                       '.'.join([args.outprefix, 'movies', 'csv']))
    args.sampler(args.alignmentset, args.nreads)
    log.info("Starting zmw analysis")
    eol_qc_zmw_stats(args.alignmentset,
                     '.'.join([args.outprefix, 'zmws', 'csv']))
    return 0

def get_parser():
    sample_picker = {'first': sample_first,
                     'random': sample_random,
                     'uniform': sample_uniform}
    p = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__,
        default_level="WARN")
    p.add_argument('subreadset', type=SubreadSet,
                   help="Input SubreadSet for an Internal BAM")
    p.add_argument('alignmentset', type=AlignmentSet,
                   help="Input AlignmentSet for the SubreadSet")
    p.add_argument('outprefix', type=str,
                   help="Output prefix for csvs")
    p.add_argument('--nreads', type=int,
                   help="The number of reads to process")
    p.add_argument('--sampler', default='uniform',
                   type=lambda x: sample_picker.get(x),
                   choices=sample_picker.keys(),
                   help="Read sampling mechanism")
    p.add_argument('--search', type=int, default=25,
                   help=('Limit the number of hns to search for a '
                         'local hit'))
    return p

def main(argv=sys.argv):
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=run,
        alog=log,
        setup_log_func=setup_log)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

