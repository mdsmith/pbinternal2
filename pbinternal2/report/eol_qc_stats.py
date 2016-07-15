#!/usr/bin/env python

import time
import sys
import logging
import numpy as np
from pbcommand.models.report import Report, Attribute, PlotGroup, Plot
from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from functools import partial
from multiprocessing import Pool
from pbinternal2 import get_version

from pbcore.io import openDataFile, SubreadSet, AlignmentSet
from pbcore.io.dataset.utils import sampleHolesUniformly

log = logging.getLogger(__name__)
__version__ = get_version()

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

# READ ACCUMULATORS:

def holenumber(read):
    return read.holeNumber

def moviename(read):
    return read.movieName

def concordance(read):
    return read.identity

def alnlen(read):
    return read.nM + read.nMM + read.nDel + read.nIns

def readspan(read):
    return read.aEnd - read.aStart

def refspan(read):
    return read.tEnd - read.tStart

def snrA(read):
    return read.hqRegionSnr[0]

def snrC(read):
    return read.hqRegionSnr[1]

def snrG(read):
    return read.hqRegionSnr[2]

def snrT(read):
    return read.hqRegionSnr[3]

def getter(attrib, x):
    return getattr(x, attrib)

def Get(key):
    func = partial(getter, key)
    func.__name__ = key
    return func

def getPkmid(read):
    return np.array(read.peer.opt('pm'), dtype=np.int)

def getPulseLabels(read):
    return np.array(list(read.peer.opt('pc')), dtype='S1')

class ReadShare(object):

    sharedstate = {}

    def __init__(self, read):
        # it will either be None if not yet populated, or perhaps populated by
        # a different read
        if not self.sharedstate.get('read') is read:
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

def pkmid_mean(read):
    return np.mean(read.get('pms', getPkmid))

def pkmid_channel_mean(channel):
    def midmean(read):
        pms = read.get('pms', getPkmid)
        pls = read.get('pls', getPulseLabels)
        assert len(pms) == len(pls)
        return np.mean(pms[pls == channel])
    midmean.__name__ = 'pkmid_{c}_mean'.format(c=channel)
    return midmean

# MOVIE ACCUMULATORS

def movienames(sset, aset):
    return sset.readGroupTable.MovieName


def write_csv(fname, header, csv):
    with open(fname, 'w') as fh:
        fh.write(','.join(header))
        fh.write('\n')
        for row in csv:
            fh.write(','.join(map(str, row)))
            fh.write('\n')

def process_read(accs, read):
    row = []
    for fun in accs:
        row.append(fun(read))
    return row

def eol_qc_zmw_stats(aset, outcsv, nproc=1):
    start = time.clock()
    acc = [moviename, holenumber, Get('qStart'), Get('qEnd'), concordance,
           snrA, snrC, snrG, snrT, pkmid_mean,
           pkmid_channel_mean('A'),
           pkmid_channel_mean('C'),
           pkmid_channel_mean('G'),
           pkmid_channel_mean('T'),
          ]
    csv = []
    for read in aset:
        row = []
        sread = ReadShare(read)
        for fun in acc:
            row.append(fun(sread))
        csv.append(row)
    log.info("ZMW info processing time: {:}".format(time.clock() - start))
    write_csv('.'.join([outcsv, 'zmws', 'csv']),
              [a.__name__ for a in acc], csv)
    return 0

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
    for movieName, movie in sset.movieIds.items():
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
    write_csv('.'.join([outcsv, 'movies', 'csv']),
              header, csv)
    return 0

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
    eol_qc_movie_stats(sset, aset, moviescsv)
    sample_uniform(aset, nholes)
    if zmwscsv.endswith('.csv'):
        zmwscsv = zmwscsv[:-4]
    if moviescsv.endswith('.csv'):
        moviescsv = moviescsv[:-4]
    eol_qc_zmw_stats(aset, zmwscsv)

def run(args):
    log.info("Starting movie analysis")
    eol_qc_movie_stats(args.subreadset, args.alignmentset, args.outprefix)
    args.sampler(args.alignmentset, args.nreads)
    log.info("Starting zmw analysis")
    eol_qc_zmw_stats(args.alignmentset, args.outprefix)
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

