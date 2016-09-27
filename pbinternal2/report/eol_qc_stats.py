#!/usr/bin/env python

import time
import sys
import os
import logging
import numpy as np
from math import atan,pi
import re
from glob import glob
import functools

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbcore.io import openDataFile, SubreadSet, AlignmentSet
from pbcore.io.dataset.utils import (sampleHolesUniformly, quadratic_expand,
                                     xy_to_hn, hn_to_xy)

from pbinternal2 import get_version

log = logging.getLogger(__name__)
__version__ = get_version()

# ZMW CSV:

def or_empty_string(f):
    """Intended for values that will error if "internal mode" isn't used. See ITG-107. """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except (IOError, ValueError):
            return ""
    return wrapper

def Get(key, altName=None, transform=(lambda x: x)):
    def getter(x):
        return transform(getattr(x, key))
    if altName is None:
        altName = key
    getter.__name__ = altName
    return getter

def Todo(key):
    """Values that need to be calculated. See ITG-85"""
    def getter(x):
        return ''
    getter.__name__ = key
    return getter

# regex to match
pattern = re.compile("updating trc.h5 matrix with measured matrix: ([+-]?\d*\.\d*) ([+-]?\d*\.\d*) ([+-]?\d*\.\d*) ([+-]?\d*\.\d*)")

def spectralAngle(a, b):
    return atan(a/b)*180/pi

def calcSpectralAngles(movieDir):
    """Calculates the spectral angles EOL QC uses
    The story behind why these angles are calculated this way and how they came to be an important metric in the final
    chip grading is a mystery. Followup discussion is in https://jira.pacificbiosciences.com/browse/ITG-93

    :param movieDir: Directory with all of the movie files exported by the instrument
    :return: green, red and final spectral angles
    """
    for filename in [f for f in glob('%s/traceSplit/*.log' % movieDir)]:
        with open(filename) as file:
            for line in file:
                for match in re.finditer(pattern, line):
                    gG,rG,gR,rR = [float(num) for num in match.groups()]
                    greenAngle = spectralAngle(gR, gG)
                    redAngle = spectralAngle(rR, rG)
                    return greenAngle, redAngle, redAngle - greenAngle
    return '', '', ''

def formatSpectralAngles(movieDir):
    return map(lambda x: ('%.2f' % x if isinstance(x, float) else x), calcSpectralAngles(movieDir))


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

@or_empty_string
def pkmid_mean(read):
    return np.mean(read.get('pms', getPkmid))

def pkmid_channel_mean(channel):
    @or_empty_string
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
    # Order comes from Remy's list in ITG-85
    acc = [Get('movieName', 'moviename'),
           Get('holeNumber', 'holenumber'),
           Get('qStart'),
           Get('qEnd'),
           Get('identity', 'concordance'),
           Todo('nread_mapped'),
           Get('hqRegionSnr', 'snrA', lambda x: x[0]),
           Get('hqRegionSnr', 'snrC', lambda x: x[1]),
           Get('hqRegionSnr', 'snrG', lambda x: x[2]),
           Get('hqRegionSnr', 'snrT', lambda x: x[3]),
           pkmid_mean,
           pkmid_channel_mean('A'),
           pkmid_channel_mean('C'),
           pkmid_channel_mean('G'),
           pkmid_channel_mean('T'),
           Todo('BaselineLevelMean_A'),
           Todo('BaselineLevelMean_C'),
           Todo('BaselineLevelMean_G'),
           Todo('BaselineLevelMean_T'),
           Todo('pd_Empty'),
           Todo('pd_Productive'),
           Todo('pd_Other'),
           Todo('pd_Undefined'),
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

def segment(dims=(10, 8)):
    """Return a list of ndarrays of holenumbers corresponding to 2D blocks of
    zmws of the specified size

    Args:
        dims: tuple of integer dimensions of blocks, cols x rows
    Returns:
        list of ndarrays of holenumbers
    """
    rowstart = 64
    rowend = 1024
    colstart = 64
    colend = 1144

    segments = []
    for rul in range(rowstart, rowend, dims[1]):
        for cul in range(colstart, colend, dims[0]):
            block = []
            for r in range(rul, rul + dims[1]):
                for c in range(cul, cul + dims[0]):
                    if c <= colend and r <= rowend:
                        hn = xy_to_hn(c, r)
                        block.append(hn)
            segments.append(np.array(block))
    return segments


def count_aligns(aset, hns):
    """
    Args:
        aset: the AlignmentSet
        hns: a list of ndarrays of holenumbers
    Returns:
        A list of counts of holenumbers within each block that contain one or
        more alignments
    """
    hits = set(aset.index.holeNumber)
    alignments = []
    for block in hns:
        alignments.append(sum(1 if hn in hits else 0 for hn in block))
    return np.array(alignments)


def blocked_loading(aset, dims=(10, 8)):
    chunks = segment(dims=dims)
    z = count_aligns(aset, chunks)
    return z, dims

@or_empty_string
def loading_efficiency(aset):
    # this was taken from the R code:
    # but there is currently a small discrepancy
    # getLoadingEfficiency = function( z, N )
    # {
    #   pol_pM = z / (N[0] * n[1])
    #   maxConc = floor( 3 / min( pol_pM[ pol_pM > 0 ] ) )
    #   conc = seq( 1, maxConc, 1 )
    #   lambda = pol_pM %o% conc
    #   single = lambda * exp( -lambda )
    #   total = colMeans( single )  # assume uniform
    #   100 * max( total, na.rm = TRUE ) * exp(1)
    # }

    z, N = blocked_loading(aset)
    pol_pm = np.true_divide(z, np.product(N))
    maxConc = int(np.true_divide(3, min(pol_pm[pol_pm > 0])))
    conc = range(1, maxConc + 1)
    lambdaVal = np.outer(pol_pm, conc)
    single = lambdaVal * np.exp(-1.0 * lambdaVal)
    total = np.mean(single, axis=0)
    penalty = 100.0 * np.nanmax(total) * np.exp(1.0)

    # debugging output
    #print np.histogram(z, bins=range(24))
    #y = z[z >= 1]
    #print sum(y)
    #print 1.0 - np.true_divide(len(y), len(z))
    #print np.mean(y)
    #print np.var(y)
    #dispersion = np.var(y)/np.mean(y) - 1.0
    #print dispersion
    #print lambdaVal
    #print single
    #print total
    #print np.nanmax(total)
    #print penalty

    return penalty

def eol_qc_movie_stats(sset, aset, outcsv, nproc=1):
    csv = []
    start = time.clock()
    # Rearranged as per Remy's list in ITG-85: https://jira.pacificbiosciences.com/browse/ITG-85
    header = ['substrate_id',
              'substrate_barcode',
              'substrate_lot_number',
              'moviename',
              'movie_length',
              'movie_index_in_cell',
              'coupler_laser_power',
              'templateprepkit_barcode',
              'templateprepkit_lot_number',
              'bindingkit_barcode',
              'bindingkit_lot_number',
              'sequencingkitplate_barcode',
              'sequencingkitplate_lot_number',
              'concordance',
              'nreads',
              'nreads_mapped',
              'nsubreads',
              'nsubreads_mapped',
              'loading_uniformity',
              # Note that these next 4 are added from an extend call
              'pd_Empty',
              'pd_Productive',
              'pd_Other',
              'pd_Undefined',
              'polrl_mean',
              'polrl_std',
              'insert_len_mean',
              'insert_len_std',
              'BaselineLevelMean_A',
              'BaselineLevelMean_C',
              'BaselineLevelMean_G',
              'BaselineLevelMean_T',
              'BaselineLevelStdMean_A',
              'BaselineLevelStdMean_C',
              'BaselineLevelStdMean_G',
              'BaselineLevelStdMean_T',
              'HqBasPkMidMean_A',
              'HqBasPkMidMean_C',
              'HqBasPkMidMean_G',
              'HqBasPkMidMean_T',
              'HqBasPkMidMed_A',
              'HqBasPkMidMed_C',
              'HqBasPkMidMed_G',
              'HqBasPkMidMed_T',
              'SnrDist_A',
              'SnrDist_C',
              'SnrDist_G',
              'SnrDist_T',
              'Green Angle',
              'Red Angle',
              'Spectral Angle',
              'ICS Version',
              'Signal Processing Version',
              'pbinternal2 Version',
              'Instrument',
              'Sample Well Name',
              'sts.xml Windows',
              'sts.xml POSIX',
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
            moviesincell = glob(pattern)
            movieincell = moviesincell.index(cellpath)
        except ValueError:
            # VAlueError: the path isn't as expected, fname not in list
            movieincell = -1

        row = []
        # substrate id
        row.append(sset.metadata.collections[0].cellPac.barcode[:8])
        # substrate barcode
        row.append(sset.metadata.collections[0].cellPac.barcode)
        # substrate lot number
        row.append(sset.metadata.collections[0].cellPac.lotNumber)
        # moviename
        row.append(movieName)
        # movie length (minutes)
        row.append(sset.metadata.collections[
                       0].automation.automationParameters['MovieLength'].value)
        # movie_index_in_cell
        row.append(movieincell)
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
        # loading uniformity:
        row.append(loading_efficiency(aset))
        # totalloading - aka
        row.extend(sset.metadata.summaryStats.prodDist.bins)
        # polrl mean
        row.append(sset.metadata.summaryStats.readLenDist.sampleMean)
        # polrl stdev
        row.append(sset.metadata.summaryStats.readLenDist.sampleStd)
        # insert len
        row.append(sset.metadata.summaryStats.insertReadLenDist.sampleMean)
        row.append(sset.metadata.summaryStats.insertReadLenDist.sampleStd)
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
        # pkmedian
        row.append(sset.metadata.summaryStats.channelDists[
                       'HqBasPkMidDist']['A'][0].sampleMed)
        row.append(sset.metadata.summaryStats.channelDists[
                       'HqBasPkMidDist']['C'][0].sampleMed)
        row.append(sset.metadata.summaryStats.channelDists[
                       'HqBasPkMidDist']['G'][0].sampleMed)
        row.append(sset.metadata.summaryStats.channelDists[
                       'HqBasPkMidDist']['T'][0].sampleMed)
        # SNR
        row.append(sset.metadata.summaryStats.channelDists[
                       'SnrDist']['A'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
                       'SnrDist']['C'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
                       'SnrDist']['G'][0].sampleMean)
        row.append(sset.metadata.summaryStats.channelDists[
                       'SnrDist']['T'][0].sampleMean)
        row.extend(formatSpectralAngles(os.path.dirname(sset.toExternalFiles()[0])))
        # inst control version
        row.append(sset.metadata.collections[0].instCtrlVer)
        # signal processing version
        row.append(sset.metadata.collections[0].sigProcVer)
        # version of this code
        row.append(__version__)
        # instrument used. e.g. 54097
        row.append(sset.metadata.collections[0].instrumentName)
        # sample name. e.g. VP52047-14_4754_BA020367_54097_2kLambda_200pM_LPTitration_Cell-5_12mW_Magbead
        row.append(sset.metadata.collections[0].wellSample.name)
        # file path to sts.xml. Windows then POSIX styles
        sts_path = '%s.sts.xml' % sset.fileNames[0].replace('.subreadset.xml', '')
        row.append('\\%s' % sts_path.replace('/', '\\'))
        row.append(sts_path)

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
        log.info("Input Subreads BAM file must be PacBio Internal Bam. Related metrics will appear as blank values.")
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

