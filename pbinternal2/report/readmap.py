#!/usr/bin/env python
"""Generates a 'Read Map', graphically depicting the contents of each read.
Uses a relatively intelligent sorting algorithm to highlight interesting
patterns."""

from __future__ import division
import argparse
import random
import sys
import os
import numpy as np
#from scipy.stats import truncnorm, t, chi2
import time
import math
import logging
import matplotlib
import datetime
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log

from pbinternal2 import get_version
from pbcore.io import openIndexedAlignmentFile, DataSet
from pbinternal2.util.Range import Range, Ranges
from pbcommand.validators import (validate_file,
                                  validate_fofn,
                                  validate_output_dir)

log = logging.getLogger(__name__)
__version__ = get_version()

# Constants (cmap color values)
CMAP = matplotlib.cm.RdGy

TYPE_HQ = 160
TYPE_SUBREAD = 30
TYPE_ADAPTER = 255
TYPE_EMPTY = 128

COLOR_BLACK = 255

COLOR_MIN = 0
COLOR_MAX = 255

#MAX_ZMW_SIZE = 5000
MAX_NUM_ZMW = 20000
# for the resampling procedure, spurious artifacts can arise if
# the binning is too fine-grained
NUM_RESAMPLE_BINS = 50
MIN_RESAMPLE_BIN_SIZE = 5

PLOT_WIDTH = 1600
PLOT_HEIGHT = 1000

FIGURE_WIDTH_R = 0.5
FIGURE_HEIGHT_R = 0.89
FIGURE_TOP_R = 0.05
FIGURE_BOTTOM_R = 0.06
FIGURE_LEFT_R = 0.05
FIGURE_RIGHT_R = 0.05
FIGURE_X_TITLE = 0.00

COLUMN_WIDTH_R = 0.1
COLUMN_HEIGHT_R = FIGURE_HEIGHT_R
COLUMN_SPACE_R = 0.025

LEGEND_X_SPACING = 5
LEGEND_X = PLOT_WIDTH * (FIGURE_WIDTH_R + FIGURE_LEFT_R) + LEGEND_X_SPACING
LEGEND_Y = PLOT_HEIGHT * (FIGURE_BOTTOM_R)
LEGEND_WIDTH = 30
LEGEND_HEIGHT = 120
LEGEND_Y_OFFSET = 8
LEGEND_LABELS = ["No/Bad Sequence", "High Quality", "Adapter", "Subread"]
LEGEND_TYPES = [TYPE_EMPTY, TYPE_HQ, TYPE_ADAPTER, TYPE_SUBREAD]

BAR_PADDING_R = 0.05
BAR_LEFT_R = FIGURE_LEFT_R + FIGURE_RIGHT_R + FIGURE_WIDTH_R + BAR_PADDING_R
BAR_WIDTH_R = 0.35
BAR_BOTTOM_R = BAR_PADDING_R + 2 * FIGURE_HEIGHT_R / 3
BAR_HEIGHT_R = FIGURE_HEIGHT_R / 3
BAR_COLOR = '0.3'

INSERTS_PAD_R = BAR_PADDING_R
INSERTS_LEFT_R = (FIGURE_LEFT_R + FIGURE_WIDTH_R + FIGURE_RIGHT_R +
                  INSERTS_PAD_R)
INSERTS_BOTTOM_R = FIGURE_BOTTOM_R + 0.25 * FIGURE_HEIGHT_R
INSERTS_WIDTH_R = 0.33
INSERTS_HEIGHT_R = FIGURE_HEIGHT_R / 3

HELP_TEXT_PADDING = 16
HELP_TEXT_X = LEGEND_X + LEGEND_WIDTH * 7
HELP_TEXT_Y = PLOT_HEIGHT * FIGURE_BOTTOM_R
HELP_TEXT_HEIGHT = 15
HELP_TEXT_WIDTH = 520
HELP_TEXT_SIZE = 14
HELP_TEXT_COLOR = 135

MIN_BAR_WIDTH = 7

X_TICK_ROTATION = 30
X_TICK_SIZE = 12
Y_TICK_SIZE = 12
NUM_X_TICKS = 4

X_LABEL_SIZE = 14
Y_LABEL_SIZE = 14
PCT_LABEL_SIZE = 14
BAR_LABEL_SIZE = 14
LGD_LABEL_SIZE = 14

Y_TICK_INC = 1000

DPI = 72

MAX_SECONDS_PER_ZMW = 0.01


class ReadMapReport(object):

    def __init__(self, infn, outfn, iterations=100,
                 max_num_zmw=MAX_NUM_ZMW, resample=False, link=None):
        self.infn = infn
        self.outfn = outfn
        #self.fofn = fofn
        self.iterations = iterations
        self.max_num_zmw = max_num_zmw
        self.resample = resample
        self.link = link

        self._max_zmw_size = 5000
        self._regions = {}
        self._hq = {}
        self._downsample_rate = 1.0
        self._empty_zmws = 0
        self.fig = plt.figure(figsize=(PLOT_WIDTH / 72, PLOT_HEIGHT / 72))

    def run(self):
        """Executes the body of the script."""
        self._load_infile()
        self._load_rgns()
        self._downsample_rgns()
        self._split_subreads()
        self._generate_insert_distributions()
        self._generate_zmw_map()

        if self._map.shape[0] > 0:
            self._sort_zmw_map()
            self._plot_zmw_map()
            self._plot_legend()
            self._plot_loading_bars()
            self._plot_insert_distributions()
        else:
            self.fig.text(0.2, 0.2, "Error: No Region information detected.")

        self._plot_help()
        self._write_out_figure()

        return 0

    def _load_rgns(self):
        """Generates _regions, a list of (range, type) keyed on (movie,
        holeNumber) and _scores, ReadScores keyed on (movie, holeNumber)"""
        rgn_fns = []
        if self.infn.endswith('xml'):
            rgn_fns.extend(DataSet(self.infn).toFofn(uri=False))
        else:
            rgn_fns = [self.infn]
        for rgn_fn in rgn_fns:
            log.debug("Loading region information from file %s" % rgn_fn)
            self._load_rgn(rgn_fn.strip())

    def _load_rgn(self, region_fn):
        with openIndexedAlignmentFile(region_fn) as in_file:
            for row in in_file:
                hole_number = row.HoleNumber
                movie_name = row.movieName
                zmw = (movie_name, hole_number)
                #if rgnType == "HQRegion":
                if row.aEnd == row.aStart:
                    self._empty_zmws += 1
                    continue
                if zmw in self._hq:
                    self._hq[zmw].union(Range(row.aStart, row.aEnd))
                else:
                    self._hq[zmw] = Range(row.aStart, row.aEnd)

    def _downsample_rgns(self):
        """We randomly downsample the regions to a manageable number to avoid
        using too much memory. The amount is recorded to adjust summary stats
        later."""
        n_regions = len(self._hq)
        if n_regions > self.max_num_zmw:
            self._downsample_rate = n_regions / self.max_num_zmw
            log.info("Downsampling Reads from %d to %d (%.2f)" %
                     (n_regions, self.max_num_zmw, self._downsample_rate))
            for key in random.sample(self._hq.keys(),
                                     n_regions - self.max_num_zmw):
                del self._hq[key]
                if key in self._regions:
                    del self._regions[key]

    def _split_subreads(self):
        """In order to support LIMS Template workflows, this step identifies
        subreads which overlap with Adapters. It then splits the subread by the
        adapter, generating two subreads to replace the original."""
        log.info("Splitting subreads by adapters (only relevant for unrolled "
                 "templates).")
        self._is_unrolled = False
        for zmw in self._regions:
            rgns = self._regions[zmw]
            subread_ranges = Ranges()
            new_subreads = []
            subs = [rgn for rgn in rgns if rgn.type == TYPE_SUBREAD]
            for subread in subs:
                subread_ranges.add_range(subread)
            n_subreads = len(subread_ranges._ranges)
            for adapter in [rgn for rgn in rgns if rgn.type == TYPE_ADAPTER]:
                subread_ranges.remove_range(adapter)
            if len(subread_ranges._ranges) == n_subreads:
                new_subreads.extend(subs)
            else:
                self._is_unrolled = True
                logging.debug("Detected unrolled template")
                for r in subread_ranges:
                    #r.advance = len(r)
                    tqr = 1.0
                    # super hacky for now since there's an object model that
                    # should be here
                    for r2 in subs:
                        if r.intersects(r2):
                            r.strand = r2.strand
                            tqr = float(r2.advance) / float(len(r2))
                            break
                    r.advance = int(tqr * float(len(r)))
                    new_subreads.append(r)

            for nsr in new_subreads:
                nsr.type = TYPE_SUBREAD
            self._regions[zmw] = (new_subreads + [rgn for rgn in rgns
                                                  if rgn.type == TYPE_ADAPTER])

    def _load_infile(self):
        """Loads in the subread information from the .bam.  Stores it in
        _subreads as a list of ( Range, RefRange ) indexed by ( movieName,
        holeNumber )."""
        infns = []
        if self.infn.endswith('xml'):
            infns.extend(DataSet(self.infn).toFofn(uri=False))
        else:
            infns = [self.infn]
        for infn in infns:
            with openIndexedAlignmentFile(infn) as infile:
                for aln in infile:
                    movie_name = aln.movieName
                    zmw = (movie_name, aln.HoleNumber)
                    self._regions.setdefault(zmw, [])
                    pbrange = Range(aln.aStart, aln.aEnd)
                    pbrange.type = TYPE_SUBREAD
                    pbrange.strand = 1 if aln.isReverseStrand else 0
                    try:
                        pbrange.advance = aln.tEnd - aln.tStart
                    except TypeError:
                        pbrange.advance = 0
                    self._regions[zmw].append(pbrange)
                    # Add adapter info, inferred from cx tag data (presence of
                    # adapter before and adapter after)
                    # This all assumes that reads for this ZMW are processed in
                    # order and are non-overlapping
                    try:
                        adapterFlag = aln.peer.opt("cx")
                        adapterBefore = adapterFlag & 1
                        adapterAfter = adapterFlag & 2
                        prevRegions = self._regions[zmw]
                        if adapterBefore:
                            # At this point this subread has already been
                            # added, so the last 'after' adapter is two back in
                            # the list
                            connected = False
                            if len(prevRegions) > 2:
                                # try to connect the two
                                prev = prevRegions[-2]
                                if prev.type == TYPE_ADAPTER:
                                    prev.end = aln.aStart - 1
                                    connected = True
                            if not connected:
                                # add a stub
                                pbrange = Range(aln.qStart, aln.aStart)
                                pbrange.type = TYPE_ADAPTER
                                # add the stub before the subread range so it
                                # isn't picked up as the after adapter if the
                                # read doesn't have a real after adapter
                                self._regions[zmw].insert(-1, pbrange)
                        if adapterAfter:
                            # add a stub
                            pbrange = Range(aln.aEnd, aln.qEnd)
                            pbrange.type = TYPE_ADAPTER
                            self._regions[zmw].append(pbrange)
                    except KeyError:
                        # Adapter information not found.
                        pass

    def _generate_insert_distributions(self):
        #
        # We need to rebalance the ZMW distribution based on the
        # the survival curve to get the best results out of this procedure
        #
        zmws = self._regions.keys()
        max_seconds = MAX_SECONDS_PER_ZMW * len(zmws) + 5
        start_time = datetime.datetime.now()
        if self.resample:
            unrolled = np.zeros(len(self._regions), 'i')
            for i in xrange(len(zmws)):
                rgns = self._regions[zmws[i]]
                rgns.sort(key=lambda r: r.start)
                unrolled[i] = rgns[-1].end - rgns[0].start
                log.debug("Unrolled: %d", unrolled[i])
            umax = max(unrolled)
            bin_size = max(50, int(umax / NUM_RESAMPLE_BINS))
            log.debug("Using bin size of %d", bin_size)
            bins = np.arange(0, umax + bin_size, bin_size)
            h = [[] for i in xrange(len(bins))]
            for i in xrange(len(zmws)):
                k = bins.searchsorted(unrolled[i]) - 1
                if k > 0:
                    h[k].append(i)
            filled_bins = np.array([bins[i] for i, h2 in enumerate(h) if
                                    len(h2) > 0], 'i')
            h = [h2 for h2 in h if len(h2) > 0]
            for b, h2 in zip(filled_bins, h):
                log.debug("[%d] : %s" % (b, str(h2)))
            resampled_regions = []
            last_rrl = -1
            while len(resampled_regions) < len(zmws):
                nrr = len(resampled_regions)
                if nrr == last_rrl and (datetime.datetime.now() -
                                        start_time).seconds > max_seconds:
                    log.warn("Bailing out after hitting resampling time "
                             "limit (%.2fs)" % max_seconds)
                    break
                last_rrl = nrr
                j = random.randint(0, len(h) - 1)
                # random.choice(bins)
                random_length = bins[random.randint(0, len(bins) - 1)]
                j = filled_bins.searchsorted(random_length)
                if j >= len(h) or len(h[j]) < MIN_RESAMPLE_BIN_SIZE:
                    continue
                i = random.randint(0, len(h[j]) - 1)
                resampled_regions.append(h[j][i])
                #resampled_regions.append( random.choice( h[j] ) )
                #logging.debug("Resampled: %d" % unrolled[h[j][i]] )
        else:
            resampled_regions = range(0, len(zmws))
            for rgns in self._regions.values():
                rgns.sort(key=lambda r: r.start)

        censored, uncensored = [], []
        #for zmw, rgns in self._regions.iteritems( ):
        for i in resampled_regions:
            rgns = self._regions[zmws[i]]
            pattern = [r.strand if r.type ==
                       TYPE_SUBREAD else 'A' for r in rgns]
            log.debug("Pattern: %s" % "".join(map(str, pattern)))
            log.debug("Lengths: %s" % ",".join(map(
                str, [len(rgn) for rgn in rgns if rgn.type == TYPE_SUBREAD])))
            log.debug("Advances: %s" % ",".join(map(
                str, [rgn.advance for rgn in rgns
                      if rgn.type == TYPE_SUBREAD])))

            s_pattern = "".join(map(str, pattern))
            if s_pattern in ["0", "1", "0A", "1A", "A0", "A1", "0A1", "1A0"]:
                censored.append(int(max(
                    [rgn.advance for rgn in rgns
                     if rgn.type == TYPE_SUBREAD])))
                log.debug("Censored: %s" % str(int(max(
                    [rgn.advance for rgn in rgns
                     if rgn.type == TYPE_SUBREAD]))))
            #
            # special case unrolled references since their (perceived)
            # orientations don't end up fitting the canonical model
            #
            elif self._is_unrolled:
                if s_pattern in ["0A0", "1A1"]:
                    censored.append(int(max(
                        [rgn.advance for rgn in rgns
                         if rgn.type == TYPE_SUBREAD])))
                    log.debug("Censored: %s" % str(int(max(
                        [rgn.advance for rgn in rgns
                         if rgn.type == TYPE_SUBREAD]))))
                else:
                    observed = []
                    for i, rgn in enumerate(rgns):
                        if (i > 0 and (i + 1) < len(pattern) and
                                pattern[i - 1] == 'A' and
                                pattern[i + 1] == 'A' and
                                rgn.type == TYPE_SUBREAD):
                            observed.append(rgn.advance)
                    log.debug("Observed: %s" % ",".join(map(str, observed)))
                    if len(observed) > 0:
                        uncensored.append(int(np.mean(observed)))
                        log.debug("Uncensored: %s" % str(
                            int(np.mean(observed))))
            elif s_pattern in ["A0A", "0A1A", "1A0A", "A0A1", "A1A0"]:
                observed = []
                for i, rgn in enumerate(rgns):
                    if (i > 0 and i + 1 < len(pattern) and
                            pattern[i - 1] == 'A' and pattern[i + 1] == 'A' and
                            rgn.type == TYPE_SUBREAD):
                        observed.append(rgn.advance)
                log.debug("Observed: %s" % ",".join(map(str, observed)))
                if len(observed) > 0:
                    uncensored.append(int(np.mean(observed)))
                    log.debug("Uncensored: %s", str(int(np.mean(observed))))
            elif len(pattern) >= 5:
                observed = []
                for i, rgn in enumerate(rgns):
                    if (rgn.type == TYPE_SUBREAD and i - 2 >= 0 and
                            pattern[i - 2] == (1 - pattern[i]) and
                            (i + 2) < len(pattern) and
                            pattern[i + 2] == (1 - pattern[i])):
                        observed.append(rgn.advance)
                log.debug("Observed: %s", ",".join(map(str, observed)))
                if len(observed) > 0:
                    uncensored.append(int(np.mean(observed)))
                    log.debug("Uncensored: %s", str(int(np.mean(observed))))

        self._censored = np.array(censored)
        self._uncensored = np.array(uncensored)
        log.info("Censored: %.2f, %.2f, %d" % (
            np.mean(self._censored), np.std(self._censored),
            len(self._censored)))
        log.info("Uncensored: %.2f, %.2f, %d" % (
            np.mean(self._uncensored), np.std(self._uncensored),
            len(self._uncensored)))
        if len(self._uncensored) == 0 or len(self._censored) == 0:
            self._predicted = self._uncensored
            return

        if self.link == 'exp':
            predicted = (
                1000.0 * np.log(self._gibbs_sample(
                    np.exp(self._censored / 1000.0),
                    np.exp(self._uncensored / 1000.0),
                    int(self.iterations))))
        elif self.link == 'log':
            predicted = np.exp(self._gibbs_sample(
                np.log(self._censored), np.log(self._uncensored),
                int(self.iterations)))
        elif self.link == 'identity':
            predicted = self._gibbs_sample(self._censored, self._uncensored,
                                           int(self.iterations))
        elif self.link == 'square':
            predicted = np.sqrt(1000.0 * self._gibbs_sample(
                np.power(self._censored, 2.0) / 1000.0,
                np.power(self._uncensored, 2.0) / 1000.0,
                int(self.iterations)))
        else:
            raise SystemExit("Don't know link function %s" % self.link)

        self._predicted = np.concatenate((predicted, uncensored))
        log.info("Predicted: %.2f, %.2f, %d" % (
            np.mean(self._predicted), np.std(self._predicted),
            len(self._predicted)))

    def _gibbs_sample(self, censored, uncensored, n):
        """Taken almost directly from JMS's R code"""

        n_samples = len(censored) + len(uncensored)
        n_censored = len(censored)

        mu = np.empty(shape=(n, ), dtype=np.float32)
        mu.fill(0)
        sg = np.empty(shape=(n, ), dtype=np.float32)
        sg.fill(1)
        xc = np.empty(shape=(n_censored, n), dtype=np.float32)
        xc.fill(0)

        mu[0] = np.mean(np.concatenate((censored, uncensored)))
        sg[0] = np.std(np.concatenate((censored, uncensored)))
        xc[:, 0] = censored

        for i in range(1, n):
            xb = np.mean(np.concatenate((uncensored, xc[:, i - 1])))

            # sample from p( mu | ... ) - t distribution
            mu[i] = xb + sg[i - 1] / np.sqrt(n_samples) * t.rvs(n_samples - 1)

            # sample from p( sigma | ... ) - inverse-chi2 distribution
            ss = (np.sum((xc[:, i - 1] - mu[i]) ** 2) +
                  np.sum((uncensored - mu[i]) ** 2))
            sg[i] = np.sqrt(ss / chi2.rvs(n_samples - 1))

            # sample from p( x | ... )
            for j in range(n_censored):
                #p = uniform.rvs()
                #z = ( censored[j] - mu[i] ) / sg[i]
                #p = p + ( 1 - p ) * norm.pdf( z )
                #xc[j,i] = norm.ppf( p, mu[i], sg[i] )
                try:
                    xc[j, i] = truncnorm.rvs((censored[j] - mu[i]) / sg[i],
                                             np.inf, loc=mu[i], scale=sg[i])
                    if xc[j, i] > 99999:
                        raise ValueError
                except ValueError:  # Precision problems
                    log.debug("Failed to compute TN for %.2f (%.2fm, %.2fs)" %
                              (censored[j], mu[i], sg[i]))
                    #import pdb; pdb.set_trace()
                    xc[j, i] = xc[j, i - 1]

        return xc[:, -1].copy()

    def _plot_insert_distributions(self):

        INSERT_ALPHA = 0.5
        ax_rect = [INSERTS_LEFT_R, INSERTS_BOTTOM_R, INSERTS_WIDTH_R,
                   INSERTS_HEIGHT_R]
        ax = self.fig.add_axes(ax_rect, frame_on=True)
        #, title="Insert Length Distribution" )

        legend_rects, labels = [], []
        for data, color, label in [
                (self._censored, "#990000", "Censored"),
                (self._uncensored, "#000000", "Uncensored"),
                (self._predicted, "#CCCCCC", "Predicted")]:

            if len(data) == 0:
                log.warning("Zero-length data for %s" % label)
                continue
            c, b = np.histogram(data, bins=30)
            x = b[:-1]
            y = np.zeros(len(x))

            #ax2 = ax.twinx()
            ax.fill_between(x, y, c, alpha=INSERT_ALPHA,
                            edgecolor="#000000", facecolor=color)

            legend_rects.append(matplotlib.patches.Rectangle(
                (0.1, 0.1), 0.1, 0.1, facecolor=color, edgecolor="#000000",
                alpha=INSERT_ALPHA))
            ## XXX: JHB Martin
            try:
                if np.isinf(max(data)):
                    labels.append("%s [inf]")
                else:
                    labels.append("%s [%dbp]" % (label, int(np.mean(data))))
            except:
                log.warn("Incapable of properly setting lables due to "
                         "infintes.")
                pass

        ax.set_xlabel("Insert Length", size=X_LABEL_SIZE)
        ax.set_ylabel("# Reads", size=Y_LABEL_SIZE)

        try:
            ax.legend(legend_rects, labels)
            ax.get_legend().get_frame().set_edgecolor('#a0a0a0')
        except ZeroDivisionError:
            # Obscure bug here -don't generate a legend if we hit it.
            pass

    def _make_length_histograms(self, control_data, sample_data, outfile):
        """Predicted length histograms"""
        len_max = max(max(control_data[1, :]), max(sample_data[1, :]))
        pass

    def _generate_zmw_map(self):
        end = 0
        for i, (zmw, hq) in enumerate(self._hq.iteritems()):
            if i >= self.max_num_zmw:
                continue
            end = max(len(hq), end)
        #MAX_ZMW_SIZE = min(int(math.ceil(float(end) / 1000)) * 1000, 20000)
        self._max_zmw_size = min(
            int(math.ceil(float(end) / 1000)) * 1000, 20000)

        #shape = (min(len(self._hq), self.max_num_zmw), MAX_ZMW_SIZE)
        shape = (min(len(self._hq), self.max_num_zmw), self._max_zmw_size)

        log.info("Initializing map to size %s" % str(shape))
        self._map = np.empty(shape=shape, dtype=np.uint8)
        self._annotations = np.recarray(
            shape=shape[0],
            dtype=[('AdapterStart', np.uint16), ('Pattern', '|S30'),
                   ('SubreadLen', np.uint16), ('HQLen', np.uint16),
                   ('NAdapters', np.uint8), ('NSubreads', np.uint8)])

        self._map.fill(TYPE_EMPTY)

        for i, (zmw, hq) in enumerate(self._hq.iteritems()):
            if i >= self.max_num_zmw:
                continue
            # This doesn't actually do anything atm since we only right-trim.
            #if zmw in self._regions:
            #    map( lambda r: r.addDelta( -hq.start ), self._regions[ zmw ] )

            pattern = ''
            adapter_centroid = [0, 0]
            subread_len = 0
            n_adapters = 0
            n_subreads = 0

            #end = min(len(hq), MAX_ZMW_SIZE)
            end = min(len(hq), self._max_zmw_size)
            # All reads in subreads.bam will be high quality, this tag is
            # meaningless...
            self._map[i][0:end] = TYPE_HQ
            if zmw in self._regions:
                my_regions = self._regions[zmw]
                my_regions.sort(key=lambda r: r.start)
                for pbrange in my_regions:
                    #if pbrange.start < MAX_ZMW_SIZE:
                    if pbrange.start < self._max_zmw_size:
                        #self._map[i][pbrange.start:
                             #min(pbrange.end, MAX_ZMW_SIZE)] = pbrange.type
                        self._map[i][pbrange.start:min(
                            pbrange.end, self._max_zmw_size)] = pbrange.type
                    # We're not really getting adapter ranges, so this number
                    # is meaningless. Rather we should store and pull adapter
                    # number from the pbrange associated with subreads
                    # (somehow). Also, adapter_centroid is never used, so the
                    # bounaries of the adapter aren't really valuable.
                    if pbrange.type == TYPE_ADAPTER:
                        adapter_centroid[0] += np.mean(
                            [pbrange.start, pbrange.end])
                        adapter_centroid[1] += 1
                        n_adapters += 1
                        pattern += "A"
                    if pbrange.type == TYPE_SUBREAD:
                        subread_len += len(pbrange)
                        n_subreads += 1
                        pattern += "S"

            self._annotations[i]['AdapterStart'] = (
                0 if adapter_centroid[1] == 0 else (
                    adapter_centroid[0] / adapter_centroid[1]))
            self._annotations[i]['Pattern'] = pattern
            self._annotations[i]['HQLen'] = len(hq)
            self._annotations[i]['SubreadLen'] = subread_len
            self._annotations[i]['NAdapters'] = n_adapters
            self._annotations[i]['NSubreads'] = n_subreads

        # Free up some memory
        del self._regions
        del self._hq

    def _sort_zmw_map(self):
        self._sort()

    def _sort(self):
        log.info("Sorting and Reordering Map...")
        order = np.argsort(self._annotations,
                           order=['NSubreads', 'NAdapters', 'Pattern',
                                  'AdapterStart', 'SubreadLen', 'HQLen'])
        log.info("Sorting Complete.")

        self._map = self._map[order]
        self._annotations = self._annotations[order]
        log.info("Reordering Complete.")

    def _plot_zmw_map(self):

        nsr = self._annotations['NSubreads']
        nad = self._annotations['NAdapters']
        subsets = [("0 Subreads", (nsr == 0)),
                   (">0 Subreads, 0 Adapters", (nsr > 0) & (nad == 0)),
                   (">0 Subreads, 1 Adapter", (nsr > 0) & (nad == 1)),
                   (">0 Subreads, >1 Adapter", (nsr > 0) & (nad > 1))]

        subset_sizes = [len(np.nonzero(subset)[0])
                        for title, subset in subsets]

        column_pixels = FIGURE_HEIGHT_R * PLOT_HEIGHT
        zmws_per_pixel = max(subset_sizes) / column_pixels
        zmws_per_pixel = max(1.0, zmws_per_pixel)

        row_pixels = COLUMN_WIDTH_R * PLOT_WIDTH
        #MAX_ZMW_SIZE = 20000
        self._max_zmw_size = 20000
        #basesPerPixel = MAX_ZMW_SIZE / rowPixels
        bases_per_pixel = self._max_zmw_size / row_pixels

        for i, (title, condition) in enumerate(subsets):

            left = (FIGURE_LEFT_R + COLUMN_SPACE_R * i * 0.5 +
                    COLUMN_SPACE_R * (i + 1) * 0.5 + COLUMN_WIDTH_R * i)
            width = COLUMN_WIDTH_R

            height = ((1 / zmws_per_pixel) * subset_sizes[i]) / PLOT_HEIGHT
            bottom = 1.0 - FIGURE_TOP_R - height

            ax_rect = [left, bottom, width, height]
            ax = self.fig.add_axes(ax_rect, frame_on=True,
                                   title=title.replace(",", ",\n"))

            # Axis Beautification
            #ax.set_xticks(range(0, MAX_ZMW_SIZE + 1,
            #int(math.floor(MAX_ZMW_SIZE / NUM_X_TICKS))))
            ax.set_xticks(range(0, self._max_zmw_size + 1,
                                int(math.floor(
                                    self._max_zmw_size / NUM_X_TICKS))))
            map(lambda x: x.set_rotation(X_TICK_ROTATION),
                ax.get_xticklabels())
            map(lambda x: x.set_size(X_TICK_SIZE), ax.get_xticklabels())

            ax.set_yticks([])

            log.info("Converting Image for '%s' to Figure" % title)
            if subset_sizes[i] > 0:
                ax.imshow(self._map[condition], aspect='auto', cmap=CMAP,
                          vmin=COLOR_MIN, vmax=COLOR_MAX,
                          interpolation='bicubic')

        global_axis = self.fig.add_axes([FIGURE_LEFT_R, FIGURE_BOTTOM_R,
                                         FIGURE_WIDTH_R, FIGURE_HEIGHT_R],
                                        frame_on=False)
        global_axis.axvline(0, 0, 1, color="black")
        global_axis.set_xticks([])
        max_y = int(max(subset_sizes) * self._downsample_rate)
        yti = Y_TICK_INC
        while max_y / yti > 50:
            yti *= 2
        yticks = [float(x) for x in range(0, max_y, yti)] + [max_y]
        global_axis.set_yticks(yticks)
        global_axis.set_yticklabels(["%.0f" % tick for tick in yticks],
                                    size=Y_TICK_SIZE)
        global_axis.set_ylabel("# Reads", size=Y_LABEL_SIZE)
        self.fig.text(FIGURE_LEFT_R + FIGURE_WIDTH_R / 3, FIGURE_X_TITLE,
                      "Read Position (bp)", size=Y_LABEL_SIZE)
        global_axis.invert_yaxis()

    def _plot_legend(self):
        legend_data = np.empty(shape=(LEGEND_HEIGHT, LEGEND_WIDTH))
        unit_height = LEGEND_HEIGHT / len(LEGEND_LABELS)
        for y, label in enumerate(LEGEND_LABELS):
            legend_data[y * unit_height: (y + 1) * unit_height,
                        0: LEGEND_WIDTH] = LEGEND_TYPES[y]
            legend_data[y * unit_height, 0: LEGEND_WIDTH] = COLOR_BLACK
            self.fig.text(
                (LEGEND_X + LEGEND_WIDTH + LEGEND_X_SPACING) / PLOT_WIDTH,
                (LEGEND_Y + LEGEND_Y_OFFSET + y * unit_height) / PLOT_HEIGHT,
                LEGEND_LABELS[-(y + 1)], size=LGD_LABEL_SIZE)
        legend_data[:, 0] = COLOR_BLACK
        legend_data[:, -1] = COLOR_BLACK
        legend_data[0, :] = COLOR_BLACK
        legend_data[-1, :] = COLOR_BLACK

        self.fig.figimage(legend_data, xo=LEGEND_X, yo=LEGEND_Y, cmap=CMAP,
                          vmin=COLOR_MIN, vmax=COLOR_MAX)

    def _plot_loading_bars(self):

        nsr = self._annotations['NSubreads']
        nad = self._annotations['NAdapters']

        def count(condition):
            return len(np.nonzero(condition)[0])

        subset_names = ["HQ, With\nAlignments", "HQ, No\nAlignments",
                        "No HQ\nSequence"]
        subset_values = [count(nsr > 0) * self._downsample_rate,
                         count(nsr == 0) * self._downsample_rate,
                         self._empty_zmws]
        log.debug("Located values for subset lengths: %s" % str(subset_values))
        total_zmws = sum(subset_values)
        subset_values = [100.0 * sv / total_zmws for sv in subset_values]

        bar_rect = [BAR_LEFT_R, BAR_BOTTOM_R, BAR_WIDTH_R, BAR_HEIGHT_R]
        ax = self.fig.add_axes(bar_rect, frame_on=False)

        posns = np.arange(len(subset_names)) + 0.5
        bars = ax.barh(posns, subset_values, align='center', height=0.5,
                       color=BAR_COLOR)

        for bar in bars:
            bar_width = bar.get_width()
            pct = "%.1f%%" % bar_width
            xloc = bar_width - 1
            clr = "white"
            align = "right"
            if bar_width < MIN_BAR_WIDTH:  # Too small to fit
                xloc = bar_width + 1
                clr = "black"
                align = "left"
            yloc = bar.get_y() + bar.get_height() / 2
            ax.text(xloc, yloc, pct, horizontalalignment=align,
                    verticalalignment='center', color=clr, weight='bold',
                    size=PCT_LABEL_SIZE)

        ax.axis([0, max(subset_values) + 1, 0, len(subset_names)])
        ax.set_yticks(posns)
        ax.set_yticklabels(subset_names, size=BAR_LABEL_SIZE)

        ax.set_xlabel("% of Reads", size=X_LABEL_SIZE)
        ax.set_xticks([])

        ax.axvline(0, 0, 1, color="black")

    def _plot_help(self):

        lines = [
            "The figure on the left depicts the analysis results of all Reads",
            "which passed filtering and have at least one high quality base.",
            "Reads are ordered by (in decreasing importance) #subreads,",
            "#adapters, subread/adapter pattern, average adapter start "
            "position,",
            "subread length, # HQ bases.",
            "",
            "The graph in the upper right shows the % of Reads which fall "
            "into three",
            "major categories:",
            "No HQ Sequence - Reads which have no good pulses or were "
            "filtered out",
            "HQ, No Alignments - Some high quality sequence but no alignments",
            "HQ, With Alignments - Reads that contain alignable sequence"]

        box = np.empty(
            shape=(HELP_TEXT_HEIGHT * len(lines) + HELP_TEXT_PADDING * 2,
                   HELP_TEXT_WIDTH + HELP_TEXT_PADDING * 2), dtype=np.uint8)
        box.fill(HELP_TEXT_COLOR)
        self.fig.figimage(box, xo=HELP_TEXT_X - HELP_TEXT_PADDING,
                          yo=HELP_TEXT_Y - HELP_TEXT_PADDING, cmap=CMAP,
                          vmin=COLOR_MIN, vmax=COLOR_MAX)
        for i in range(len(lines)):
            self.fig.text(HELP_TEXT_X / PLOT_WIDTH,
                          (HELP_TEXT_Y + i * HELP_TEXT_HEIGHT) / PLOT_HEIGHT,
                          lines[-(i + 1)], size=HELP_TEXT_SIZE)

    def _write_out_figure(self):
        # Free up some memory before we make this last call, which takes a lot.
        del self._map
        del self._annotations
        time.sleep(2)
        log.info("Writing Image to %s" % self.outfn)
        plt.savefig(self.outfn, format="png", dpi=DPI)
        log.info("Writing Complete.")


def make_report(infn, outfn, iterations=100, max_num_zmw=MAX_NUM_ZMW,
                resample=False, link=None):
    app = ReadMapReport(infn, outfn, iterations=iterations,
                        max_num_zmw=max_num_zmw, resample=resample, link=link)
    app.run()
    return 0


def arg_runner(args):
    """Callback for pbsystem cmdline functions"""
    log.info(
        "Starting {f} version {v} report generation".format(f=__file__,
                                                            v=__version__))
    app = ReadMapReport(args.infn, args.outfn,
                        iterations=args.iterations,
                        max_num_zmw=args.max_num_zmw, resample=args.resample,
                        link=args.link)
    app.run()
    return 0


def _get_parser():
    """Handle command line argument parsing"""

    description = ReadMapReport.__doc__
    parser = get_default_argparser_with_base_opts(
        version=__version__,
        description=description)
    parser.add_argument("infn", type=validate_file,
                        help="The input bam file")
    parser.add_argument("--outfn", type=str, default="readmap.png",
                        help="The name of the output png file")
    parser.add_argument("-n", "--iterations", default=100, type=int,
                        help=("The number of iterations to perform when using "
                              "gibbs sampling to estimate the insert length "
                              "distribution"))
    parser.add_argument("--max_num_zmw", type=int, default=MAX_NUM_ZMW,
                        dest='max_num_zmw',
                        help="Maximum number of ZMWs to use for downsampling")
    parser.add_argument("--resample", action="store_true", default=False,
                        help=("Resample the distribution before Gibbs "
                              "sampling unbias the read distribution"))
    parser.add_argument("--link", default='exp', type=str,
                        help=("Specify a link function [exp,identity,log] "
                              "(default=exp)"))
    parser.set_defaults(func=arg_runner)
    return parser


def main(argv=sys.argv):
    """Main point of entry"""
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=arg_runner,
        alog=log,
        setup_log_func=setup_log)
