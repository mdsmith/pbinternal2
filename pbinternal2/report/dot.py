#!/usr/bin/env python

""" Dot Report

Make a dot report from input reads and reference fasta files.

"""
import argparse
import logging
import os
import sys
import tempfile
import shutil

from pbcommand.cli import (pacbio_args_runner, get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcommand.models.report import Report, PlotGroup, Plot
from pbcore.io import FastaReader
from pbinternal2 import get_version
from pbinternal2.external_tools import gepard
from pbinternal2.util.fastaShim import ensure_fasta
from pbreports.io.validators import (validate_file,
                                     validate_output_dir)

log = logging.getLogger(__name__)
__version__ = get_version()

TEMP_DIR = tempfile.mkdtemp(suffix="dot_report")


def _add_dot_plot(plot_group, entry, idx, out_dir, max_width, word_length,
                  ref_fn):
    """Use gepard to generate dot plots

    Args:
        plot_group: The plot group to which dot plot Plot objects are added
        entry: The FastaReader entry in question
        idx: The number of the entry from an enumerated list of entries
        out_dir: The output directory for dot plot .pngs
        max_width: The maximum width of a dot plot
        word_length: Substring length for alignment
        ref_fn: The file name of the reference Fasta

    Emits:
        temporary fastas, one per entry, as gepard input
        dot plot pngs, one per entry (gepard output)

    Affects:
        plot_group: plots are added to plot group, one per read.
    """
    name = entry.name
    prefix = (entry.name.replace(os.sep, "_").replace("|", "_").split(
        " ")[0])
    fasta_temp = "%s/%s.fasta" % (TEMP_DIR, prefix)
    fasta_fh = open(fasta_temp, "w")
    fasta_fh.write(str(entry))
    fasta_fh.close()

    png_fn = os.path.join(out_dir, '%s.png' % prefix)
    gepard(max_width, word_length, ref_fn, fasta_temp, png_fn)

    plot = Plot(
        "dotplot_" + str(idx),
        image=os.path.basename(png_fn),
        caption=("Dot plot for %s (contig %i)" % (name.split(" ")[0],
                                                  idx)))

    plot_group.add_plot(plot)
    os.unlink(fasta_temp)


def make_report(read_fn, ref_fn, out_dir, min_length=500, word_length=10,
                debug=False):
    """Creates a report using gepard to create a dot plot comparing a
    sequence (usually from de novo assembly) to a reference.

    Args:
        read_fn: The name of a fasta, bam, or xml file containing reads
        ref_fn: The name of a fasta or xml file containing the reference (only
                the first sequence in this file is used)
        out_dir: The output directory for the report json and affiliated image
                 files.
        min_length: The minimum length of a read. Will be at least 100
                    for gepard's benefit (500).
        word_length: The length of substrings aligned by gepard (10).
        debug: If set, temp files are retained.
    """
    log.debug("Temporary files can be found here: {s}".format(s=TEMP_DIR))
    absolute_min_dot_length = 100
    reference_to_plot_width_ratio = 500
    minimum_plot_width = 750
    if int(min_length < absolute_min_dot_length):
        print >> sys.stderr, ("Setting minDotLength to %i because %i was "
                              "too low, can cause gepard dotplotter to "
                              "error out" % (absolute_min_dot_length,
                                             min_length))
        min_length = absolute_min_dot_length

    # The reference fasta is necessary for gepard input, therefore just convert
    # early and use the same old code paths.
    ref_fn = ensure_fasta(ref_fn)
    reader = FastaReader(ref_fn)
    for entry in reader:
        ref_length = len(entry.sequence)
        ref_name = entry.name.split(" ")[0]
        break
    max_width = int(ref_length / reference_to_plot_width_ratio)
    max_width = max(minimum_plot_width, max_width)

    plot_group = PlotGroup('dot_plot_group',
                           title=('Per read dot plots for %s' % ref_name))
    read_fn = ensure_fasta(read_fn)
    reader = FastaReader(read_fn)
    for idx, entry in enumerate(reader):
        if len(entry.sequence) < int(min_length):
            continue
        _add_dot_plot(plot_group, entry, idx + 1, out_dir, max_width,
                      word_length, ref_fn)

    report = Report('dot_report', plotgroups=[plot_group])
    if not debug:
        shutil.rmtree(TEMP_DIR)
    return report


def arg_runner(args):
    log.info("Starting {f} version {v} report generation".format(
        f=__file__, v=__version__))
    report = make_report(args.readFile, args.refFile,
                         args.output, args.minLength, args.wordLength,
                         args.debug)
    report.write_json(os.path.join(args.output, args.jsonOut))
    return 0


def get_parser():
    description = __doc__
    parser = argparse.ArgumentParser(version=__version__,
                                     description=description)
    parser.add_argument('readFile',
                        help=('Input file (fasta, bam or DataSet xml) '
                              'containing the reads'),
                        type=validate_file)
    parser.add_argument('refFile',
                        help=('Input file (fasta, or ReferenceSet xml) '
                              'containing the reference (only the first will '
                              'be used)'),
                        type=validate_file)
    parser.add_argument('--output', type=validate_output_dir,
                        help='Output directory for associated files',
                        default='.')
    parser.add_argument('--minLength', help='Minimum length for dot plot',
                        type=int, default=500)
    parser.add_argument('--wordLength', help='Word length for dot plot',
                        type=int, default=10)
    parser.add_argument('--jsonOut', help='json output file', type=str,
                        default="dot_report.json")
    parser.add_argument('--debug', help='Save intermediate files',
                        default=False, action='store_true')
    parser.set_defaults(func=arg_runner)
    return parser


def main(argv=sys.argv):
    """Main point of Entry"""
    return pacbio_args_runner(
        argv=argv[1:],
        parser=get_parser(),
        args_runner_func=arg_runner,
        alog=log,
        setup_log_func=setup_log)
