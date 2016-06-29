#!/usr/bin/env python

"""Missing Adapter Report

Process those reads missing adapters

"""
import argparse
import logging
import os
import sys
import numpy as np

from pbcommand.models.report import Report, Attribute, PlotGroup, Plot
from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcore.io import FastaReader

from pbreports.io.validators import validate_file, validate_output_dir
from pbreports.plot.helper import (get_fig_axes_lpr, apply_histogram_data,
                                   save_figure_with_thumbnail)

from pbinternal2 import get_version
from pbinternal2.util.fastaShim import ensure_fasta, ensure_fastas

log = logging.getLogger(__name__)
__version__ = get_version()


def make_report(read_fn, missed_adapter_fns, out_dir):
    """Generate the Report

    Make the PlotGroup, populate with Plots, add to the report.

    Args:
        read_fn: The general population of reads in fasta, bam or DataSet XML
        missed_adapter_fns: a list of one or more files (fasta, bam, or DataSet
                            XML) of reads missing adapters
        output: an optional specified output file (default='.',
                ignored and output to stdout)
    """
    length_thresholds = [0, 500, 1000, 2000]

    # Build the histogram of read lengths for the general population:
    read_fn = ensure_fasta(read_fn)
    initial_lengths = _get_lengths(read_fn)
    plot_group = PlotGroup('subread_plot_group', title=('All subreads'))
    png_file_name = os.path.join(out_dir, 'read_lengths.png')
    _make_histogram_png_theme_lpr(initial_lengths, png_file_name,
                                  ('Subread length', 'Count'))
    plot = Plot('subread_plot', image=os.path.basename(png_file_name))
    plot_group.add_plot(plot)
    report = Report('missing_adapter_report', plotgroups=[plot_group])

    # Build histograms for each of the files of reads with missed adapters:
    missed_adapter_fns = ensure_fastas(missed_adapter_fns)
    for i, fasta in enumerate(missed_adapter_fns):
        plot_group = PlotGroup(
            'fasta_plot_group',
            title=('Subread length distributions for adapter phenotypes'))
        prefix = os.path.basename(fasta).split(".")[0].lower()
        png_file_name = os.path.join(out_dir, '%s.png' % prefix)
        png_thumb_name = os.path.join(out_dir, '%s_thumb.png' % prefix)
        missed_adapter_lengths = _get_lengths(fasta)
        _make_histogram_png_theme_lpr(missed_adapter_lengths, png_file_name,
                                      ('Read length', 'Count'))
        plot = Plot("{s}_{i}_plot".format(s=prefix, i=i),
                    os.path.basename(png_file_name),
                    thumbnail=os.path.basename(png_thumb_name),
                    caption='%s adapter subreads' % prefix)

        if len(missed_adapter_lengths):
            for min_length in length_thresholds:
                num_missed = _num_gt_min(missed_adapter_lengths, min_length)
                num_all = _num_gt_min(initial_lengths, min_length)
                attr_id = '{p}_{m}_plot_attr'.format(p=prefix, m=min_length)
                if num_all == 0:
                    continue
                attribute = Attribute(
                    attr_id,
                    ('%.2f = %i / %i' % (
                        (100.0 * float(num_missed) / num_all),
                        num_missed,
                        num_all)),
                    name='%s adapter (length > %i):' % (prefix, min_length))
                report.add_attribute(attribute)
        plot_group.add_plot(plot)
        report.add_plotgroup(plot_group)
    return report


def _num_gt_min(values, minimum):
    return float(len([x for x in values if x > minimum]))


def _get_lengths(fasta):
    return np.array([len(e.sequence) for e in FastaReader(fasta)])


def _make_histogram_png_theme_lpr(initial_lengths, png_file_name, ax_names,
                                  bins=30):
    fig, axes = get_fig_axes_lpr()
    apply_histogram_data(axes, initial_lengths, bins, ax_names)
    save_figure_with_thumbnail(fig, png_file_name)


def arg_runner(args):
    """Callback for pbsystem cmdline functions"""
    log.info(
        "Starting {f} version {v} report generation".format(f=__file__,
                                                            v=__version__))
    report = make_report(args.read_fn, args.missed_adapter_fns,
                         args.out_dir)
    report.write_json(os.path.join(args.out_dir,
                                   'missing_adapter.json'))
    return 0


def _get_parser():
    """Returns an argparse parser for the positional and optional inputs"""
    parser = get_default_argparser_with_base_opts(
        version=__version__,
        description=__doc__)
    parser.add_argument('read_fn',
                        help=('Name of file (fasta, bam, xml) containing '
                              'subreads'),
                        type=validate_file)
    parser.add_argument('--output', dest='out_dir',
                        help=('Output directory for report json and '
                              'associated files'),
                        type=validate_output_dir,
                        default='.')
    parser.add_argument('missed_adapter_fns',
                        help=('Name of file(s) (fasta, bam, xml) containing '
                              'reads missing adapters'),
                        type=validate_file, nargs='+')
    return parser


def main(argv=sys.argv):
    """Main point of entry"""
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=arg_runner,
        alog=log,
        setup_log_func=setup_log)
