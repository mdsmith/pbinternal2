#!/usr/bin/env python
"""
Script to generate a json report of per-reference rainbow plots and summary
metrics.
Example:
    $ makeRainbowReport.py --output=resultDir/ infile.bam
"""

import sys
import os
import logging
import numpy as np
from pbcore.io import openIndexedAlignmentFile
from pbcommand.models.report import (Report, Table, Column, PlotGroup, Plot,
                                     Column)
from pbcommand.utils import setup_log
from pbcommand.cli import pacbio_args_runner, get_default_argparser
from pbreports.io.validators import (validate_file,
                                     validate_output_dir)
from pbinternal2 import get_version
from pbinternal2.util import alignment_to_png
from pbcore.io import DataSet

log = logging.getLogger(__name__)
__version__ = get_version()


def _data_by_reference(in_fn):
    """Returns a dictionary of numpy arrays, keyed on reference"""
    data = {}
    names = DataSet(in_fn).externalResources.resourceIds
    for name in names:
        # Name here is a URI. There are a lot of ways to solve this problem. As
        # an internal report, this is fine:
        with openIndexedAlignmentFile(name.split(':')[-1]) as in_file:
            for row in in_file:
                #try as a bam, fall back to cmph5
                try:
                    reference = row.referenceName
                    n_ins = row.aEnd - row.aStart - row.nM - row.nMM
                    n_del = row.tEnd - row.tStart - row.nM - row.nMM
                    read_length = row.aEnd - row.aStart
                    accuracy = (1.0 - (row.nMM + n_ins + n_del) /
                                float(read_length))
                except (AttributeError, IndexError):
                    reference = in_file.referenceInfo(row.RefGroupID).FullName
                    read_length = row.rEnd - row.rStart
                    accuracy = (1.0 - (row.nMM + row.nIns + row.nDel) /
                                float(read_length))
                data.setdefault(reference, [])
                data[reference].append((read_length, accuracy))

    for reference in data.keys():
        data[reference] = np.array(data[reference],
                                   dtype=[("Length", int),
                                          ("Accuracy", float)])
    return data


def make_report(in_fn, out_dir='.'):
    """Make a rainbow report from the bam, cmp.h5 or DataSet XML file

    Args:
        in_fn: input sequences in bam, cmp.h5 or DataSet xml format
        out_dir: the location in which to safe the figures
    Returns:
        A pbreports Report object with the table and plot
    Emissions:
        out_dir/rainbow{i}.png one rainbow plot for each reference
    """
    data = _data_by_reference(in_fn)
    if len(data.keys()) == 0:
        log.error("No results found in input file, aborting.")
        return 1

    report = Report('rainbow_report')
    columns = (Column('reference_name_column',
                      header="Reference"),
               Column('aligned_subreads_column',
                      header="Aligned Subreads"),
               Column('mean_accuracy_column',
                      header="Mean Accuracy"),
               Column('mean_subread_length_column',
                      header='Mean Subread Readlength'))
    table = Table('stats_table',
                  title='Per-Reference Summary Stats',
                  columns=columns)
    report.add_table(table)

    for i, reference in enumerate(data.keys()):
        columns[0].values.append(reference[:20])
        columns[1].values.append("%d" % len(data[reference]["Accuracy"]))
        columns[2].values.append('%.2f%%' % (
            100.0 * np.mean(data[reference]["Accuracy"])))
        columns[3].values.append('%.0f bp' % (
            np.mean(data[reference]["Length"])))

        png_fn = "rainbow%d.png" % i
        alignment_to_png.make_report(in_fn, out_dir=out_dir, name=png_fn,
                                     dpi=100, reference=reference)
        plot = Plot('reference_{i}_plot'.format(i=i), png_fn)
        plot_group = PlotGroup('reference_{i}_plot_group'.format(i=i),
                               title='Acc vs. RL for {s}'.format(s=reference))
        plot_group.add_plot(plot)
        report.add_plotgroup(plot_group)
    return report


def arg_runner(args):
    """Callback for pbsystem cmdline functions"""
    report = make_report(args.infile, out_dir=args.out_dir)
    if report == 1:
        return 1
    report.write_json(os.path.join(args.out_dir, 'rainbow.json'))
    return 0


def _get_parser():
    """Handle command line argument parsing"""

    description = __doc__
    #parser = argparse.ArgumentParser(version=__version__,
    #                                 description=description)
    p = get_default_argparser(__version__, description)

    p.add_argument("infile", type=validate_file,
                   help=('The input sequence file (bam, cmp.h5 or '
                         'DataSet XML)'))
    p.add_argument("--output", default='.', type=validate_output_dir,
                   dest='out_dir',
                   help="The report output directory")
    p.set_defaults(func=arg_runner)
    return p


def main(argv=sys.argv):
    """Main point of entry"""
    log.info("Starting {f} version {v} report generation".format(
        f=__file__, v=__version__))
#    return main_runner_default(argv[1:], _get_parser(), log)
    return pacbio_args_runner(argv[1:], _get_parser(), arg_runner, log, setup_log)
