#!/usr/bin/env python

"""AlignmentToPng Report

Convert an input bam or DataSet XML file (used to be a cmp.h5 file, thus the
name) to a figure of Accuracy vs. Subread Length. Each point on the graph
represents the accuracy and length of a single subread as measured by a local
alignment to the reference. The points are colored by qv-score.
"""

import sys
import os
import argparse
import logging

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbcommand.validators import (validate_file, validate_output_dir)
from pbreports.plot.rainbow import make_report
from pbcore.io import openIndexedAlignmentFile, DataSet
from pbcommand.models.report import Report, PlotGroup, Plot


from pbinternal2 import get_version

log = logging.getLogger(__name__)
__version__ = get_version()


def arg_runner(args):
    """Callback for pbsystem cmdline functions"""
    log.info(
        "Starting {f} version {v} report generation".format(f=__file__,
                                                            v=__version__))
    report = make_report(args.infile, out_dir=args.out_dir, bounds=args.bounds,
                         nolegend=args.nolegend, reference=args.reference,
                         dpi=args.dpi)
    stem = os.path.splitext(os.path.basename(args.infile))[0]
    outjson = os.path.join(args.out_dir, '%s.json' % (stem))
    report.write_json(outjson)
    return 0


def _get_parser():
    description = __doc__
    parser = get_default_argparser_with_base_opts(
        version=__version__,
        description=description)
    parser.add_argument("infile", type=validate_file,
                        help="Input alignment file (bam or DataSet XML)")
    parser.add_argument("--output", type=validate_output_dir, default='.',
                        dest='out_dir',
                        help="Output directory for associated files")
    parser.add_argument("--json", action="store_true", default=False,
                        help="Output an JSON description of the resulting "
                             "graph")
    parser.add_argument("--bounds", type=str, default=None,
                        help="The xmin:xmax:ymin:ymax for the plot")
    parser.add_argument("--nolegend", action="store_true", default=False,
                        help="Don't put a legend on the plot")
    parser.add_argument("--reference", type=str, default=None,
                        help="Specify the name of a reference to plot values "
                             "for that reference only.")
    parser.add_argument("--dpi", type=int, default=60,
                        help="Specify the DPI of the resulting image.")
    return parser


def main(argv=sys.argv):
    """Main point of entry"""
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=arg_runner,
        alog=log,
        setup_log_func=setup_log)
