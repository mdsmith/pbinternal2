#!/usr/bin/env python
"""
Given a fasta file and a reference, align the subreads in the fasta to
the reference and determine which appear to have issues with adapters.

Usage: findAdapterIssues [--help] reads.fsta reference.fsta
"""

import sys
import logging
import argparse
import tempfile
import os
import shutil

#from pbinternal2.tools import blasr
from pbinternal2.util.bamFixer import rectify
from pbinternal2 import get_version
from pbinternal2.util.Range import Range
from pbalign.pbalignrunner import PBAlignRunner
from pbcore.io import openAlignmentFile
from pbcommand.cli import (pacbio_args_runner,
                           get_default_argparser_with_base_opts)
from pbcommand.utils import setup_log
from pbreports.io.validators import (validate_file,
                                     validate_output_dir)

log = logging.getLogger(__name__)
if not log.handlers:
    log.addHandler(logging.StreamHandler())
__version__ = get_version()

ADAPTER_MIN_LENGTH = 30
CHIMERA_QUERY_MAX_OVERLAP = 30
CHIMERA_TARGET_DISTANCE = 1000
ANNOTATE_NON_ISSUES = True
CHIMERA_END_TOLERANCE = 100
ADAPTER_MAX_LENGTH = 60
OVERHANG_MIN_DISTANCE = 50


def _get_longest_on_strand(hits, strand="+"):
    hits.sort(key=len, reverse=True)
    strand_hits = [hit for hit in hits
                   if (hit.isForwardStrand and (strand == "+"))
                   or (hit.isReverseStrand and (strand == "-"))]
    return strand_hits[0] if strand_hits else None


def _distance_query(positive, negative):
    distance = min(abs(positive.aEnd - negative.aStart),
                   abs(negative.aEnd - positive.aStart))
    return distance


def _distance_target(positive, negative):
    distance = min(abs(positive.tEnd - negative.tEnd),
                   abs(positive.tStart - negative.tStart))
    return distance


def _is_missed(positive, negative):
    distance = _distance_query(positive, negative)
    if ADAPTER_MIN_LENGTH <= distance and distance <= ADAPTER_MAX_LENGTH:
        log.debug("Missed distance for %s = %i",
                  positive.qName, distance)
        return True
    return False


def _is_missing(positive, negative):
    distance = _distance_query(positive, negative)
    if distance < ADAPTER_MIN_LENGTH:
        log.debug("Missing distance for %s = %i" % (positive.qId,
                                                    distance))
        return True
    return False


def _is_overhang(positive, negative):
    distance = _distance_target(positive, negative)
    if OVERHANG_MIN_DISTANCE < distance:
        log.debug("Overhang distance for %s = %i" %
                  (positive.qName, distance))
        return True
    return False


def _is_end_straddler(hits):
    hits.sort(key=len, reverse=True)
    if len(hits) < 2:
        return False
    refLength = hits[0].referenceInfo['Length']
    target_beginning = Range(0, CHIMERA_END_TOLERANCE)
    target_ending = Range(refLength - CHIMERA_END_TOLERANCE,
                          refLength)
    target_range1 = Range(hits[0].tStart, hits[0].tEnd)
    target_range2 = Range(hits[1].tStart, hits[1].tEnd)
    return ((target_range1.intersects(target_beginning) and
             target_range2.intersects(target_ending)) or
            (target_range1.intersects(target_ending) and
             target_range2.intersects(target_beginning)))


def _is_chimeric(hits):
    if len(hits) < 2:
        return False
    hits.sort(key=len, reverse=True)
    query_range1 = Range(hits[0].aStart, hits[0].aEnd)
    query_range2 = Range(hits[1].aStart, hits[1].aEnd)
    if len(query_range1.intersect(query_range2)) > CHIMERA_QUERY_MAX_OVERLAP:
        return False
    target_range1 = Range(hits[0].tStart,
                          hits[0].tEnd + CHIMERA_TARGET_DISTANCE)
    target_range2 = Range(hits[1].tStart,
                          hits[1].tEnd + CHIMERA_TARGET_DISTANCE)
    if target_range1.intersects(target_range2):
        return False
    return True


def _to_annotation(hits, min_intersect_length):
    if _is_end_straddler(hits):
        return None
    if _is_chimeric(hits):
        return "Chimera"
    positive = _get_longest_on_strand(hits, strand="+")
    negative = _get_longest_on_strand(hits, strand="-")
    if not negative or not positive:
        return None
    range1 = Range(positive.tStart, positive.tEnd)
    range2 = Range(negative.tStart, negative.tEnd)
    if len(range1.intersect(range2)) < min_intersect_length:
        return None

    if _is_overhang(positive, negative):
        return "Overhang"
    if _is_missed(positive, negative):
        return "Undetected"
    if _is_missing(positive, negative):
        return "Missing"
    return None


def align(read_fn, ref_fn, tmp_dir, nproc=16):
    """Find hits using a modified pbalign run (use pbalign to wrap blasr)"""
    out_fn = os.path.join(tmp_dir, "blasr_hits.bam")
    # using pbalign gives us DataSet XML IO, as well as blasr
    # dependency resolution
    args = [read_fn, ref_fn, out_fn, '--hitPolicy', 'all',
            '--algorithmOptions',
            '-minMatch 10 -maxLCPLength 12 -maxScore -100 -noSplitSubreads '
            '-minFrac 0.001 -useGuidedAlign -nCandidates 10 '
            '-minSubreadLength 0 -nproc {n}'.format(n=nproc)]
    PBAlignRunner(args).start()
    out_fixed = os.path.join(tmp_dir, "blasr_hits.pbver.bam")
    rectify(out_fn, out_fixed)
    return out_fixed


def find_adapter_issues(read_fn, ref_fn, tmp_dir, outfn=None,
                        min_overlap_length=100, nproc=16,
                        min_intersect_length=50):
    log.info("Starting findAdapterIssues")
    hits_by_query = {}
    issues = []
    with openAlignmentFile(align(read_fn, ref_fn, tmp_dir,
                                 nproc=nproc)) as infile:
        for row in infile:
            length = abs(row.aEnd - row.aStart)
            if length < int(min_overlap_length):
                continue
            hits_by_query.setdefault(row.qName, []).append(row)

        issues.append("id,adapter_annotation\n")
        for query, hits in hits_by_query.items():
            annotation = _to_annotation(hits, min_intersect_length)
            if annotation:
                issues.append("%s,%s\n" % (query, annotation))
            elif annotation is None and ANNOTATE_NON_ISSUES:
                issues.append("%s,%s\n" % (query, "AlignedNoIssue"))

    if outfn is None:
        for issue in issues:
            log.info(issue)
        for issue in issues:
            print issue
    else:
        with open(outfn, 'w') as outf:
            outf.writelines(issues)
    return issues


def arg_runner(args):
    """Callback for pbsystem cmdline functions"""
    log.info(
        "Starting {f} version {v} report generation".format(f=__file__,
                                                            v=__version__))
    tmp_dir = tempfile.mkdtemp(suffix="find_adapter_issues_pbalign")
    clean = True
    if args.tmp_dir:
        tmp_dir = args.tmp_dir
        # This doesn't have to be the case, but it makes life simpler.
        clean = False
    if args.debug:
        clean = False
    find_adapter_issues(args.reads, args.reference, tmp_dir, args.output_file,
                        args.min_overlap_length, args.nproc,
                        args.min_intersect_length)
    if clean:
        shutil.rmtree(tmp_dir)
    return 0


def _get_parser():
    """Handle command line argument parsing"""
    description = __doc__
    parser = get_default_argparser_with_base_opts(
        version=__version__,
        description=description)
    parser.add_argument("reads", type=validate_file,
                        help="input reads as fasta, DataSet XML")
    parser.add_argument("reference", type=validate_file,
                        help="input reference as fasta or ReferenceSet XML")
    parser.add_argument("--output_file", type=str,
                        help="output csv file for found issues (default=log)")
    parser.add_argument("--tmp_dir", type=validate_output_dir,
                        help="output directory for temporary files",
                        default=None)
    parser.add_argument("--nproc", default=8, type=int,
                        help="Number of processes to use")
    parser.add_argument("--debug", default=False, action='store_true',
                        help="Keep temporary files")
    parser.add_argument("--alignmentThreshold", type=float, default=7.0,
                        dest='alignment_threshold',
                        help="Minimum alignment score of alignments to accept")
    parser.add_argument("--minOverlapLength", type=int, default=100,
                        dest='min_overlap_length',
                        help="Minimum length of alignments to accept")
    parser.add_argument("--minIntersectLength", type=int, default=50,
                        dest='min_intersect_length',
                        help="Minimum length of intersection to accept")
    return parser


def main(argv=sys.argv):
    """Main point of entry"""
    return pacbio_args_runner(
        argv=argv[1:],
        parser=_get_parser(),
        args_runner_func=arg_runner,
        alog=log,
        setup_log_func=setup_log)
