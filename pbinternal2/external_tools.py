#!/usr/bin/env python

"""An interface for external tools and silent dependencies"""

import logging
from pbcore.util.Process import backticks
import os
import sys

log = logging.getLogger(__name__)


def gepard(max_width, word_length, ref_fn, fasta_fn, png_fn):
    """Use gepard to generate the dotplot for dot.py, or makeDotReport"""

    matrix = ""

    # Build
    if os.getenv("ASSEMBLY_HOME"):
        matrix = os.getenv("ASSEMBLY_HOME") + '/etc/gepard-1.21/edna.mat'
        cmd = ("gepardcmd.sh "
               "-maxwidth %i -word %i -matrix "
               "%s -seq1 %s -seq2 %s -outfile %s" % (
                   max_width, int(word_length),
                   matrix, ref_fn, fasta_fn,
                   png_fn))
    # Local testing (GEPARD must be defined)
    elif os.getenv("GEPARD"):
        matrix = os.getenv("GEPARD") + '/matrices/edna.mat'
        cmd = ("java -cp %s/lib/gepard.jar "
               "org.krumsiek.gepard.client.cmdline.CommandLine "
               "-maxwidth %i -word %i -matrix "
               "%s -seq1 %s -seq2 %s -outfile %s" % (
                   os.getenv("GEPARD"), max_width, int(word_length),
                   matrix, ref_fn, fasta_fn,
                   png_fn))
    # Jenkins testing using Chef build
    else:
        matrix = '/usr/local/bin/gepard-1.30/matrices/edna.mat'
        cmd = ("java -cp /usr/local/bin/gepard-1.30/lib/gepard.jar "
               "org.krumsiek.gepard.client.cmdline.CommandLine "
               "-maxwidth %i -word %i -matrix "
               "%s -seq1 %s -seq2 %s -outfile %s" % (
                   max_width, int(word_length),
                   matrix, ref_fn, fasta_fn,
                   png_fn))
    log.error("Running %s" % cmd)
    backticks(cmd)


def bamPbi(in_bam, referenceFasta):

    inPath = False
    test = backticks("makePbi")
    if test[0] or test[1] != 127:
        inPath = True

    # Build
    if os.getenv("SMRT_ROOT"):
        makePbi_exec = os.path.join(os.getenv("SMRT_ROOT"), "bin", "makePbi")
        cmd = ("%s %s --referenceFasta %s" % (makePbi_exec, in_bam,
                                              referenceFasta))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Local testing (blasr must be in path)
    if inPath:
        makePbi_exec = 'makePbi'
        cmd = ("%s %s --referenceFasta %s" % (makePbi_exec, in_bam,
                                              referenceFasta))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Jenkins testing using Chef build
    makePbi_exec = '/usr/local/bin/makePbi'
    cmd = ("%s %s --referenceFasta %s" % (makePbi_exec, in_bam,
                                          referenceFasta))
    log.debug("Running %s" % cmd)
    results = backticks(cmd)
    if not results[1]:
        return results[0]
    else:
        log.debug("%s Returned non-zero" % cmd)
    sys.exit("makePbi not found in PATH, SMRT_ROOT/bin/ or /usr/local/bin/, "
             "exiting")


def bamIndex(in_bam):

    inPath = False
    test = backticks("samtools")
    if test[0] or test[1] == 1:
        inPath = True

    # Build
    if os.getenv("SMRT_ROOT"):
        samtools_exec = os.path.join(os.getenv("SMRT_ROOT"), "bin", "samtools")
        cmd = ("%s index %s" % (samtools_exec, in_bam))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Local testing (blasr must be in path)
    if inPath:
        samtools_exec = 'samtools'
        cmd = ("%s index %s" % (samtools_exec, in_bam))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Jenkins testing using Chef build
    samtools_exec = '/usr/local/bin/samtools'
    cmd = ("%s index %s" % (samtools_exec, in_bam))
    log.debug("Running %s" % cmd)
    results = backticks(cmd)
    if not results[1]:
        return results[0]
    else:
        log.debug("%s Returned non-zero" % cmd)
    sys.exit("samtools not found/successful in PATH, SMRT_ROOT/bin/ or "
             "/usr/local/bin/, " "exiting")


def fasIndex(in_fas):

    inPath = False
    test = backticks("samtools")
    if test[0] or test[1] == 1:
        inPath = True

    # Build
    if os.getenv("SMRT_ROOT"):
        samtools_exec = os.path.join(os.getenv("SMRT_ROOT"), "bin", "samtools")
        cmd = ("%s faidx %s" % (samtools_exec, in_fas))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Local testing (blasr must be in path)
    if inPath:
        samtools_exec = 'samtools'
        cmd = ("%s faidx %s" % (samtools_exec, in_fas))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Jenkins testing using Chef build
    samtools_exec = '/usr/local/bin/samtools'
    cmd = ("%s faidx %s" % (samtools_exec, in_fas))
    log.debug("Running %s" % cmd)
    results = backticks(cmd)
    if not results[1]:
        return results[0]
    else:
        log.debug("%s Returned non-zero" % cmd)
    sys.exit("samtools not found/successful in PATH, SMRT_ROOT/bin/ or "
             "/usr/local/bin/, " "exiting")


def blasr(read_fn, reference_fn, out_fn, maxMatch=12, minMatch=10,
          minPctIdentity=0.001, nproc=14):

    inPath = False
    test = backticks("blasr")
    if test[0] or test[1] != 127:
        inPath = True

    # Build
    if os.getenv("SMRT_ROOT"):
        blasr_exec = os.path.join(os.getenv("SMRT_ROOT"), "bin", "blasr")
        cmd = ("%s %s %s -maxLCPLength %i -maxScore -100 -minFrac 0.001 "
               "-noSplitSubreads -useGuidedAlign -nCandidates 10 -minMatch %i "
               "-minPctIdentity %i -nproc %i -bam -out %s" % (
                   blasr_exec, read_fn, reference_fn, maxMatch, minMatch,
                   minPctIdentity, nproc, out_fn))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Local testing (blasr must be in path)
    if inPath:
        blasr_exec = 'blasr'
        cmd = ("%s %s %s -maxLCPLength %i -maxScore -100 -minFrac 0.001 "
               "-noSplitSubreads -useGuidedAlign -nCandidates 10 -minMatch %i "
               "-minPctIdentity %i -nproc %i -bam -out %s" % (
                   blasr_exec, read_fn, reference_fn, maxMatch, minMatch,
                   minPctIdentity, nproc, out_fn))
        log.debug("Running %s" % cmd)
        results = backticks(cmd)
        if not results[1]:
            return results[0]
        else:
            log.debug("%s Returned non-zero" % cmd)

    # Jenkins testing using Chef build
    blasr_exec = '/usr/local/bin/blasr'
    cmd = ("%s %s %s -maxLCPLength %i -maxScore -100 -minFrac 0.001 "
           "-noSplitSubreads -useGuidedAlign -nCandidates 10 -minMatch %i "
           "-minPctIdentity %i -nproc %i -bam -out %s" % (
               blasr_exec, read_fn, reference_fn, maxMatch, minMatch,
               minPctIdentity, nproc, out_fn))
    log.debug("Running %s" % cmd)
    results = backticks(cmd)
    if not results[1]:
        return results[0]
    else:
        log.debug("%s Returned non-zero" % cmd)
    sys.exit("blasr not found in PATH, SMRT_ROOT/bin/ or /usr/local/bin/, "
             "exiting")
