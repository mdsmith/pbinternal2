"""Task to Find the Trace file from the SubreadSet"""

import os
import logging

from pbcore.io.dataset import SubreadSet

log = logging.getLogger(__name__)


def subreadset_to_trace_file(sset_or_path):
    """

    Attempts to resolve the trace file from the SubreadSet.

    Note, this only works with non-merged SubreadSets.

    :param sset_or_path: Path to SubreadSet or
    :type sset_or_path: basestring | SubreadSet
    :return: Path to trace file

    :raises Assertion, IOError and any DataSet IO related errors
    """
    log.info("Attempting to resolve Trace file from {}".format(sset_or_path))

    subread_bam_ext = ".subreads.bam"
    trace_ext = ".trc.h5"

    if isinstance(sset_or_path, basestring):
        sset = SubreadSet(sset_or_path)
    elif isinstance(sset_or_path, SubreadSet):
        sset = sset_or_path
    else:
        raise TypeError("Expected path or SubreadSet")

    fs = sset.toExternalFiles()
    # not a merged dataset
    assert fs[0].endswith(subread_bam_ext)
    assert len(fs) == 1
    # MK. In future version of Primary Analysis, the SubreadSet will be written
    # with an explicit path to the Trace file as an external Resource
    # within the SubreadSet.
    trace_file = fs[0].replace(subread_bam_ext, trace_ext)
    assert os.path.isfile(trace_file)

    log.info("Successfully resolved trace file {}".format(trace_file))
    return os.path.abspath(trace_file)