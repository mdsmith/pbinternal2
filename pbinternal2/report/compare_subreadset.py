"""Primary Analysis Util for Comparing SubreadSets

Specific Usecase is for comparing refarmed data.

"""
import os
import logging

from pbcore.io.dataset import SubreadSet
from pbcommand.models.report import Report, Attribute

log = logging.getLogger(__name__)

__version__ = "0.1.0"


def compare_subreadset_to_report(subreadset_a_path, subreadset_b_path):
    """Subreadset A is assumed to be the baseline value

    :type subreadset_a_path: basestring
    :type subreadset_b_path: basestring

    :rtype Report
    """
    # MK. Adding the scaffolding to get the pipeline implemented

    report_id = "subreadset_compare"
    report_title = "SubreadSet Comparison"

    a = SubreadSet(subreadset_a_path)
    b = SubreadSet(subreadset_b_path)

    delta_nrecords = b.numRecords - a.numRecords
    delta_nlength = b.totalLength - a.totalLength
    dataset_uuids = [a.uuid, b.uuid]

    log.info("SubreadSet A {}".format(a))
    log.info("SubreadSet B {}".format(b))

    attrs = [
        ("compare_version", __version__, "SubreadSet Compare Version"),
        ("delta_nrecords", delta_nrecords, "Delta NumRecords (B - A)"),
        ("delta_ntotal_length", delta_nlength, "Delta TotalLength (B -A)")
    ]

    attributes = [Attribute(i, v, name=n) for i, v, n in attrs]

    r = Report(report_id, report_title,
               attributes=attributes,
               dataset_uuids=dataset_uuids)

    log.info("Report {}".format(r))
    return r


def run_compare_subreadset_compare(path_a, path_b, report_output):
    r = compare_subreadset_to_report(path_a, path_b)
    r.write_json(report_output)
    return 0