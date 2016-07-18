"""Primary Analysis Util for Comparing SubreadSets

Specific Usecase is for comparing refarmed data.

"""
import os
import logging

from pbcore.io.dataset import SubreadSet
from pbcommand.models.report import Report, Attribute

log = logging.getLogger(__name__)


def compare_subreadset_to_report(subreadset_a_path, subreadset_b_path):
    """Subreadset A is assumed to be the baseline value

    :type subreadset_a_path: basestring
    :type subreadset_b_path: basestring

    :rtype Report
    """
    # MK. Adding the scaffolding to get the pipeline implemented

    a = SubreadSet(subreadset_a_path)
    b = SubreadSet(subreadset_b_path)

    log.info("SubreadSet A {}".format(a))
    log.info("SubreadSet B {}".format(b))

    mock_attribute = Attribute("test", "Not Implemented Yet")

    r = Report("subreadset_compare", "SubreadSet Comparison",
               attributes=[mock_attribute])

    log.info("Report {}".format(r))
    return r


def run_compare_subreadset_compare(path_a, path_b, report_output):
    r = compare_subreadset_to_report(path_a, path_b)
    r.write_json(report_output)
    return 0