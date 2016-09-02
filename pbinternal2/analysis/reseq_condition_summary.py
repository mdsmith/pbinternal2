"""Hello World-ish util to create a report from a Condition JSON

The Report Table will show a summary of the conditions
"""
import os
import logging

from pbcommand.models.report import Report, Table, Column, Attribute
from pbcommand.models import ReseqConditions
from pbcommand.pb_io import load_reseq_conditions_from

log = logging.getLogger(__name__)


def to_report(reseq_conditions):
    """
    :type reseq_conditions: ReseqConditions
    :rtype: Report
    """

    cond_ids = []
    aln_paths = []
    subread_paths = []
    path_exists = []

    for c in reseq_conditions.conditions:
        cond_ids.append(c.cond_id)
        subread_paths.append(c.subreadset)
        aln_paths.append(c.alignmentset)
        path_exists.append(os.path.exists(c.alignmentset))

    c1 = Column("cond_id", "Condition Id", cond_ids)
    c2 = Column("file_exists", "File Exists", path_exists)
    c3 = Column("sub_path", "SubreadSet Path", subread_paths)
    c4 = Column("aln_path", "AlignmentSet Path", aln_paths)
    columns = [c1, c2, c3, c4]

    table = Table("cond_summary",
                  title="Condition Summary",
                  columns=columns)

    cond_ids = set(cond_ids)
    attr_nrecords = Attribute("nrecords", len(reseq_conditions.conditions), name="Number of Records")
    attr_nconditions = Attribute("nconditions", len(cond_ids), name="Number of Conditions")

    attributes = [attr_nrecords, attr_nconditions]

    r = Report("cond_summary_rpt",
               title="Condition Summary",
               attributes=attributes,
               tables=(table, ))

    log.info(r)
    for attr in r.attributes:
        log.info(attr)

    for t in r.tables:
        log.info("\n{t}".format(t=t))
    return r


def run_cond_to_report(condition_json, output_report):

    c = load_reseq_conditions_from(condition_json)
    r = to_report(c)
    r.write_json(output_report)

    return 0
