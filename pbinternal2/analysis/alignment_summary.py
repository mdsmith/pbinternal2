"""POC Dev report of the AlignmentSets by Condition.

Generates a report of the number of Mapped Subread and mapped Subread bases
by condition
"""
import os
from collections import namedtuple, defaultdict
import logging

from pbcommand.models.report import Report, Table, Column, Attribute
from pbcommand.models import ReseqConditions
from pbcommand.pb_io import load_reseq_conditions_from

from pbcore.io.dataset import AlignmentSet

log = logging.getLogger(__name__)

Stats = namedtuple("Status", "nbases nsubreads")


def get_nmapped_subreads(alignmentset_paths):

    nbases = 0
    nsubreads = 0

    for file_name in alignmentset_paths:
        log.info("Analyzing {f}".format(f=file_name))
        # Remove this after more testing
        try:
            a = AlignmentSet(file_name)
            nbases += a.totalLength
            nsubreads += a.numRecords
        except Exception as e:
            log.error("unable to process {f}".format(f=file_name), exc_info=True)

    return Stats(nbases, nsubreads)


def get_stats(conditions):
    """
    :type conditions: List[ReseqCondition]
    """

    # {cond-id : List[ReseqCondition]
    cond_d = defaultdict(list)

    for c in conditions:
        cond_d[c.cond_id].append(c)

    # {cond-id : Stats}
    results = {}

    for cond_id, reseq_cond_list in cond_d.iteritems():
        alignmentset_paths = [a.alignmentset for a in reseq_cond_list]
        stats = get_nmapped_subreads(alignmentset_paths)
        results[cond_id] = stats

    return results


def _get_or(a_list, default):
    try:
        return max(a_list)
    except ValueError:
        return default


def to_report(reseq_conditions):
    """
    :type reseq_conditions: ReseqConditions
    :rtype: Report
    """
    log.info("Analyzing conditions {a}".format(a=reseq_conditions))

    # {cond-id: Stats}
    results = get_stats(reseq_conditions.conditions)

    # create table
    cond_ids = []
    nbases = []
    nsubreads = []

    for c_id, stat in results.iteritems():
        cond_ids.append(c_id)
        nbases.append(stat.nbases)
        nsubreads.append(stat.nsubreads)

    c1 = Column("cond_id", "Condition Id", cond_ids)
    c2 = Column("nbases", "Number of Mapped Bases", nbases)
    c3 = Column("nsubreads", "Number of Mapped Subreads", nsubreads)
    columns = [c1, c2, c3]

    table = Table("cond_alignmentset_summary",
                  title="Condition AlignmentSet Summary",
                  columns=columns)

    nsubreads = [x.nsubreads for x in results.values()]
    nbases = [x.nbases for x in results.values()]

    # this should probably be None
    max_nsubreads = _get_or(nsubreads, 0)
    max_nbases = _get_or(nbases, 0)

    attr_nrecords = Attribute("max_nbases", max_nbases, name="Max number of Bases (by Condition)")
    attr_nconditions = Attribute("max_nsubreads", max_nsubreads, name="Max number of Subreads (by Condition)")

    attributes = [attr_nrecords, attr_nconditions]

    r = Report("cond_aln_summary_rpt",
               title="Condition AlignmentSet Summary",
               attributes=attributes,
               tables=(table, ))

    log.info(r)
    for attr in r.attributes:
        log.info(attr)

    for t in r.tables:
        log.info("\n{t}".format(t=t))
    return r


def run_cond_alignmentsets_to_report(condition_json, output_report):

    c = load_reseq_conditions_from(condition_json)
    r = to_report(c)
    r.write_json(output_report)

    return 0
