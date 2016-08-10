#!/usr/bin/env python

import sys
import logging
from pbcommand.models import FileTypes, OutputFileType
from pbcommand.cli import registry_builder, registry_runner
from pbinternal2.report.eol_qc_stats import eol_qc_stats
from pbinternal2 import TOOL_NAMESPACE

__version__ = "0.1.0"
__author__ = "Martin Smith"

log = logging.getLogger(__name__)

registry = registry_builder(TOOL_NAMESPACE, "python -m pbinternal2.tasks.eol_qc")


def _eol_qc_outputs():
    csv_file_type_id = FileTypes.CSV.file_type_id
    f1 = OutputFileType(csv_file_type_id, "csv_0", "Per Zmw Stats ", "Per Zmw Statistics", "file_per_zmw")
    f2 = OutputFileType(csv_file_type_id, "csv_1", "Per Movie Stats", "Per Movie Statistics", "file_per_movie")
    return f1, f2


@registry("eol_qc", __version__,
          (FileTypes.DS_SUBREADS, FileTypes.DS_ALIGN),
          _eol_qc_outputs(),
          nproc=1, is_distributed=True, options=dict(nreads=32768))
def run_rtc(rtc):
    """
    Run an EOL-QC analysis on an subreadset and alignmentset.
    """
    num_holes_id = '{}.task_options.nreads'.format(TOOL_NAMESPACE)
    num_holes = rtc.task.options[num_holes_id]

    return eol_qc_stats(rtc.task.input_files[0], rtc.task.input_files[1],
                        rtc.task.output_files[0], rtc.task.output_files[1],
                        num_holes)


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
