#!/usr/bin/env python

import sys
import logging
from pbcommand.models import FileTypes
from pbcommand.cli import registry_builder, registry_runner
from pbinternal2.report.eol_qc_stats import eol_qc_stats
from pbinternal2 import TOOL_NAMESPACE

__version__ = "0.1.0"
__author__ = "Martin Smith"

log = logging.getLogger(__name__)

registry = registry_builder(TOOL_NAMESPACE, "python -m pbinternal2.tasks.eol_qc")

@registry("eol_qc", __version__, (FileTypes.DS_SUBREADS, FileTypes.DS_ALIGN), (FileTypes.CSV, FileTypes.CSV),
          nproc=1, is_distributed=False, options=dict(nreads=32768))
def run_rtc(rtc):
    """
    Run an EOL-QC analysis on an subreadset and alignmentset.
    """

    return eol_qc_stats(rtc.task.input_files[0], rtc.task.input_files[1],
                        rtc.task.output_files[0], rtc.task.output_files[1],
                        rtc.task.options['{:}.task_options.nreads'.format(TOOL_NAMESPACE)])


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
