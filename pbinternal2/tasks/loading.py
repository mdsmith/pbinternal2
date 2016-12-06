"""
Primary/Post-Primary analysis tasks for internal testing with pbsmrtpipe.
"""

import logging
import sys

from pbcommand.cli import registry_builder, registry_runner
from pbcommand.models import FileTypes

from pbinternal2 import TOOL_NAMESPACE
from pbinternal2.report.LoadingEvaluator import loading_vs_poisson

__author__ = "Martin Smith"

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler()) # To avoid warning messages


class Constants(object):

    DRIVER_BASE = "python -m pbinternal2.tasks.loading "

    LOADING_DIST_ID = "loading_dist"
    LOADING_DIST = 20


def _get_id(base_opt):
    return "{n}.task_options.{o}".format(n=TOOL_NAMESPACE,
                                         o=base_opt)

registry = registry_builder(TOOL_NAMESPACE, Constants.DRIVER_BASE)


@registry("loading_vs_poisson_report", "0.1.0",
          FileTypes.STS_H5,
          FileTypes.REPORT,
          is_distributed=True,
          nproc=1,
          options={Constants.LOADING_DIST_ID: Constants.LOADING_DIST})
def task_loading_vs_poisson_report(rtc):
    return loading_vs_poisson(
        rtc.task.input_files[0],
        rtc.task.output_files[0],
        rtc.task.nproc,
        rtc.task.options[_get_id(Constants.LOADING_DIST_ID)])

if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:],
                             default_log_level=logging.INFO))

