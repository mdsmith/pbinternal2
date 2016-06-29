"""Central Location for TC/RTC for pbinternal2.analysis.* """
import sys
import logging

from pbcommand.cli import registry_runner, registry_builder
from pbcommand.models import FileTypes

from pbinternal2.analysis.alignment_summary import run_cond_alignmentsets_to_report
from pbinternal2.analysis.reseq_condition_summary import run_cond_to_report
from pbinternal2.analysis.hello_world import run_hello_world
from pbinternal2.analysis.yield_plots import run_yield_plots

TOOL_NAMESPACE = 'pbinternal2'
DRIVER_BASE = 'python -m pbinternal2.analysis_tools '

registry = registry_builder(TOOL_NAMESPACE, DRIVER_BASE)


@registry("hello_world_analysis", "0.1.0", FileTypes.FASTA, FileTypes.REPORT, is_distributed=True, nproc=1)
def run_rtc(rtc):
    return run_hello_world(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("yield_plot_analysis", "0.1.1", FileTypes.COND_RESEQ, FileTypes.REPORT, is_distributed=True, nproc=1)
def run_rtc(rtc):
    return run_yield_plots(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("cond_to_report", "0.2.0", FileTypes.COND_RESEQ, FileTypes.REPORT, nproc=1, is_distributed=False)
def run_rtc(rtc):
    return run_cond_to_report(rtc.task.input_files[0], rtc.task.output_files[0])


@registry("cond_to_alignmentsets_report", "0.2.0", FileTypes.COND_RESEQ, FileTypes.REPORT, nproc=1, is_distributed=True)
def run_rtc(rtc):
    return run_cond_alignmentsets_to_report(rtc.task.input_files[0], rtc.task.output_files[0])


if __name__ == '__main__':
    default_log_level = logging.DEBUG
    sys.exit(registry_runner(registry,
                             sys.argv[1:],
                             default_log_level=default_log_level))


if __name__ == '__main__':
    sys.exit(registry_runner(registry, sys.argv[1:]))
