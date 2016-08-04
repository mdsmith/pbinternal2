"""
Primary/Post-Primary analysis tasks for internal testing with pbsmrtpipe.
"""

import subprocess
import tempfile
import logging
import shutil
import re
import os.path as op
import sys

from pbcommand.engine import run_cmd
from pbcommand.cli import registry_builder, registry_runner
from pbcommand.models import FileTypes, SymbolTypes
from pbcore.io import SubreadSet

from pbinternal2 import TOOL_NAMESPACE
from pbinternal2.report.compare_subreadset import run_compare_subreadset_compare
from pbinternal2.util.subreadset_to_trace import subreadset_to_trace_file

__author__ = "Nat Echols"

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler()) # To avoid warning messages


class Constants(object):

    DRIVER_BASE = "python -m pbinternal2.pa_tasks "

    BASECALLER_OPTIONS = "--internal --method=TA_DME_FFHmm_P2B"
    BASECALLER_OPTIONS_ID = "basecaller_options"
    BASECALLER_EXE = "basecaller-console-app"
    BASECALLER_EXE_ID = "basecaller_exe"

    BASE_CALLER_OPTIONS = {BASECALLER_EXE_ID: BASECALLER_EXE,
                           BASECALLER_OPTIONS_ID: BASECALLER_OPTIONS}

    BAZ2BAM_EXE = "baz2bam"
    BAZ2BAM_EXE_ID = "baz2bam_exe"
    MIN_SUBREAD_LENGTH_ID = "min_subread_length"
    MIN_SUBREAD_LENGTH = 50
    BAM2BAM_EXE = "bam2bam"
    BAM2BAM_EXE_ID = "bam2bam_exe"
    BLASR_OPTIONS_ID = "blasr_options"
    BLASR_OPTIONS = "--bestn 1 --forwardOnly --fastMaxInterval --maxAnchorsPerPosition 30000 --ignoreHQRegions --minPctIdentity 60"


def _get_id(base_opt):
    return "{n}.task_options.{o}".format(n=TOOL_NAMESPACE,
                                         o=base_opt)

registry = registry_builder(TOOL_NAMESPACE, Constants.DRIVER_BASE)


def run_basecaller(trc_file, baz_file,
                   nproc=1,
                   stdout=sys.stdout,
                   stderr=sys.stderr,
                   basecaller_exe=Constants.BASECALLER_EXE,
                   basecaller_options=Constants.BASECALLER_OPTIONS):
    """
    Run the offline basecaller on a trace file.
    """
    args = [
        basecaller_exe,
        "--inputfile={i}".format(i=trc_file),
        "--outputbazfile={o}".format(o=baz_file),
        "--numthreads={n}".format(n=nproc),
    ]
    if basecaller_options != "":
        args.extend(basecaller_options.split(' '))
    logging.info(' '.join(args))
    result = run_cmd(' '.join(args), stdout_fh=stdout, stderr_fh=stderr)
    assert op.isfile(baz_file)
    return result.exit_code


def _adapters_to_str(left, right):
    # This is similar to how Paws works. Paws will call toUppercase for
    # reasons that are unclear to me. Not doing that here.
    if left == right:
        return ">LeftAdapter\n{}".format(left)
    else:
        return ">LeftAdapter\n{l}\n>RightAdapter\n{r}".format(l=left, r=right)


def write_adapters(left, right, output_file):
    with open(output_file, 'w') as out:
        out.write(_adapters_to_str(left, right))
    return output_file


def write_adapters_from_subreadset(subreadset_path, output_file):
    with SubreadSet(subreadset_path) as ds:
        right = ds.metadata.collections[0].templatePrepKit.rightAdaptorSequence
        left = ds.metadata.collections[0].templatePrepKit.leftAdaptorSequence
        write_adapters(left, right, output_file)
    log.info("wrote adapters to {}".format(output_file))



def run_baz2bam(baz_file, adapter_fa, metadata_xml, output_file,
                nproc=1,
                min_subread_length=Constants.MIN_SUBREAD_LENGTH,
                baz2bam_exe=Constants.BAZ2BAM_EXE,
                stdout=sys.stdout, stderr=sys.stderr,
                dataset_name_suffix=None):
    """
    Convert the .baz file from the basecaller to a SubreadSet.

    Note, the emitted SubreadSet will have a new UUID

    :param output_file: Base prefix for output files

    :param dataset_name_suffix: Will update the dataset name with the supplied suffix
    :type dataset_name_suffix: str | None
    """

    assert output_file.endswith(".subreadset.xml")
    output_base = re.sub(".subreadset.xml", "", output_file)
    output_dir = op.dirname(output_file)

    args = [
        baz2bam_exe,
        baz_file,
        "--silent",
        "--minSubLength", str(min_subread_length),
        "--metadata={x}".format(x=metadata_xml),
        "--adapter={f}".format(f=adapter_fa),
        "-o", output_base,
        "-j", str(nproc)
    ]
    logging.info(" ".join(args))
    result = run_cmd(' '.join(args), stdout, stderr)
    assert result.exit_code == 0, \
        "Failed with exit code {c}".format(c=result.exit_code)
    subreads_file = output_base + ".subreads.bam"
    scraps_file = output_base + ".scraps.bam"
    assert op.isfile(subreads_file), subreads_file
    assert op.isfile(scraps_file), scraps_file
    subreadset_file = output_base + ".subreadset.xml"
    assert op.isfile(subreadset_file)
    tmp_ds = tempfile.NamedTemporaryFile(suffix=".subreadset.xml").name

    # Must copy the adapters file to new SubreadSet output dir
    # otherwise, the file will be invalid
    new_adapters = op.join(output_dir, op.basename(adapter_fa))
    shutil.copy(adapter_fa, new_adapters)

    # FIXME, This should really update the PA version (SigProcVer) or at a minimum,
    # augment the version

    with SubreadSet(subreadset_file) as ds:
        ds.makePathsAbsolute()
        if dataset_name_suffix is not None:
            name = ds.name
            new_ds_name = "_".join([name, dataset_name_suffix])
            ds.name = new_ds_name
        ds.newUuid(setter=True)
        ds.write(tmp_ds)
        log.info("Wrote new SubreadSet {u} to {p}".format(u=ds.uuid, p=subreads_file))

    shutil.move(tmp_ds, subreadset_file)
    return 0


def run_bam2bam(subreads_ds, bam_file,
                bam2bam_exe=Constants.BAM2BAM_EXE):
    """
    """
    prefix = re.sub(".bam$", "", bam_file)
    ds = SubreadSet(subreads_ds)
    if len(ds.externalResources) != 1:
        raise ValueError("This operation is only supported for single-BAM " +
                         "SubreadSets.")
    subreads_bam = ds.externalResources[0].bam
    scraps_bam = ds.externalResources[0].scraps
    args = [
        bam2bam_exe,
        "--polymerase",
        "-o", prefix,
        subreads_bam,
        scraps_bam,
    ]
    assert subprocess.call(args) == 0
    bam_out = prefix + ".polymerase.bam"
    shutil.move(bam_out, bam_file)
    return 0


def run_bax2bam_polymerase(hdfsubreads_ds, bam_file):
    prefix = re.sub(".bam$", "", bam_file)
    args = [
        "bax2bam",
        "--polymeraseread",
        "-o", prefix,
        "--xml",
        hdfsubreads_ds,
    ]
    assert subprocess.call(args) == 0
    bam_out = prefix + ".polymerase.bam"
    shutil.move(bam_out, bam_file)
    return 0


def run_pbalign_unrolled(bam_file, reference, alignment_set,
                         nproc=1,
                         blasr_options=Constants.BLASR_OPTIONS):
    args = [
        "pbalign",
        "--noSplitSubreads",
        "--hitPolicy=leftmost",
        "--nproc", str(nproc),
        "--algorithmOptions", blasr_options,
        bam_file,
        reference,
        alignment_set,
    ]
    assert subprocess.call(args) == 0
    assert op.isfile(alignment_set)
    return 0


@registry("basecaller", "0.1.0",
          FileTypes.TRC,
          FileTypes.BAZ,
          is_distributed=True,
          nproc=SymbolTypes.MAX_NPROC,
          options=Constants.BASE_CALLER_OPTIONS)
def task_basecaller(rtc):
    return run_basecaller(
        trc_file=rtc.task.input_files[0],
        baz_file=rtc.task.output_files[0],
        nproc=rtc.task.nproc,
        basecaller_options=rtc.task.options[
            _get_id(Constants.BASECALLER_OPTIONS_ID)],
        basecaller_exe=rtc.task.options[_get_id(Constants.BASECALLER_EXE_ID)])


@registry("baz2bam", "0.1.0",
          (FileTypes.BAZ, FileTypes.FASTA, FileTypes.XML),
          FileTypes.DS_SUBREADS,
          is_distributed=True,
          nproc=SymbolTypes.MAX_NPROC,
          options={
              Constants.BAZ2BAM_EXE_ID: Constants.BAZ2BAM_EXE,
              Constants.MIN_SUBREAD_LENGTH_ID: Constants.MIN_SUBREAD_LENGTH,
          })
def task_baz2bam(rtc):
    return run_baz2bam(
        baz_file=rtc.task.input_files[0],
        adapter_fa=rtc.task.input_files[1],
        metadata_xml=rtc.task.input_files[2],
        output_file=rtc.task.output_files[0],
        min_subread_length=rtc.task.options[
            _get_id(Constants.MIN_SUBREAD_LENGTH_ID)],
        nproc=rtc.task.nproc,
        baz2bam_exe=rtc.task.options[_get_id(Constants.BAZ2BAM_EXE_ID)])


@registry("bam2bam", "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.BAM,
          is_distributed=True,
          nproc=1,
          options={Constants.BAM2BAM_EXE_ID: Constants.BAM2BAM_EXE})
def task_bam2bam(rtc):
    return run_bam2bam(
        subreads_ds=rtc.task.input_files[0],
        bam_file=rtc.task.output_files[0],
        bam2bam_exe=rtc.task.options[_get_id(Constants.BAM2BAM_EXE_ID)])


@registry("bax2bam_polymerase", "0.1.0",
          FileTypes.DS_SUBREADS_H5,
          FileTypes.BAM,
          is_distributed=True,
          nproc=1,
          options={})
def task_bax2bam_polymerase(rtc):
    return run_bax2bam_polymerase(
        hdfsubreads_ds=rtc.task.input_files[0],
        bam_file=rtc.task.output_files[0])


@registry("pbalign_unrolled", "0.1.0",
          (FileTypes.BAM, FileTypes.DS_REF),
          FileTypes.DS_ALIGN,
          is_distributed=True,
          nproc=SymbolTypes.MAX_NPROC,
          options={Constants.BLASR_OPTIONS_ID: Constants.BLASR_OPTIONS})
def task_pbalign_unrolled(rtc):
    return run_pbalign_unrolled(
        bam_file=rtc.task.input_files[0],
        reference=rtc.task.input_files[1],
        alignment_set=rtc.task.output_files[0],
        nproc=rtc.task.nproc,
        blasr_options=rtc.task.options[_get_id(Constants.BLASR_OPTIONS_ID)])


@registry("basecaller_from_subreadset",
          "0.1.0",
          FileTypes.DS_SUBREADS,
          FileTypes.BAZ,
          is_distributed=True,
          nproc=1,
          options=Constants.BASE_CALLER_OPTIONS,
          name="SubreadSet/Trace Refarm")
def run_base_caller_from_subreadset(rtc):
    """Attempt to Resolve the Trace file from the SubreadSet, then call basecaller"""

    subreadset_path = rtc.task.input_files[0]
    trace_path = subreadset_to_trace_file(subreadset_path)
    baz_output_file = rtc.task.output_files[0]

    base_caller_options = rtc.task.options[_get_id(Constants.BASECALLER_OPTIONS_ID)]
    basecaller_exe = rtc.task.options[_get_id(Constants.BASECALLER_EXE_ID)]

    return run_basecaller(trace_path, baz_output_file,
                          nproc=rtc.task.nproc,
                          basecaller_options=base_caller_options,
                          basecaller_exe=basecaller_exe)


@registry("baz2subreadset", "0.2.1",
          (FileTypes.DS_SUBREADS, FileTypes.BAZ),
          (FileTypes.DS_SUBREADS, ),
          is_distributed=True, name="Baz To SubreadSet")
def run_baz2subreadset(rtc):
    """Convert the Baz File to a SubreadSet.
    The Output is the Refarmed SubreadSet with a new UUID and the name
    will have a suffix of "refarmed"

    inputs: (Original SubreadSet, Input Baz)
    outputs: (Output SubreadSet, )
    """

    original_subreadset = rtc.task.input_files[0]
    baz = rtc.task.input_files[1]
    output_subreadset = rtc.task.output_files[0]

    output_dir = op.dirname(output_subreadset)

    nproc = rtc.task.nproc

    # Crude resolving model
    def fx(ext):
        return original_subreadset.replace(".subreadset.xml", ext)

    adapters_fasta = fx(".adapters.fasta")
    run_metadata_xml = fx(".run.metadata.xml")

    # If the adapters file isn't present (i.e., it wasn't transferred off
    # the Inst), then write the adapter file to the new output dir
    if not op.exists(adapters_fasta):
        adapters_fasta = op.join(output_dir, op.basename(adapters_fasta))
        write_adapters_from_subreadset(original_subreadset, adapters_fasta)

    return run_baz2bam(baz, adapters_fasta, run_metadata_xml,
                       output_subreadset,
                       nproc=nproc,
                       dataset_name_suffix="refarmed")


@registry("compare_subreadsets_report", "0.2.0",
          (FileTypes.DS_SUBREADS, FileTypes.DS_SUBREADS), FileTypes.REPORT,
          is_distributed=True, nproc=1, name="SubreadSet Compare")
def run_compare_subreadsets(rtc):
    """Compare Two SubreadSets. The first value is assumed to be the baseline"""
    return run_compare_subreadset_compare(rtc.task.input_files[0],
                                          rtc.task.input_files[1],
                                          rtc.task.output_files[0])


if __name__ == '__main__':
    # FIXME. This logging should be able to be removed
    logging.basicConfig(level=logging.INFO)
    sys.exit(registry_runner(registry, sys.argv[1:], default_log_level=logging.INFO))
