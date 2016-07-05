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

__author__ = "Nat Echols"


class Constants(object):

    DRIVER_BASE = "python -m pbinternal2.pa_tasks "

    BASECALLER_OPTIONS = "--internal --method=TA_DME_FFHmm_P2B"
    BASECALLER_OPTIONS_ID = "basecaller_options"
    BASECALLER_EXE = "basecaller-console-app"
    BASECALLER_EXE_ID = "basecaller_exe"
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


def run_baz2bam(baz_file, adapter_fa, metadata_xml, output_file,
                nproc=1,
                min_subread_length=Constants.MIN_SUBREAD_LENGTH,
                baz2bam_exe=Constants.BAZ2BAM_EXE,
                stdout=sys.stdout, stderr=sys.stderr):
    """
    Convert the .baz file from the basecaller to a SubreadSet.
    """
    assert output_file.endswith(".subreadset.xml")
    output_base = re.sub(".subreadset.xml", "", output_file)
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
    with SubreadSet(subreadset_file) as ds:
        ds.makePathsAbsolute()
        ds.write(tmp_ds)
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
          options={
              Constants.BASECALLER_EXE_ID: Constants.BASECALLER_EXE,
              Constants.BASECALLER_OPTIONS_ID: Constants.BASECALLER_OPTIONS,
          })
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


if __name__ == '__main__':
    # FIXME. This logging should be able to be removed
    logging.basicConfig(level=logging.INFO)
    sys.exit(registry_runner(registry, sys.argv[1:], default_log_level=logging.INFO))