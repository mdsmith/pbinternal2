#!/usr/bin/env python
"""Automation util to submit Analysis Jobs to SMRT Link from a SubreadSet Path, Rset UUID and pipeline id.

This util will Run (and block) an import-dataset job to import the SubreadSet (if it's not already imported).
Then it will check and see if an analysis job has already been run with the SubreadSet (to avoid submitting a duplicate job)
and it will create a new Analysis Job (if one hasn't already been run).

Only Resequencing based pipelines (i.e., pipelines that start from only a subreadset and reference set) are supported.

"""
# Template pulled from https://github.com/PacificBiosciences/pbcommand/blob/master/pbcommand/cli/examples/template_simple.py

import sys
import logging
from pbcommand.models import FileTypes

from pbcommand.validators import validate_file
from pbcommand.utils import setup_log
from pbcommand.cli import get_default_argparser_with_base_opts, pacbio_args_runner
from pbcommand.services import ServiceAccessLayer, ServiceEntryPoint, JobStates

from pbcore.io.dataset import SubreadSet

log = logging.getLogger(__name__)

__author__ = "M. Kocher"
__version__ = "0.2.1"


class Constants(object):
    DEFAULT_PIPELINE = "pbsmrtpipe.pipelines.sa3_ds_resequencing"
    DEFAULT_RSET_UUID = "c2f73d21-b3bf-4b69-a3a3-6beeb7b22e26"


class ImportJobNotFoundException(BaseException):
    pass


class AnalysisJobNotFoundException(BaseException):
    pass


def get_parser():
    p = get_default_argparser_with_base_opts(__version__, __doc__)
    f = p.add_argument

    f('subreadset_path', type=validate_file, help="Path to SubreadSet XML")
    f('--name', type=str, default="Auto Job Name", help="Job Name")
    f('--host', type=str, default="smrtlink-beta", help="SMRT Link host")
    f('--port', type=int, default=8081, help="SMRT Link port")
    f('-r', '--reference-set-uuid', type=str, default=Constants.DEFAULT_RSET_UUID, help="ReferenceSet UUID")
    f('-p', '--pipeline-id', type=str, default=Constants.DEFAULT_PIPELINE, help="Pipeline Id")
    f('--block', default=False, action="store_true", help="Block and Wait for job to complete")

    return p


def or_raise_not_found(result, msg):
    if result is None:
        raise ImportJobNotFoundException(msg)
    return result


def get_job_by_subreadset_uuid(sal, uuid):
    """
    Will resolve the first to the most recent Analysis Job or raise AnalysisJobNotFoundException

    The SubreadSet must already be imported and import-dataset job must be
    successfully completed (or will raise ImportJobNotFoundException)

    :type sal: ServiceAccessLayer
    """

    ds = sal.get_subreadset_by_id(uuid)

    import_job_id = or_raise_not_found(ds['jobId'], "Unable to find {}".format(uuid))

    # for the Report Metrics
    import_job = or_raise_not_found(sal.get_job_by_id(import_job_id), "Unable to find Import Job {}".format(import_job_id))

    log.debug("Found import job for SubreadSet")
    log.debug("UUID          : {}".format(uuid))
    log.debug("Import Job Id : {}".format(ds['jobId']))
    log.debug("context       : {}".format(ds['metadataContextId']))
    log.debug("path          : {}".format(ds['path']))

    log.info("Import Job for SubreadSet {}".format(uuid))
    log.info(import_job)

    # This is really brutal. This requires getting *all* the jobs, iterating over them
    # and make a separate call to the entry-points.
    # Assuming the the newest Analysis Job is the expected result (and to potentially)
    # return the quickest result, therefore, reversing the order of the jobs
    all_jobs = sal.get_analysis_jobs()
    all_jobs.reverse()

    log.info("Found {} total analysis jobs".format(len(all_jobs)))
    for job in all_jobs:
        if job.state == JobStates.SUCCESSFUL:
            epoints = sal.get_analysis_job_entry_points(job.id)
            for ep in epoints:
                if ep.dataset_uuid == uuid:
                    return job

    raise AnalysisJobNotFoundException("Unable to find SUCCESSFUL analysis job for SubreadSet UUID {}".format(uuid))


def get_job_by_subreadset_uuid_or_none(sal, sset_uuid):
    try:
        job = get_job_by_subreadset_uuid(sal, sset_uuid)
        return job
    except AnalysisJobNotFoundException:
        return None


def run_main(path, host, port, job_name, pipeline_id, referenceset_uuid, block=False, custom_options=None):
    """

    :param path: Path to SubreadSet XML will be imported (if it's not already been imported)
    :param host: SL Host
    :param port: SL Port
    :param job_name:  Job name
    :param pipeline_id:  Pipeline Id (e.g, pbsmrtpipe.pipelines.my_pipeline
    :param referenceset_uuid: UUID of Rset. This *must* already be imported
    :param block: To block and poll for the analysis job to complete

    :param custom_options: Dictionary of task options for the provided
    Pipeline in the form
    {"pbalign.task_options.concordant":True}


    :type custom_options: dict | None
    :rtype: int
    """

    # look up the reference set UUID from pbservice CLI util or
    # http://smrtlink-beta:8081/secondary-analysis/datasets/references
    # TODO. 1. Import SubreadSet if it's not already imported
    # TODO. 2. Check and see if the Job with the SubreadSet UUID was already submitted
    # TODO. 3. Add option to force a new submission to override (2)
    # TODO. 4. Enable custom pipeline options json file at the CLI

    # sanity test
    sset = SubreadSet(path)
    log.info("Loaded SubreadSet {}".format(sset))

    sal = ServiceAccessLayer(host, port)
    # Sanity Check
    _ = sal.get_status()

    # Step 1. Import SubreadSet (and block) if it's not imported already
    service_sset = sal.get_subreadset_by_id(sset.uuid)
    # TODO. Add check to see if Job was successful
    if service_sset is None:
        log.info("Running Import-DataSet job with {}".format(path))
        sset_import_job = sal.run_import_dataset_subread(path)
        log.info("Import-DataSet job {}".format(sset_import_job))
    else:
        log.info("Found already imported SubreadSet {}".format(service_sset))

    # Step 2. Check and See if an previous analysis job has already been run
    # Immediately exit if an analysis job is found
    analysis_job = get_job_by_subreadset_uuid_or_none(sal, sset.uuid)
    if analysis_job is not None:
        log.info("Found exiting job {} for SubreadSet {}".format(analysis_job, sset))
        return 0

    # Step 3. Create a new Analysis job with custom task options (if provided)
    task_options = {} if custom_options is None else custom_options

    # Get the already Successfully imported DataSets
    service_sset_d = sal.get_dataset_by_uuid(sset.uuid)
    service_rset_d = sal.get_dataset_by_uuid(referenceset_uuid)

    f = sal.run_by_pipeline_template_id if block else sal.create_by_pipeline_template_id

    # The API takes the Int id of the DataSet
    epoints = (ServiceEntryPoint("eid_subread", FileTypes.DS_SUBREADS.file_type_id, service_sset_d['id']),
               ServiceEntryPoint("eid_ref_dataset", FileTypes.DS_REF.file_type_id, service_rset_d['id']))

    job = f(job_name, pipeline_id, epoints, task_options=task_options)

    log.info("Analysis Job {}".format(job))

    if block:
        exit_code = 0 if job.state == JobStates.SUCCESSFUL else 1
    else:
        # the job is in the created state
        exit_code = 0

    return exit_code


def args_runner(args):
    log.debug("Raw args {a}".format(a=args))
    return run_main(args.subreadset_path, args.host, args.port, args.name,
                    args.pipeline_id, args.reference_set_uuid, block=args.block)


def main(argv):
    return pacbio_args_runner(argv[1:],
                              get_parser(),
                              args_runner,
                              log,
                              setup_log_func=setup_log)


if __name__ == '__main__':
    sys.exit(main(sys.argv))