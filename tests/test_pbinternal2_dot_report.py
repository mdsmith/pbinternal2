import os
import logging
import traceback
import unittest
from unittest.case import SkipTest
import tempfile

#import pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.report import dot

log = logging.getLogger(__name__)

EXE = 'makeDotReport'
_DATA_NAME = 'DotReport'
_DATA_DIR = os.path.join(ROOT_DATA_DIR, _DATA_NAME)
_SAMPLE_FASTA = os.path.join(_DATA_DIR, 'head_ecoli.fasta')
_SMALL_CHIMERA = os.path.join(PACKAGE_DATA_DIR,
                              'dot_report_Chimera_small.fasta')
_SMALL_CHIMERA_BAM = os.path.join(PACKAGE_DATA_DIR,
                                  'dot_report_Chimera_small.bam')
_SMALL_CHIMERA_XML = os.path.join(PACKAGE_DATA_DIR,
                                  'dot_report_Chimera_small.xml')
_LAMBDA_REFERENCE = os.path.join(PACKAGE_DATA_DIR, 'lambdaNEB.fasta')
_LAMBDA_REFERENCE_XML = os.path.join(PACKAGE_DATA_DIR, 'lambdaNEB.xml')


class TestDotReport(BaseTestCase):

    def setUp(self):
        """
        Before *every* test
        """
        try:
            BaseTestCase.setUp(self)
        except Exception as err:
            log.error(err)
            tb = traceback.format_exc()
            log.error(tb)
            raise
        log.debug("In setUp()")
        self.output = tempfile.mkdtemp(suffix="dot_reports")

    def tearDown(self):
        """
        After *every* test
        """
        try:
            BaseTestCase.tearDown(self)
        except Exception as err:
            log.error(err)
            tb = traceback.format_exc()
            log.error(tb)
            raise

    def simple_testcase(self):
        """Test console entry point on 2xfasta in"""
        args = dict(e=EXE, reads=_SMALL_CHIMERA, ref=_LAMBDA_REFERENCE,
                    out=self.output)
        cmd_str = '{e} {reads} {ref} --output {out}'.format(**args)

        rcode = run_backticks(cmd_str)
        self.assertEqual(
            rcode, 0,
            "Exit code '{o}' for command '{c}'".format(o=rcode, c=cmd_str))

    #Now that pbalign won't take fasta files, turning off until a sample bam is
    #available
    @SkipTest
    def test_bam_io(self):
        """Test console entry point on 1xbam 1xfasta in"""
        args = dict(e=EXE, reads=_SMALL_CHIMERA_BAM, ref=_LAMBDA_REFERENCE,
                    out=self.output)
        cmd_str = '{e} {reads} {ref} --output {out}'.format(**args)

        rcode = run_backticks(cmd_str)
        self.assertEqual(
            rcode, 0,
            "Exit code '{o}' for command '{c}'".format(o=rcode, c=cmd_str))

    #Now that pbalign won't take fasta files, turning off until a sample bam is
    #available
    @SkipTest
    def test_xml_io(self):
        """Test console entry point on 2xXML in"""
        args = dict(e=EXE, reads=_SMALL_CHIMERA_XML, ref=_LAMBDA_REFERENCE_XML,
                    out=self.output)
        cmd_str = '{e} {reads} {ref} --output {out}'.format(**args)

        rcode = run_backticks(cmd_str)
        self.assertEqual(
            rcode, 0,
            "Exit code '{o}' for command '{c}'".format(o=rcode, c=cmd_str))

    def test_report(self):
        args = dict(e=EXE, reads=_SMALL_CHIMERA, ref=_LAMBDA_REFERENCE,
                    out=self.output)

        report = dot.make_report(args["reads"], args["ref"], args["out"])
        self.assertTrue(report)

