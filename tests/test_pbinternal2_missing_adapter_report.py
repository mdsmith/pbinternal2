import os
import logging
import traceback
import unittest
from unittest.case import SkipTest
import tempfile
import shlex

#from pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.report.missing_adapter import (make_report, _num_gt_min,
                                                _get_lengths, main,
                                                _make_histogram_png_theme_lpr)

log = logging.getLogger(__name__)

EXE = 'makeMissingAdapterReport'
_INITIAL_FASTA = os.path.join(PACKAGE_DATA_DIR,
                              'missing_adapter_filtered_subreads.fasta')
_INITIAL_BAM = os.path.join(PACKAGE_DATA_DIR,
                            'missing_adapter_filtered_subreads.bam')
_INITIAL_XML = os.path.join(PACKAGE_DATA_DIR,
                            'missing_adapter_filtered_subreads.xml')
_MISSING_FASTA = os.path.join(PACKAGE_DATA_DIR,
                              'missing_adapter_Missing.fasta')
_MISSING_BAM = os.path.join(PACKAGE_DATA_DIR,
                            'missing_adapter_Missing.bam')
_MISSING_XML = os.path.join(PACKAGE_DATA_DIR,
                            'missing_adapter_Missing.xml')

class TestMissingAdapterReport(BaseTestCase):
    """Unit and integration tests for the Missing Adapter Report class and \
    associated module functions"""

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
        self.output = tempfile.mkdtemp(suffix="missing_adapter_reports")

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

    def test_list_input(self):
        report = make_report(_INITIAL_FASTA, [_MISSING_FASTA],
                                           self.output)

    def test_main(self):
        cmd = "{e} --output {o} {i} {m}".format(e=EXE, i=_INITIAL_FASTA,
                                                m=_MISSING_FASTA,
                                                o=self.output)
        log.info(cmd)
        rcode = main(shlex.split(cmd))
        self.assertEqual(rcode, 0)

    def test_console_script(self):
        """Smoke test 2xfasta"""
        cmd = "{e} --output {o} {i} {m}".format(e=EXE, i=_INITIAL_FASTA,
                                                m=_MISSING_FASTA,
                                                o=self.output)
        log.info(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

    #Now that pbalign won't take fasta files, turning off until a sample bam is
    #available
    @SkipTest
    def test_bam_io(self):
        """Smoke test 2xbam"""
        cmd = "{e} --output {o} {i} {m}".format(e=EXE, i=_INITIAL_BAM,
                                                m=_MISSING_BAM,
                                                o=self.output)
        log.info(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

    #Now that pbalign won't take fasta files, turning off until a sample bam is
    #available
    @SkipTest
    def test_xml_io(self):
        """Smoke test 2xXML"""
        cmd = "{e} --output {o} {i} {m}".format(e=EXE, i=_INITIAL_XML,
                                                m=_MISSING_XML,
                                                o=self.output)
        log.info(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

    def test_get_lengths(self):
        def checkList(list1, list2):
            return len(list1) == len(list2) and sorted(list1) == sorted(list2)
        self.assertTrue(
            checkList([5984, 6570, 8218, 6455, 11151, 1891, 740],
                      _get_lengths(_INITIAL_FASTA)))

    def test_make_histogram_png_theme_lpr(self):
        # This is testing a very simple composition of pbreports functions that
        # should already be tested in pbreports
        data = list(range(0, 1000)) + list(range(100, 300))
        png_file_name = os.path.join(self.output, "test.png")
        axes_names = ("Number", "Lengths")
        _make_histogram_png_theme_lpr(data, png_file_name, axes_names)
        self.assertTrue(os.path.exists(png_file_name))

