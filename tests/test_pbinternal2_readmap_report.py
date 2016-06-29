import os
import logging
import traceback
import unittest
import tempfile

#from pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.report import readmap
from pbcore.data import datasets as data

log = logging.getLogger(__name__)

EXE = 'makeReadMapReport'

class TestReadMapReport(BaseTestCase):
    """Unit and integrationt tests for the Missing Adapter Report class and \
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
        self.output = tempfile.mkdtemp(suffix="read_map_reports")

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

    def test_constructor(self):
        log.info("Writing output files to {s}".format(s=self.output))
        outPngFn = os.path.join(self.output, 'construTest.png')
        ret_code = readmap.make_report(data.getBam(), outPngFn)
        self.assertEqual(ret_code, 0)

    def test_console_script(self):
        """Smoke test 1xbam in"""
        outPngFn = os.path.join(self.output, 'consoleTest.png')
        args = dict(e=EXE, inbam=data.getBam(), out=outPngFn)
        cmd_str = '{e} {inbam} --outfn {out}'.format(**args)

        rcode = run_backticks(cmd_str)
        self.assertEqual(
            rcode, 0,
            "Exit code '{o}' for command '{c}'".format(o=rcode, c=cmd_str))

        outPngFn = os.path.join(self.output, 'consoleTest.png')
        args = dict(e=EXE, inbam=data.getXml(8), out=outPngFn)
        cmd_str = '{e} {inbam} --outfn {out}'.format(**args)

        rcode = run_backticks(cmd_str)
        self.assertEqual(
            rcode, 0,
            "Exit code '{o}' for command '{c}'".format(o=rcode, c=cmd_str))

