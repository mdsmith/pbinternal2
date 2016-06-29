import os
import logging
import traceback
import tempfile

#from pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.report import rainbow
from pbcore.io import DataSet
from pbcore.data import datasets as data

log = logging.getLogger(__name__)

EXE = 'makeRainbowReport'

_DATA_NAME = 'RainbowReport'
_DATA_DIR = os.path.join(ROOT_DATA_DIR, _DATA_NAME)
_SMALL_INPUT_CMP = os.path.join(PACKAGE_DATA_DIR, 'cmph5_mapping.cmp.h5')

class TestRainbowReport(BaseTestCase):
    """Unit and integration tests for the RainbowReport class"""

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
        self.output = tempfile.mkdtemp(suffix="rainbow_report")

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

    def test_data_import(self):
        ds = DataSet(data.getBam())
        self.assertEqual(ds.numExternalResources, 1)

    def test_constructor(self):
        report = rainbow.make_report(_SMALL_INPUT_CMP, out_dir=self.output)
        self.assertIsNotNone(report)

        report = rainbow.make_report(data.getBam(), out_dir=self.output)
        self.assertIsNotNone(report)

        report = rainbow.make_report(data.getXml(8), out_dir=self.output)
        self.assertIsNotNone(report)

    def test_main(self):
        report = rainbow.make_report(_SMALL_INPUT_CMP, out_dir=self.output)
        self.assertTrue(report)
        self.assertNotEqual(report, 1)

        report = rainbow.make_report(data.getBam(), out_dir=self.output)
        self.assertTrue(report)
        self.assertNotEqual(report, 1)

        report = rainbow.make_report(data.getXml(8), out_dir=self.output)
        self.assertTrue(report)
        self.assertNotEqual(report, 1)

    def test_console_script(self):
        cmd = "{e} {b} --output {o}".format(e=EXE, b=_SMALL_INPUT_CMP,
                                            o=self.output)
        log.debug(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

        cmd = "{e} {b} --output {o}".format(e=EXE, b=data.getBam(),
                                            o=self.output)
        log.info(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

        cmd = "{e} {b} --output {o}".format(e=EXE, b=data.getXml(8),
                                            o=self.output)
        log.info(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

    def test_data_by_reference(self):
        _data = rainbow._data_by_reference(_SMALL_INPUT_CMP)
        expected = [476, 0.9432773109243697]
        self.assertListEqual(list(_data[_data.keys()[0]][0]), expected)

        _data = rainbow._data_by_reference(data.getBam())
        expected = [705, 0.89787234042553188]
        self.assertListEqual(list(_data[_data.keys()[0]][0]), expected)

        _data = rainbow._data_by_reference(data.getXml(8))
        expected = [705, 0.89787234042553188]
        self.assertListEqual(list(_data[_data.keys()[0]][0]), expected)
