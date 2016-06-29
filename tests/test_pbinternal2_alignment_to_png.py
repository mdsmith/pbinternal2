import os
import logging
import traceback
import unittest
import tempfile
import shlex
import numpy as np

#from pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.util import alignment_to_png
from pbcore.data import datasets as data

log = logging.getLogger(__name__)

EXE = 'alignmentToPng'

# TODO:
_DATA_NAME = 'RainbowReport'
_DATA_DIR = os.path.join(ROOT_DATA_DIR, _DATA_NAME)
#_SMALL_INPUT_BAM = os.path.join(PACKAGE_DATA_DIR, 'bam_mapping.bam')
_SMALL_INPUT_CMP = os.path.join(PACKAGE_DATA_DIR, 'cmph5_mapping.cmp.h5')

class TestAlignmentToPng(BaseTestCase):
    """Unit and integration tests for the AlignmentToPng class"""

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
        self.output = tempfile.mkdtemp(suffix="alignment_to_png")

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

    def test_make_report(self):
        report = alignment_to_png.make_report(_SMALL_INPUT_CMP,
                                            out_dir=self.output)
        self.assertIsNotNone(report)

        report = alignment_to_png.make_report(data.getBam(),
                                              out_dir=self.output)
        self.assertIsNotNone(report)

    def test_main(self):
        cmd = "{e} {b} --output {o}".format(e=EXE, b=_SMALL_INPUT_CMP,
                                            o=self.output)
        log.info(cmd)
        rcode = alignment_to_png.main(shlex.split(cmd))
        self.assertEqual(rcode, 0)

        cmd = "{e} {b} --output {o}".format(e=EXE, b=data.getBam(),
                                            o=self.output)
        log.info(cmd)
        rcode = alignment_to_png.main(shlex.split(cmd))
        self.assertEqual(rcode, 0)

    def test_console_script(self):
        cmd = "{e} {b} --output {o}".format(e=EXE, b=_SMALL_INPUT_CMP,
                                            o=self.output)
        log.info(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)

    def test_readInFile(self):
        from pbreports.plot.rainbow import _read_in_file
        _data = _read_in_file(_SMALL_INPUT_CMP)
        expected = np.array([4.76000000e+02, 0.9432773109243697,
                             2.54000000e+02])
        self.assertListEqual(_data[0, :].tolist(), expected.tolist())

        _data = _read_in_file(data.getBam())
        expected = np.array([1551.0,
                             0.8852353320438426,
                             41.0])
        self.assertListEqual(_data[0, :].tolist(), expected.tolist())
