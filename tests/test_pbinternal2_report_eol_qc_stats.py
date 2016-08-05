import os
import logging
import traceback
import tempfile

#from pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.report.eol_qc_stats import loading_efficiency
from pbcore.io import SubreadSet, AlignmentSet
from pbcore.data import datasets as data

log = logging.getLogger(__name__)

class TestEOLQCStats(BaseTestCase):

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
        self.output = tempfile.mkdtemp(suffix="eolqcstats")

    def tearDown(self):
        """
        After *every* test
        """
        try:
            BaseTestCase.tearDown(self)
        except Exception as err:
            log.exception("EOLQCStats setUp failed.")
            raise


    def test_loading_efficiency(self):
        sset = ('/pbi/collections/312/3120170/r54009_20160728_220413/'
                '1_A01/m54009_160728_220732.subreadset.xml')
        aset = ('/pbi/dept/secondary/siv/smrtlink/smrtlink-beta/'
                'smrtsuite_166987/userdata/jobs_root/018/018404/'
                'tasks/pbcoretools.tasks.gather_alignmentset-1/'
                'file.alignmentset.xml')
        sset = SubreadSet(sset)
        aset = AlignmentSet(aset)
        exp = 84.1558498903621
        obs = loading_efficiency(aset)
        self.assertEqual(obs, exp)
