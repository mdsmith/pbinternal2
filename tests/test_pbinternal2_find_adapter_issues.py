
import os
import logging
import traceback
import itertools
import tempfile
import csv

from nose.plugins.skip import SkipTest
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.util.find_adapter_issues import find_adapter_issues

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

EXE = 'findAdapterIssues'
_DATA_NAME = 'FindAdapterIssues'
_DATA_DIR = os.path.join(ROOT_DATA_DIR, _DATA_NAME)
_INPUT_FASTA = os.path.join(_DATA_DIR, 'filtered_subreads.fasta')
_INPUT_REFERENCE = os.path.join(_DATA_DIR, 'ecoliK12_pbi_March2013.fasta')
_OUTPUT_CSV = os.path.join(_DATA_DIR, 'adapter_issues.csv')
_SMALL_INPUT_FASTA = os.path.join(PACKAGE_DATA_DIR, 'find_adapter_issues',
                                  'lambda',
                                  'filtered_subreads.small.fasta')
_SMALL_INPUT_REFERENCE = os.path.join(PACKAGE_DATA_DIR, 'find_adapter_issues',
                                      'lambda',
                                      'lambdaNEB.fasta')
_SMALL_INPUT_XML = os.path.join(_DATA_DIR, 'xmlTest', 'reads.xml')
_SMALL_INPUT_BAM = os.path.join(_DATA_DIR, 'xmlTest', 'adapter.subreads.bam')
_SMALL_INPUT_REFERENCE_XML = os.path.join(_DATA_DIR,
                                          'xmlTest',
                                          'ref.xml')
_SMALL_OUTPUT_CSV = os.path.join(PACKAGE_DATA_DIR, 'find_adapter_issues',
                                 'lambda',
                                 'small_issues.csv')
_SMALL_XML_OUTPUT_CSV = os.path.join(_DATA_DIR,
                                     'xmlTest',
                                     'adapter_issues.csv')

class TestFindAdapterIssues(BaseTestCase):
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
        self.output = tempfile.mkdtemp(suffix="test_find_adapter_issues")

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

    # Get appropriate bam source file and turn back on....
    @SkipTest
    def test_on_small_real_data(self):
        outfile = os.path.join(self.output, "test_adapter_issues.csv")
        cmd = "{e} {f} {r} --output_file {o}".format(e=EXE,
                                                     f=_SMALL_INPUT_FASTA,
                                                     r=_SMALL_INPUT_REFERENCE,
                                                     o=outfile)
        log.debug(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)
        with open(outfile, 'r') as test_file:
            with open(_SMALL_OUTPUT_CSV, 'r') as golden:
                test_lines = test_file.readlines()
                golden_lines = golden.readlines()
                for test_l in test_lines:
                    found = False
                    for gold_l in golden_lines:
                        if test_l.split(',')[0] == gold_l.split(',')[0]:
                            found = True
                            if (test_l.split(',')[1] !=
                                    gold_l.split(',')[1]):
                                print test_l
                            self.assertEqual(test_l.split(',')[1],
                                             gold_l.split(',')[1])
                    if found:
                        golden_lines.remove(test_l)
                    self.assertTrue(found)
                # The reverse of this test: all lines in golden should have
                # been removed (as they are found)
                self.assertEqual(len(golden_lines), 0)

    # Skip until proper bam version is available (fasta mapping?)
    @SkipTest
    def test_on_small_xml_data(self):
        outfile = os.path.join(self.output, "test_adapter_issues.csv")
        cmd = "{e} {f} {r} --output_file {o}".format(
            e=EXE, f=_SMALL_INPUT_XML, r=_SMALL_INPUT_REFERENCE_XML, o=outfile)
        log.debug(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)
        with open(outfile, 'r') as test_file:
            with open(_SMALL_XML_OUTPUT_CSV, 'r') as golden:
                notSame = 0
                test_lines = test_file.readlines()
                golden_lines = golden.readlines()
                for test_l in test_lines:
                    found = False
                    for gold_l in golden_lines:
                        if test_l.split(',')[0] == gold_l.split(',')[0]:
                            found = True
                            if (test_l.split(',')[1] !=
                                    gold_l.split(',')[1]):
                                notSame += 1
                    if found:
                        try:
                            golden_lines.remove(test_l)
                        except ValueError:
                            pass
                    self.assertTrue(found)
                    self.assertTrue(notSame/len(test_lines) < 0.025)
                # The reverse of this test: all lines in golden should have
                # been removed (as they are found)
                self.assertEqual(len(golden_lines), notSame)

    @SkipTest
    def test_on_small_bam_data(self):
        outfile = os.path.join(self.output, "test_adapter_issues.csv")
        cmd = "{e} {f} {r} --output_file {o}".format(
            e=EXE, f=_SMALL_INPUT_BAM, r=_SMALL_INPUT_REFERENCE_XML, o=outfile)
        log.debug(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)
        with open(outfile, 'r') as test_file:
            with open(_SMALL_XML_OUTPUT_CSV, 'r') as golden:
                notSame = 0
                test_lines = test_file.readlines()
                golden_lines = golden.readlines()
                for test_l in test_lines:
                    found = False
                    for gold_l in golden_lines:
                        if test_l.split(',')[0] == gold_l.split(',')[0]:
                            found = True
                            if (test_l.split(',')[1] !=
                                    gold_l.split(',')[1]):
                                notSame += 1
                    if found:
                        try:
                            golden_lines.remove(test_l)
                        except ValueError:
                            pass
                    self.assertTrue(found)
                    self.assertTrue(notSame/len(test_lines) < 0.025)
                # The reverse of this test: all lines in golden should have
                # been removed (as they are found)
                self.assertEqual(len(golden_lines), notSame)

    # Skip until proper bam version is available (fasta mapping?)
    @SkipTest
    def test_on_large_real_data(self):
        outfile = os.path.join(self.output, "test_adapter_issues.csv")
        cmd = "{e} {f} {r} --output_file {o} --nproc 14".format(
            e=EXE, f=_INPUT_FASTA, r=_INPUT_REFERENCE, o=outfile)
        log.debug(cmd)
        rcode = run_backticks(cmd)
        self.assertEqual(rcode, 0)
        differences = []
        with open(outfile, 'r') as test_file:
            with open(_OUTPUT_CSV, 'r') as golden:
                golden_lines = golden.readlines()
                #num_found = 0
                for test_l in test_file:
                    found = False
                    for gold_l in golden_lines:
                        if test_l.split(',')[0] == gold_l.split(',')[0]:
                            #num_found += 1
                            found = gold_l
                            if test_l.split(',')[1] != gold_l.split(',')[1]:
                                differences.append(test_l + " " + gold_l)
                            #self.assertEqual(test_l.split(',')[1],
                                             #gold_l.split(',')[1])
                    #if (not found and
                            #test_l.split(',')[1].strip() != "AlignedNoIssue"):
                        #print test_l
                    if found:
                        golden_lines.remove(found)
                    #self.assertTrue(found)
                #print "Num found: {i}".format(i=num_found)
                for difference in differences:
                    print difference
                self.assertEqual(len(golden_lines), 0)
                #for test_l, gold_l in itertools.izip(test_file, golden):
                    #self.assertEqual(test_l, gold_l)
