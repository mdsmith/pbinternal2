import os
import logging
import traceback
import unittest
import tempfile

#from pbcore.util.Process import backticks
from base_test_case import PACKAGE_DATA_DIR, ROOT_DATA_DIR, run_backticks
from base_test_case import BaseTestCase
from pbinternal2.util.range import Range, Ranges, OverlappingRanges

log = logging.getLogger(__name__)

class TestPBRange(BaseTestCase):
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
        self.output = tempfile.mkdtemp(suffix="range_tests")

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

    def test_range(self):
        r1 = Range(5, 10)
        r2 = Range(10, 15)
        r3 = Range(7, 15)
        r4 = Range(5, 10)
        self.assertTrue(r1.contains(5))
        self.assertFalse(r1.contains(10))
        self.assertFalse(r1.intersects(r2))
        self.assertTrue(r1.intersects(r3))
        self.assertTrue(Range(0, 0) == r1.intersect(r2))
        self.assertTrue(Range(7, 10) == r1.intersect(r3))
        self.assertTrue(r1 == r4)

    def test_ranges(self):
        r = Ranges()
        r.add_range(Range(1, 3))
        r.add_range(Range(9, 12))
        self.assertEqual("Rs {[1,3) [9,12)}", str(r))
        r.add_range(Range(2, 5))
        self.assertEqual("Rs {[1,5) [9,12)}", str(r))
        r.add_range(Range(1, 12))
        self.assertEqual("Rs {[1,12)}", str(r))
        r.add_range(Range(14, 15))
        self.assertEqual("Rs {[1,12) [14,15)}", str(r))
        r.add_range(Range(20, 25))
        self.assertEqual("Rs {[1,12) [14,15) [20,25)}", str(r))
        r.add_range(Range(11, 22))
        self.assertEqual("Rs {[1,25)}", str(r))

        r1 = Range(5, 10)
        r.remove_range(r1)
        self.assertEqual("Rs {[1,5) [10,25)}", str(r))
        r5 = Range(3, 12)
        r.remove_range(r5)
        self.assertEqual("Rs {[1,3) [12,25)}", str(r))
        r.remove_range(Range(-1, 3))
        self.assertEqual("Rs {[12,25)}", str(r))
        r.remove_range(Range(1, 3))
        self.assertEqual("Rs {[12,25)}", str(r))
        r.remove_range(Range(1, 25))
        self.assertEqual("Rs {}", str(r))

        r5 = Ranges()
        r5.add_range(Range(1, 25))
        r5.add_range(Range(27, 29))
        r5.add_range(Range(35, 40))
        r6 = Ranges()
        r6.add_range(Range(2, 5))
        r6.add_range(Range(20, 30))
        r6.add_range(Range(42, 45))
        r5.merge(r6)
        self.assertEqual("Rs {[1,30) [35,40) [42,45)}", str(r5))

    def test_overlapping_ranges(self):
        r = OverlappingRanges()
        r1 = Range(0, 15)
        r2 = Range(1, 3)
        r3 = Range(9, 12)
        r.add_range(r1)
        r.add_range(r2)
        r.add_range(r3)
        self.assertEqual("Rs {[0,15) [1,3) [9,12)}", str(r))

        query_ranges = [Range(0, 1), Range(2, 4), Range(13, 15), Range(15, 17)]
        answers = [[r1, r2], [r1, r2], [r1], []]
        for q_range, expected in zip(query_ranges, answers):
            for actual, exp in zip(r.overlapping_ranges(q_range), expected):
                self.assertTrue(actual == exp)
