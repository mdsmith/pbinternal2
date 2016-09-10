
"""
Tests for end-to-end tool contract support.  Any report tool that has emit and
resolve functions should have a class here.
"""

import unittest
import os.path
import shutil
import pbtestdata
import pbcommand.testkit

class TestEolQcReport(pbcommand.testkit.PbTestApp):
    MAX_NPROC = 12
    RESOLVED_NPROC = 1
    SSET = pbtestdata.get_file('internal-subreads')
    ASET = pbtestdata.get_file('aligned-internal-subreads')
    DRIVER_BASE = "python -m pbinternal2.tasks.eol_qc "
    DRIVER_EMIT = (DRIVER_BASE +
                   " emit-tool-contract pbinternal2.tasks.eol_qc")
    DRIVER_RESOLVE = DRIVER_BASE + " run-rtc "
    REQUIRES_PBCORE = True
    INPUT_FILES = [SSET, ASET]
    TASK_OPTIONS = {}

    def run_after(self, rtc, output_dir):
        shutil.rmtree(output_dir)    

if __name__ == "__main__":
    unittest.main()
