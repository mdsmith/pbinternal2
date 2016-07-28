
"""
Tests for end-to-end tool contract support.  Any report tool that has emit and
resolve functions should have a class here.
"""

import unittest
import os.path
import pbtestdata
import pbcommand.testkit

class TestEolQcReport(pbcommand.testkit.PbTestApp):
    MAX_NPROC = 12
    RESOLVED_NPROC = 1
    SSET = pbtestdata.get_file('internal-subreads')
    ASET = pbtestdata.get_file('aligned-internal-subreads')
    DRIVER_BASE = "python -m pbinternal2.tasks.eol_qc "
    DRIVER_EMIT = (DRIVER_BASE +
                   " emit-tool-contract pbsmrtpipe_internal.tasks.eol_qc")
    DRIVER_RESOLVE = DRIVER_BASE + " run-rtc "
    REQUIRES_PBCORE = True
    INPUT_FILES = [SSET, ASET]
    TASK_OPTIONS = {}

if __name__ == "__main__":
    unittest.main()
