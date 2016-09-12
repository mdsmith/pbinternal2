
"""
Tests for end-to-end tool contract support.  Any report tool that has emit and
resolve functions should have a class here.
"""

import unittest
import os.path
import pbtestdata
import pbcommand.testkit
from pbcommand.models import ReseqCondition, ReseqConditions
import json
import shutil
import tempfile

class TestReseqCondReport(pbcommand.testkit.PbTestApp):
    MAX_NPROC = 12
    RESOLVED_NPROC = 1
    # Make a temporary file with the condition summary
    # In the future if we take a Java dependency we can make this directly with
    # See: https://github.com/PacificBiosciences/smrtflow/tree/master/smrt-server-analysis-internal#resequencing-condition-type 
    SSET = pbtestdata.get_file('internal-subreads')
    ASET = pbtestdata.get_file('aligned-internal-subreads')
    RSET = pbtestdata.get_file('lambdaNEB')
    rc = ReseqCondition("TestCond", SSET, ASET, RSET)
    rcs = ReseqConditions([rc])
    as_dict = rcs.to_dict()
    as_dict["_condition_doc"] = "Example of a 'Resequencing' Condition Type"
    json_out = json.dumps(as_dict)
    # We need to delete this when done.
    hndl, tmpname = tempfile.mkstemp()
    f = os.fdopen(hndl, "w")
    f.write(json_out)
    f.close()
    INPUT_FILES = [tmpname]
    
    DRIVER_BASE = "python -m pbinternal2.analysis_tools "
    DRIVER_EMIT = (DRIVER_BASE +
                   " emit-tool-contract pbinternal2.tasks.cond_to_report")
    DRIVER_RESOLVE = DRIVER_BASE + " run-rtc "
    REQUIRES_PBCORE = True

    TASK_OPTIONS = {}
    
    def run_after(self, rtc, output_dir):
        for f in self.INPUT_FILES:
            os.remove(f)
        shutil.rmtree(output_dir)

if __name__ == "__main__":
    unittest.main()
