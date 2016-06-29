import logging
import os
import traceback
import tempfile
from base_test_case import BaseTestCase
from pbinternal2.analysis.yield_plots import run_yield_plots, run_alignment_yield_plot
from pbinternal2.models import AnalysisConditions
from matplotlib.figure import Figure as fig
from pbinternal2.util.range import Range, Ranges, OverlappingRanges
from base_test_case import get_temp_file, get_pkg_data_file

log = logging.getLogger(__name__)

class TestYieldPlots(BaseTestCase):
    """Unit and integration tests for analysis yield_plots"""

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
        self.output = tempfile.mkdtemp(suffix="yield_plots")

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
    def test_run_e2e_yield_plots(self):
        conditions = get_pkg_data_file("conditions_hello_world.json")
        output = os.path.join(self.output, "yield_plot_report.json")
        rcode = run_yield_plots(conditions, output)
        self.assertIs(rcode, 0)

    def test_yield_plots(self):
        conditions = get_pkg_data_file("conditions_hello_world.json")
        c = AnalysisConditions.load_conditions_from(conditions)
        figure = run_alignment_yield_plot(c)
        self.assertIsInstance(type(figure), type(fig))