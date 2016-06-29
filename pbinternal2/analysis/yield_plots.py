import logging


from pbcore.io.dataset.DataSetIO import AlignmentSet
from pbcore.io import openIndexedAlignmentFile
from pbreports.plot.helper import save_figure_with_thumbnail

from pbcommand.models.report import Report, PlotGroup, Plot
from pbcommand.pb_io import load_reseq_conditions_from
from pbcommand.models import ReseqConditions

# FIXME
# import pandas as pd

log = logging.getLogger(__name__)
log.addHandler(logging.StreamHandler())

log.setLevel(logging.DEBUG)


def run_yield_plots(cond_json_path, output_report):
    """ Run analyses that generate plots and return fig objects for
    subsequent report
    :param cond_json_path:
    :param output_report:
    :return:
    """
    # FIXME
    import pandas as pd

    log.info("Running yield plot analysis")
    plots = []
    c = load_reseq_conditions_from(cond_json_path)
    log.info("Generating alignment yield plot")
    yield_plot_alignment = run_alignment_yield_plot(c)
    plots.append(('yield_plot_alignment', yield_plot_alignment))

    plot_list = []
    for plot_id, figure in plots:
        plot_path = '{i}.png'.format(i=plot_id)
        thumb_path = '{i}_thumbs.png'.format(i=plot_id)
        save_figure_with_thumbnail(figure, plot_path)
        plot = Plot(plot_id, plot_path, thumbnail=thumb_path)
        plot_list.append(plot)

    # Format the output for pbreports style Report.json
    plot_group = PlotGroup('yield_plots', plots=plot_list)
    report = Report('yield_plots', tables=(), plotgroups=[plot_group],
                    attributes=())
    log.info("Writing report: {f}".format(f=output_report))
    report.write_json(output_report)
    log.info("Completed yield plot analysis")
    return 0


def _get_colors(num_conditions):
    colorlist = ['cyan', 'magenta', 'yellow', 'fuchsia', 'k']
    return colorlist[:num_conditions]


def run_alignment_yield_plot(analysis_conditions):
    """
    :type analysis_conditions: ReseqConditions

    :param output_report: path to output_report.json
    :return: Returns matplotlib fig object
    """
    import pandas as pd

    conditions = analysis_conditions.conditions
    log.info("Generating Plot for conditions")
    log.info(analysis_conditions)

    c_dict = {}
    num_conditions = len(conditions)
    colorlist = _get_colors(num_conditions)

    for i, val in enumerate(conditions):
        dataset_path = conditions[i].files
        id = conditions[i].id

        # eventually we will might want to merge bamfiles from the same
        # condition so this will need to handle that
        dataset = dataset_path[0]

        alignment_dataset = AlignmentSet(dataset)
        alignment_bam = alignment_dataset.toExternalFiles()[0]
        indexed_alignments = openIndexedAlignmentFile(alignment_bam)
        pbi = indexed_alignments.pbi
        c_dict[id] = pd.Series(len(pbi), index=None)

    log.debug("Loading DataFrame")
    df = pd.DataFrame(c_dict)

    log.debug("Plotting DataFrame")
    ax = df.plot(kind="bar",
                 title="Yield (Number of Alignments)",
                 legend=True,
                 color=colorlist)

    log.debug("Adjusting Plot")
    # move legend outside of box
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.set_ylabel("Number of Alignments")
    ax.yaxis.grid()
    ax.legend(loc="center left",
              bbox_to_anchor=(1, 0.5),
              title="Conditions",
              frameon=False)
    ax.tick_params(labelbottom='off')
    fig = ax.get_figure()

    return fig
