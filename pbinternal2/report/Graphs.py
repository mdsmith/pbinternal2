
"""A set of stylistically and behaviorally consisteng graphing functions from
pricompare"""

import logging
from collections import defaultdict, OrderedDict
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import matplotlib.patches as pat
from matplotlib import cm
import numpy as np
from scipy.stats import gaussian_kde
from Bio.Statistics.lowess import lowess
from pbcore.io.dataset.utils import hn_to_xy

log = logging.getLogger(__name__)

FIG_SIZE = (18, 10.75)

def lwrt(fname):
    log.info("Writing {:}".format(fname))

def savefig(fig, fname):
    legend = fig.get_axes()[0].get_legend()
    if not legend is None:
        fig.savefig(fname, bbox_extra_artists=[legend],
                    bbox_inches='tight')
    else:
        fig.savefig(fname)
    lwrt(fname)

def swap_extension(fname, new_ext):
    if new_ext.startswith('.') and len(new_ext) > 1:
        new_ext = new_ext[1:]
    num_pop = 1
    if fname.endswith('h5') or fname.endswith('xml') or fname.endswith('bam'):
        num_pop = 2
    tbr = fname.split('.')[:(-1 * num_pop)]
    if new_ext:
        tbr.append(new_ext)
    tbr = '.'.join(tbr)
    return tbr

def check_figax(fig, ax, figsize=FIG_SIZE):
    if fig is None or ax is None:
        fig = plt.figure(figsize=figsize, dpi=200)
        ax = fig.add_subplot(111)
    return fig, ax

def draw_chip(data, fig=None, ax=None, rowstart=64, colstart=64, rows=1024,
              cols=1144, figsize=FIG_SIZE, color='green', drawfull=True):
    """data is a list of hns to highlight on the chip"""
    fig, ax = check_figax(fig, ax, figsize)
    x = []
    y = []
    if len(data):
        for hn in data:
            if not isinstance(hn, tuple):
                hn = hn_to_xy(hn)
            x.append(hn[0])
            y.append(hn[1])
        fig, ax = gen_mh_scatter(x, y, fig=fig, ax=ax, color=color)
    if drawfull:
        margin = 10
        ax.set_xlim(xmin=colstart - 10, xmax=cols + 10)
        ax.set_ylim(ymin=rows - 10, ymax=rowstart + 10)
    #else:
    #    ax.set_ylim(ax.get_ylim()[::-1])
    return fig, ax

def mh_conf_matrix(matrix, title="Error Rate Matrix", norm=True, mod_diag=True,
                   even_range=False, vmin=-0.1, vmax=0.1):
    if mod_diag:
        diag_f = lambda x: x * -1.0
    else:
        diag_f = None
    return gen_mh_heatmap(matrix, list("ACTG-#"), list("ACTG-"),
                          "Observed", "Ground Truth", rowsum=False,
                          logcolor=False, title=title,
                          diag_f=diag_f,
                          norm=norm, even_range=even_range, vmin=vmin,
                          vmax=vmax)

def gen_mh_heatmap(data, xlabs=None, ylabs=None, xtitle=None, ytitle=None,
                   norm=True, rowsum=False,
                   logcolor=False, title=None, diag_f=None, even_range=True,
                   vmin=-0.1, vmax=0.1, radius=None, fig=None, ax=None,
                   figsize=FIG_SIZE):
    fig, ax = check_figax(fig, ax, figsize)
    # override for now:
    vmin = -0.1
    vmax = 0.1
    radius = None

    fontsize = 14
    rowsums = []
    data = np.array(data, dtype=float)
    if norm:
        if rowsum:
            for i, row in enumerate(data):
                tot = np.sum(row)
                rowsums.append(tot)
        tot = np.sum(data)
        data = np.true_divide(data, tot)

    if rowsum:
        rowsums = np.array(rowsums).reshape((len(rowsums), 1))
        fig, (ax, rsax) = plt.subplots(nrows=1, ncols=2,
                                       gridspec_kw = {'width_ratios':[5, 1]})

        rsax.set_frame_on(False)

        heatmap = rsax.pcolor(rowsums, cmap=cm.get_cmap('Reds'))

        for t in rsax.xaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False
        for t in rsax.yaxis.get_major_ticks():
            t.tick1On = False
            t.tick2On = False

        xticks = np.arange(rowsums.shape[1])+0.5
        yticks = np.arange(rowsums.shape[0])+0.5
        summax = np.max(rowsums)
        for y, row in zip(yticks, rowsums):
            for x, val in zip(xticks, row):
                rsax.text(x, y + 0.05, '%.0f' % val, ha='center',
                          color='black' if val/summax < 0.7 else 'white',
                          fontsize=fontsize)
        rsax.set_xticks([])
        rsax.set_yticks([])

        rsax.invert_yaxis()

        rsax.set_xticklabels([])
        rsax.set_yticklabels([])
        rsax.set_xlabel('Row Sum', fontsize=fontsize)
    else:
        fig, ax = plt.subplots()

    fig.set_size_inches(18, 10.75)
    fig.set_dpi(200)
    if title:
        fig.suptitle(title, fontsize=20, y=0.95)

    ax.set_frame_on(False)

    color_data = np.copy(data)
    if not diag_f is None:
        # exclude the last row, which should have 0 density
        for i, row in enumerate(color_data[:-1]):
            row[i] = diag_f(row[i])

    # narrow colormap range

    # The range I like of [0,1]:
    cmin = 0.3
    cmax = 0.7
    if not radius is None:
        vmin = -1 * radius
        vmax = 1 * radius
        lvTmp = np.linspace(cmin, cmax, 256)
        cmTmp = cm.get_cmap('seismic')(lvTmp)
        newCmap = mcol.ListedColormap(cmTmp)
    else:
        cmin = 0.45

        # The range midpoint I want at 0:
        cmid = 0.5

        # The fraction of data range below 0
        vmid = 1 - vmax/(vmax + abs(vmin))

        # normalize the range of colors around 0.0, so that small negative
        # amplitudes aren't given deep blues just because they're the only
        # negatives in town.
        if even_range:
            max_range = max(abs(vmax), abs(vmin))
            below_short = 1.0 - abs(vmin)/max_range
            above_short = 1.0 - abs(vmax)/max_range
            cmin += (cmid - cmin) * below_short
            cmax -= (cmax - cmid) * above_short
        cmap = cm.get_cmap('seismic')
        name = 'shifted'

        # realign data '0' normalization with cmap '0'
        cold = defaultdict(list)

        regi = np.hstack([np.linspace(cmin, 0.5, 128, endpoint=False),
                               np.linspace(0.5, cmax, 129)])

        shifti= np.hstack([np.linspace(0.0, vmid, 128, endpoint=False),
                           np.linspace(vmid, 1.0, 129)])

        for ri, si in zip(regi, shifti):
            r, g, b, a = cmap(ri)
            cold['red'].append((si, r, r))
            cold['green'].append((si, g, g))
            cold['blue'].append((si, b, b))
            cold['alpha'].append((si, a, a))

        newCmap = mcol.LinearSegmentedColormap(name, cold)
    cmap = newCmap



    kwargs = {}
    kwargs['cmap'] = cmap
    kwargs['vmin'] = vmin
    kwargs['vmax'] = vmax
    if logcolor:
        kwargs['norm'] = mcol.LogNorm()
    heatmap = ax.pcolor(color_data, **kwargs)

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    norm = mcol.Normalize(vmin=vmin, vmax=vmax)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    xticks = np.arange(data.shape[0])+0.5
    yticks = np.arange(data.shape[1])+0.5
    for y, crow, row in zip(yticks, color_data, data):
        for x, cval, val in zip(xticks, crow, row):
            if logcolor:
                ax.text(x, (y * 1.01) + 0.04, '%.3f' % val, ha='center',
                        color=('black' if sum(mapper.to_rgba(cval)) > 2.0
                               else 'white'),
                        fontsize=fontsize)
            else:
                ax.text(x, (y * 1.01) + 0.04, '%.3f' % val, ha='center',
                        color=('black' if sum(mapper.to_rgba(cval)) > 2.5
                               else 'white'),
                        fontsize=fontsize)

    ax.set_xticks(xticks, minor=False)
    ax.set_yticks(yticks, minor=False)

    ax.invert_yaxis()

    if not xtitle is None:
        ax.set_xlabel(xtitle, fontsize=fontsize)
    if not ytitle is None:
        ax.set_ylabel(ytitle, fontsize=fontsize)
    if not xlabs is None:
        ax.set_xticklabels(xlabs, minor=False, fontsize=fontsize)
    if not ylabs is None:
        ax.set_yticklabels(ylabs, minor=False, fontsize=fontsize)

    return fig, ax

def scatter_set(labels, xs, ys):
    scatters = OrderedDict()
    for l, x, y in zip(labels, xs, ys):
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        if not isinstance(y, np.ndarray):
            y = np.array(y)
        scatters[l] = (x, y, np.array([True for _ in x]))
    return scatters

def deoverlap(boxplots, order):
    ovlpd = False
    boxes = [boxplots['boxes'][i] for i in order]
    medians = [boxplots['medians'][i] for i in order]
    if len(boxes) <= 1:
        return
    last = boxes[0]
    lastmed = medians[0]
    for i, (box, med) in enumerate(zip(boxes, medians)):
        if i == 0:
            continue
        cmed = med.get_xdata()
        lmed = lastmed.get_xdata()
        cur = box.get_xdata()
        las = last.get_xdata()
        ncur = np.copy(cur)
        nlas = np.copy(las)
        ovlp = las[1] - cur[0]
        if ovlp > 0.0:
            ovlpd = True
            ncur[0] += ovlp/2.0
            ncur[3] += ovlp/2.0
            ncur[4] += ovlp/2.0
            nlas[1] -= ovlp/2.0
            nlas[2] -= ovlp/2.0
            cmed[0] += ovlp/2.0
            lmed[1] -= ovlp/2.0
        med.set_xdata(cmed)
        medians[i-1].set_xdata(lmed)
        box.set_xdata(ncur)
        boxes[i-1].set_xdata(nlas)
        last = box
        lastmed = med
    return ovlpd

def gen_mh_multiscatter(scatters, xlabel=None, ylabel=None, one_line=False,
        colors=None, markers=None, connect=False, legend_title='conditions',
        showlegend=True, quartiles=False, trendlines=False, alpha=1.0,
        fig=None, ax=None):
    """ mask allows you to provide a mask for which values to supply to the
    quartiles and trendline analyses"""
    if colors is None:
        colors = mh_colors(len(scatters))
    if markers is None:
        markers = ['o' for _ in scatters]
    if quartiles:
        alpha=0.5
    for condname, (x, y, mask), c, m in zip(scatters.keys(), scatters.values(),
                                         colors, markers):
        fig, ax = gen_mh_scatter(x, y, color=c, xlabel=xlabel,
                                 ylabel=ylabel, one_line=one_line, fig=fig,
                                 ax=ax, marker=m, connect=connect,
                                 label=condname, trendline=trendlines,
                                 alpha=alpha, mask=mask)
        # only do this once as it adds insteads of sets...
        one_line = False
    if showlegend:
        box = ax.get_position()
        #print box.x0
        #print box.y0
        #print box.width
        #print box.height
        ax.set_position([box.x0 * 0.4, box.y0, box.width * 0.95, box.height])
        legend = ax.legend(loc='center left',
                           bbox_to_anchor=(1, 0.5), title=legend_title,
                           frameon=False, scatterpoints=5, fontsize=12,
                           markerscale=2.0)
        legend.get_title().set(fontsize='12', weight='bold',
                               horizontalalignment='left')
    if quartiles:
        new_scatters = OrderedDict()
        for condname, (x, y, mask) in scatters.iteritems():
            if len(x[mask]) >= 1:
                new_scatters[condname] = (x[mask], y[mask], [True for _ in
                    range(np.count_nonzero(mask))])
        if len(new_scatters) >= 1:
            fig, ax = gen_mh_scatter_boxplot(new_scatters, colors=colors,
                                             showlegend=False, fig=fig, ax=ax)
    return fig, ax

def gen_mh_scatter_boxplot(scatters, xlabel=None, ylabel=None, colors=None,
        legend_title='conditions', showlegend=True, fig=None, ax=None,
        fill=True, figsize=FIG_SIZE):
    if fig is None or ax is None:
        fig = plt.figure(figsize=figsize, dpi=200)
        ax = fig.add_subplot(111)
    xmin = scatters.values()[0][0][0]
    xmax = scatters.values()[0][0][0]
    vals = []
    positions = []
    for condname, (x, y, _) in zip(scatters.keys(), scatters.values()):
        if np.min(x) < xmin:
            xmin = np.min(x)
        if np.max(x) > xmax:
            xmax = np.max(x)
        vals.append(y)
        positions.append(np.mean(x))
    positions = np.array(positions)
    margin = (xmax - xmin) * 0.05

    opos = ax.get_position()
    xt = ax.get_xticks()
    tl = ax.get_xticks()
    xl = ax.get_xlim()
    box_artist = ax.boxplot(vals, positions=positions,
                            showfliers=False, showcaps=False,
                            showbox=True)
    overlap_found = deoverlap(box_artist, np.argsort(positions))
    for patch, color in zip(box_artist['boxes'], colors):
        #patch.set_color(color)
        patch.set_color('black')
        if fill:
            ax.add_patch(pat.Polygon(patch.get_xydata(), color=color))
        if overlap_found:
            patch.set_linewidth(patch.get_linewidth() * 0.5)

    ax.set_position(opos)
    ax.set_xticks(xt)
    ax.set_xticklabels(xt)
    xl = (xmin - margin, xmax + margin)
    ax.set_xlim(xl)
    for whisker in box_artist['whiskers']:
        whisker.set(ls='solid', color='#000000')
        if overlap_found:
            whisker.set_linewidth(whisker.get_linewidth() * 0.5)
    for median in box_artist['medians']:
        median.set(color='#000000')
        if overlap_found:
            median.set_linewidth(median.get_linewidth() * 0.5)
    return fig, ax

def polyfit_apply(coefs, x):
    final = 0
    for i, ce in enumerate(coefs):
        degree = len(coefs) - i - 1
        final += ce * (x ** degree)
    return final

def gen_mh_bandplot(xs, ys, color='blue',
                    fig=None, ax=None, figsize=FIG_SIZE):
    """Takes y quartiles (or whatever min, mid, max values you want to provide)
    at xs."""
    fig, ax = check_figax(fig, ax, figsize)
    #ymins, ymids, ymaxs = zip(*yquarts)
    #ax.plot(xs, ymids, '-', color=color)
    #ax.fill_between(xs, ymins, ymaxs, alpha=0.2, facecolor=color)
    ax.plot(xs, ys[1], '-', color=color)
    ax.fill_between(xs, ys[0], ys[2], alpha=0.2, facecolor=color)

    return fig, ax

def gen_mh_scatter(x, y, color='green', xlabel=None, ylabel=None,
                   one_line=False, fig=None, ax=None, marker='.',
                   connect=False, label=None, trendline=None, alpha=1.0,
                   mask=None, zero_lines=False, edgecolors='none', grid=None,
                   figsize=FIG_SIZE, linestyle='-'):
    """
    Notes:

    one_line:   False, <style, e.g. 'r--', 'r:'>
    zero_lines: T/F
    grid:       None, 'major, 'minor', 'both'
    edgecolors: 'none', 'green' etc
    zero_lines: T/F
    trendline:  False, 1, 2, ... (degree), 'lowess'

    """
    # some defaults if no style is provided, just True:
    if one_line is True:
        one_line = 'r--'

    if fig is None or ax is None:
        #fig, ax = plt.subplots()
        fig = plt.figure(figsize=figsize, dpi=200)
        ax = fig.add_subplot(111)

    if connect:
        ax.set_color_cycle([color])
        ax.plot(x, y, c=color, marker=marker, label=label, alpha=alpha)
    else:
        ax.scatter(x, y, c=color, marker=marker, edgecolors=edgecolors,
                   label=label, alpha=alpha)
    if not grid is None:
        ax.grid(b=True, which=grid)
    if xlabel:
        ax.set(xlabel=xlabel)
    if ylabel:
        ax.set(ylabel=ylabel)
    if one_line:
        lmin = max(min(x), min(y))
        lmax = min(max(x), max(y))
        ax.plot((lmin, lmax), (lmin, lmax), one_line)
    if zero_lines:
        lmin = max(min(x), min(y))
        lmax = min(max(x), max(y))
        ax.axvline(0)
        ax.axhline(0)

    if trendline:
        if not mask is None:
            x = x[mask]
            y = y[mask]
        if (len(set(x)) > 1 ) and (len(set(y)) > 1):
            if trendline == 'lowess':
                order = np.argsort(x)
                lx = x[order]
                ly = lowess(x[order], y[order])
            else:
                coefs = np.polyfit(x, y, trendline)
                lx = np.linspace(min(x), max(x), 100)
                ly = [polyfit_apply(coefs, x) for x in lx]
            ax.plot(lx, ly, '-', color=color)

    return fig, ax

def gen_density(data, bw=1, fig=None, ax=None, figsize=FIG_SIZE):
    if fig is None and ax is None:
        fig = plt.figure(figsize=figsize, dpi=200)
        ax = fig.add_subplot(111)

    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data), 200)
    density.covariance_factor = lambda: bw
    density._compute_covariance()
    ax.plot(xs, density(xs))
    return fig, ax

def stack_densities(conditions, labels, bw=1, xlabel=None, ylabel=None,
                    legend_title="condition", figsize=FIG_SIZE, fig=None,
                    ax=None, colors=None):
    fig, ax = check_figax(fig, ax, figsize)
    densities = []
    if colors is None:
        colors = mh_colors(len(labels))
    for data, color, label in zip(conditions, colors, labels):
        density = gaussian_kde(data)
        xs = np.linspace(min(data), max(data), 200)
        density.covariance_factor = lambda: bw
        density._compute_covariance()
        plot = ax.plot(xs, density(xs), color=color, label=label)
        densities.append(plot)

    #box = ax.get_position()
    #ax.set_position([box.x0 * 0.4, box.y0, box.width * 0.95, box.height])

    legend = ax.legend(loc='center left',
                       bbox_to_anchor=(1, 0.5), title=legend_title,
                       frameon=False, fontsize=12)
    legend.get_title().set(fontsize='12', weight='bold',
                           horizontalalignment='left')

    ax.margins(0.1)
    ax.get_yaxis().tick_left()
    if ylabel:
        ax.set(ylabel=ylabel)
    return fig, ax

def gen_mh_violin(data, labels=None, xlabel=None, ylabel=None, pos=None,
                  fig=None, ax=None, figsize=FIG_SIZE, title=None):
    """Data is a list of lists containing data to make into as many violin
    plots.
    """
    fig, ax = check_figax(fig, ax, figsize)
    if pos is None:
        pos = range(len(data))

    comps = ax.violinplot(data, pos, showextrema=True, showmedians=True)

    ymin, ymax = ax.get_ylim()
    margin = (ymax - ymin) * 0.05
    ax.set_ylim([ymin - margin, ymax + margin])
    ymin, ymax = ax.get_ylim()
    ymargin = (ymax - ymin) * 0.01
    xmin, xmax = ax.get_xlim()
    xmargin = (xmax - xmin) * 0.01
    for i, dset in enumerate(data):
        ax.text(i + xmargin,
                np.median(dset) + ymargin,
                str(np.median(dset)))

    if xlabel:
        ax.set(xlabel=xlabel)
    if ylabel:
        ax.set(ylabel=ylabel)
    if labels is None:
        ax.set_xticks([])
        ax.set_xticklabels([])
    else:
        ax.set_xticks(pos)
        ax.set_xticklabels(labels)
    if not title is None:
        ax.set_title(title)
    return fig, ax

def gen_histogram(data, ymin=None, ymax=None, height=None, width=None,
                  binwidth=1, bins=None, barwidth=None, fit_func=None,
                  step=False, ynorm=None, cumulative=False, color='steelblue',
                  ylabel=None, xlabel=None, flipx=False, fig=None, ax=None,
                  figsize=FIG_SIZE, linewidth=1, rethandles=False):
    """
    args:
        density is None or the bandwidth of the density chart
    """
    if not bins is None and binwidth != 0:
        log.warn("Using gen_histogram(bins=) trumps binwidth kwarg")
    fig, ax = check_figax(fig, ax, figsize)
    if not bins is None:
        pass
    elif binwidth == 1:
        bins = range(max(int(np.round(max(data))) + 1, 10))
    elif binwidth == 0:
        bins = 20
    else:
        print "sorry, bin widths not fully supported yet"
        bins = 20
    hist, bins = np.histogram(data, bins=bins)
    if barwidth is None:
        barwidth = (bins[1] - bins[0]) * 0.8
    bins = bins[:-1]
    if ynorm:
        hist = np.true_divide(hist, ynorm)
    if flipx:
        hist = hist[::-1]
        bins = bins[::-1]
        ax.invert_xaxis()
    if cumulative:
        tot = 0
        newhist = []
        #found = False
        for bar in hist:
            tot += bar
            #if tot >= 0.00006 and not found:
            #    print bins[len(newhist)]
            #    found = True
            newhist.append(tot)
        hist = newhist
    if step:
        handles = ax.plot(bins, hist, '-', color=color, linewidth=linewidth)
    else:
        handles = ax.bar(bins, hist, width=barwidth, color=color)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12)
    if fit_func:
        y = fit_func(hist)
        # The last bin value is the right edge of the last bin...
        ax.plot(bins, y, 'r--')
    if (not ymin is None) and (not ymax is None):
        log.info("Setting limits")
        ax.set_ylim([ymin, ymax])
    if rethandles:
        return fig, ax, handles
    return fig, ax

def mh_colors(number):
    if number <= 8:
        colors = [cm.get_cmap('Set1')(1. * i / 8.0)
                  for i in range(9)]
    else:
        colors = [cm.get_cmap('gist_rainbow')(1. * i / (number - 1))
                  for i in range(number)]
    return colors

def gen_mh_barplot(data, labels, ylabel=None, legend_title=None, fig=None,
                   ax=None, positions=None, groups=None, show_legend=True,
                   ret_handles=False, legend_coords=(1, 0.5),
                   legend_coord_system=None, colors=None,
                   xmin=None, xmax=None, ymin=None, ymax=None,
                   legend=None, fontsize=12, figsize=FIG_SIZE):
    """data is list of lists, labels is list of strings"""
    if not colors:
        colors = mh_colors(len(labels))
    fig, ax = check_figax(fig, ax, figsize)

    width = 1
    if positions is None:
        positions = range(len(data))
    if groups:
        width = 1.0/(max([len(group) for group in groups.values()]) + 1)
        positions = [0 for _ in data]
        curi = 0
        for conditions in groups.values():
            start = curi - (np.true_divide(len(conditions), 2) * width)
            for i, cond in enumerate(conditions):
                positions[cond] = start + i * width
            curi += 1
        if len(positions) > 10:
            fontsize=11
    box_artist = ax.bar(positions, data,
                        width=width, color=colors[:len(labels)],
                        linewidth=0)

    maxh = max([bar.get_height() for bar in box_artist.patches])
    shift = maxh * 0.01
    for bar, val in zip(box_artist.patches, data):
        plt.text(bar.get_x() + (0.5 * width), bar.get_height() + shift,
                 '{0:.3g}'.format(val),
                 ha='center',
                 fontsize=fontsize) # center didn't work

    if show_legend:
        lhandles = box_artist
        llabels = labels
        if legend:
            llabels = legend.keys()
            lhandles = [lhandles[i] for i in legend.values()]
        legend = ax.legend(lhandles, llabels, loc='center left',
                           bbox_to_anchor=legend_coords, title=legend_title,
                           frameon=False, fontsize=fontsize,
                           bbox_transform=legend_coord_system)
        legend.get_title().set(fontsize=fontsize, weight='bold',
                               horizontalalignment='left')

    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    ax.margins(0.1)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().tick_left()
    if ylabel:
        ax.set(ylabel=ylabel)
    if groups:
        ax.set_xticks(range(len(groups.keys())))
        ax.set_xticklabels(groups.keys())
        ax.get_xaxis().set_visible(True)

    if ret_handles:
        return fig, ax, box_artist
    else:
        return fig, ax

def gen_mh_boxplot(data, labels, ylabel=None, legend_title=None,
                   showfliers=False, fig=None, ax=None, groups=None,
                   show_legend=True, ret_patches=False, legend_coords=(1, 0.5),
                   legend_coord_system=None, colors=None,
                   xmin=None, xmax=None, ymin=None, ymax=None, legend=None,
                   fontsize=12, figsize=FIG_SIZE):
    """data is list of lists, labels is list of strings"""
    fig, ax = check_figax(fig, ax, figsize)
    width = 1
    positions = None
    if groups:
        width = 1.0/(max([len(group) for group in groups.values()]) + 1)
        positions = [0 for _ in data]
        curi = 0
        for conditions in groups.values():
            start = curi - (np.true_divide(len(conditions), 2) * width)
            for i, cond in enumerate(conditions):
                positions[cond] = start + i * width + width / 2.0
            curi += 1
        if len(positions) > 10:
            fontsize=11
    box_artist = ax.boxplot(data, labels=labels, patch_artist=True,
                            widths=width, positions=positions, showcaps=False,
                            showfliers=showfliers)

    if not colors:
        colors = mh_colors(len(labels))

    for patch, color in zip(box_artist['boxes'], colors):
        patch.set_facecolor(color)

    for whisker in box_artist['whiskers']:
        whisker.set(ls='solid', color='#000000')

    for line in box_artist['medians']:
        line.set(color='black')
        x1, y = line.get_xydata()[0]
        x2 = line.get_xydata()[1][0]
        ax.text(x1 + ((x2 - x1) / 2.0),
                y + 0.001,
                '%.3f' % y,
                ha='center',
                fontsize=fontsize) # center didn't work

    if show_legend:
        lhandles = box_artist['boxes']
        llabels = labels
        if legend:
            llabels = legend.keys()
            lhandles = [lhandles[i] for i in legend.values()]
        legend = ax.legend(lhandles, llabels, loc='center left',
                           bbox_to_anchor=legend_coords, title=legend_title,
                           frameon=False, fontsize=fontsize,
                           bbox_transform=legend_coord_system)
        legend.get_title().set(fontsize=fontsize, weight='bold',
                               horizontalalignment='left')

    ax.set_xlim(xmin=xmin, xmax=xmax)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    ax.get_yaxis().tick_left()
    if ylabel:
        ax.set(ylabel=ylabel)
    if groups:
        ax.get_xaxis().set_ticks(range(len(groups.keys())))
        ax.set_xticklabels(groups.keys())
        ax.get_xaxis().set_visible(True)
    else:
        ax.get_xaxis().set_visible(False)
        ax.get_xaxis().set_ticks(range(0, len(labels) + 2))
    ax.margins(0.1)

    if ret_patches:
        return fig, ax, box_artist
    else:
        return fig, ax

