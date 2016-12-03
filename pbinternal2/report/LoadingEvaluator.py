#!/usr/bin/env python

"""Evaluate the loading prediction method that produced an sts.h5 file"""

import argparse
import copy
import logging
import math
import os
from functools import partial
from random import shuffle
import numpy as np
import h5py
from matplotlib import pyplot as plt
from scipy.stats import poisson, chi2 as chi2dist
from pbcore.io.dataset.DataSetIO import _stackRecArrays
from pbinternal2.util.DataSetUtils import hn_to_xy
from pbinternal2.report.Graphs import (draw_chip, mh_colors,
                                       gen_mh_multiscatter, scatter_set,
                                       gen_histogram, gen_mh_scatter,
                                       gen_mh_boxplot, gen_mh_violin,
                                       stack_densities, savefig, lwrt)

from pbcommand.models.report import (Report, Attribute, Table, Column,
                                     PlotGroup, Plot)

log = logging.getLogger(__name__)

__version__ = "0.1.0"

class Load(object):
    EMPTY = 0
    SINGLE = 1
    MULTI = 2
    IND = 3

class Types(object):
    Empty = 0
    FullHqRead0 = 1
    FullHqRead1 = 2
    PartialHqRead0 = 3
    PartialHqRead1 = 4
    PartialHqRead2 = 5
    Indeterminate = 6

TypeMap = {Types.Empty: Load.EMPTY,
           Types.FullHqRead0: Load.SINGLE,
           Types.FullHqRead1: Load.SINGLE,
           Types.PartialHqRead0: Load.MULTI,
           Types.PartialHqRead1: Load.SINGLE,
           Types.PartialHqRead2: Load.SINGLE,
           Types.Indeterminate: Load.IND,
          }

class HQTypes(object):
    Empty = 0
    EarlyStarter = 1
    LateStarter = 2
    EarlyStarterPostSeq = 3
    LateStarterPostSeq = 4
    PreSeqLateStarter = 5
    PreSeqLateStarterPostSeq = 6
    Indeterminate = 7

HQTypeMap = {HQTypes.Empty: Load.EMPTY,
             HQTypes.EarlyStarter: Load.SINGLE,
             HQTypes.LateStarter: Load.SINGLE,
             HQTypes.EarlyStarterPostSeq: Load.MULTI,
             HQTypes.LateStarterPostSeq: Load.MULTI,
             HQTypes.PreSeqLateStarter: Load.MULTI,
             HQTypes.PreSeqLateStarterPostSeq: Load.MULTI,
             HQTypes.Indeterminate: Load.IND,
            }

class Constants(object):

    LOAD_PG = "loading_plots"
    PROD_PG = "productivity_plots"
    PROD_VS_PR_ID = "prod_vs_pulserate"
    REPORT_ID = "loading_vs_poisson"
    REPORT_TITLE = "Loading Prediction vs Poisson Expectation"

def write_csv(fname, csv, vlit=None):
    """vlit specifies one variable length iterable, any more would make even
    less sense in a csv"""
    header = ','.join(csv.dtype.names)
    if not vlit is None:
        header += ',{:}'.format(vlit.dtype.names[0])
    header += '\n'
    with open(fname, 'w') as ofh:
        ofh.write(header)
        for i, row in enumerate(csv):
            #csvrow = ','.join(map(str, list(row))) + ','
            csvrow = ','.join(map(str, list(row)))
            if not vlit is None:
                csvrow = csvrow + ','
                v = vlit[i][0]
                if v.dtype == bool:
                    v = v.astype(int)
                csvrow += ','.join(map(str, v))
            ofh.write(csvrow)
            ofh.write('\n')
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


def get_fake_perfect_loading(h5fh, lm=1.0, prox=20):
    hns = h5fh['/ZMW/HoleNumber'].value
    if lm == 'random':
        blocks = divide_holes(hns, prox)
        loads = np.zeros(len(hns))
        for block in blocks:
            lm = np.random.rand() * 5
            tloads = np.random.poisson(lm, len(block)).astype(np.float)
            loads[block] = tloads
    else:
        loads = np.random.poisson(lm, len(hns)).astype(np.float)
    loads[loads > 2] = 2.0
    return hns, np.array(loads)

get_variable_random_load = partial(get_fake_perfect_loading, lm='random')

def get_hns(h5fh):
    return h5fh['/ZMW/HoleNumber'].value

def get_prods(h5fh):
    return h5fh['/ZMWMetrics/Productivity'].value

def get_loads(h5fh):
    return h5fh['/ZMWMetrics/Loading'].value

def get_pulserates(h5fh):
    return h5fh['/ZMWMetrics/PulseRate'].value

def get_pulseratesVsT(h5fh):
    return h5fh['/ZMWMetrics/VsT/PulseRate'].value

def get_readtype_loading(h5fh):
    hns = h5fh['/ZMW/HoleNumber'].value
    labels = h5fh['/ZMWMetrics/ReadType'].value
    #print np.histogram(labels, bins=range(8))
    #print len(labels)
    for rtype, lval in TypeMap.iteritems():
        labels[labels == rtype] = lval
    return hns, labels

def get_hqreadtype_loading(h5fh):
    hns = h5fh['/ZMW/HoleNumber'].value
    labels = h5fh['/ZMWMetrics/HqReadType'].value
    #print np.histogram(labels, bins=range(8))
    #print len(labels)
    for rtype, lval in HQTypeMap.iteritems():
        labels[labels == rtype] = lval
    return hns, labels

def touches(xy, block):
    """Determine if xy, a coordinate, is contiguous with block, a list of
    coordinates. Note this is NOT the same as 'within the rectangle defined by
    the extent of coordinates within block'"""

def close(xy0, xy1, maxdist):
    if abs(xy0[0] - xy1[0]) > maxdist:
        return False
    elif abs(xy0[1] - xy1[1]) > maxdist:
        return False
    return True

def divide_holes(hns, prox):
    """Take a list of holenumbers, return a list of list of indices into the
    original list describing appropriate sets of holes to treat as
    observations.
    """
    # dumb greedy algorithm:
    xys = map(hn_to_xy, hns)
    components = {}
    last = xys[0]
    lmems = [0]
    for i, xy in enumerate(xys[1:]):
        i += 1
        found = False
        if close(xy, last, prox):
            lmems.append(i)
            found = True
        else:
            if last not in components:
                components[last] = lmems
            for comp, compmems in components.iteritems():
                if close(xy, comp, prox):
                    last = comp
                    lmems = compmems
                    lmems.append(i)
                    found = True
                if found:
                    break
            if not found:
                last = xy
                lmems = [i]
                found = True
    if last not in components:
        components[last] = lmems
    # filter components for size:
    dellist = []
    for comp, mems, in components.iteritems():
        if len(mems) < 20:
            dellist.append(comp)
    for comp in dellist:
        components.pop(comp)

    return map(np.array, components.values())

def make_dist(obs):
    return np.array([np.count_nonzero(obs == Load.EMPTY),
                     np.count_nonzero(obs == Load.SINGLE),
                     np.count_nonzero(obs == Load.MULTI),
                    ])

def make_dists(obss):
    return [make_dist(obs) for obs in obss]

def obs_chip(blocks, loadings):
    """Take blocks of indices that point to sections of loading predictions
    and produce a list of loadings for each block.
    """
    obs = []
    for block in blocks:
        vals = loadings[block].astype(np.float)
        obs.append(vals)
    return obs

def obs_hns(blocks, hns):
    """Take blocks of indices that point to sections of loading predictions
    and produce a list of loadings for each block.
    """
    if not isinstance(hns, np.ndarray):
        hns = np.array(hns)
    obs = []
    for block in blocks:
        vals = hns[block].astype(str)
        obs.append(vals)
    return obs

def remove_indeterminate_blocks(blocks, loadings, indfrac=0.2, remholes=True):
    """blocks are lists of lists"""
    newblocks = []
    for block in blocks:
        if np.count_nonzero(loadings[block] == Load.IND)/len(block) <= indfrac:
            if remholes:
                newblocks.append(block[loadings[block] != Load.IND])
            else:
                newblocks.append(block)
    return newblocks

def open_h5s(h5fns):
    return [h5py.File(stsh5, 'r') for stsh5 in h5fns]

def stsh5_to_observations(stsh5, accessor=get_readtype_loading,
                          outprefix=None, prox=20):
    """Go from an sts.h5 file, containing the output of a loading prediction
    method, to input for evalloading. This stsh5 file may be for a portion of
    the chip larger than what is appropriate for a single observation. It may
    therefore be necessary to break the stsh5 into physical chip portions, each
    with its own observed distribution of loading counts.
    """
    fh = h5py.File(stsh5, 'r')

    # get the per-zmw results:
    hns, loadings = accessor(fh)
    figs = []

    if outprefix:
        ofn = outprefix + '_chiproi.png'
        fig, ax = draw_chip(hns)
        savefig(fig, ofn)
        figs.append(('_'.join(['chip_roi_layout',
                               os.path.basename(outprefix)]),
                     ofn))

    # divide into regions of observation:
    blocks = divide_holes(hns, prox)

    blocks = remove_indeterminate_blocks(blocks, loadings)

    if outprefix:
        # indicate groups:
        ofn = outprefix + '_obsblocks.png'
        colors = mh_colors(len(blocks))
        shuffle(colors)
        fig, ax = None, None
        for block, col in zip(blocks, colors):
            # This will roll fig, ax through the list intentionally.
            fig, ax = draw_chip(hns[block], color=col, fig=fig, ax=ax)
        savefig(fig, ofn)
        figs.append(('_'.join(['observed_clusters',
                               os.path.basename(outprefix)]),
                     ofn))

        # hist of group sizes:
        ofn = outprefix + '_blocksizes.png'
        fig, ax = gen_histogram([len(bl) for bl in blocks])
        savefig(fig, ofn)
        figs.append(('_'.join(['blocksize_dist',
                               os.path.basename(outprefix)]),
                     ofn))

    # tally things up:
    obs = obs_chip(blocks, loadings)

    hnblocks = obs_hns(blocks, hns)

    return obs, hnblocks, figs

def poisson_curves(obss, names, outprefix):
    markers = 'ovs*pHDdooooooooooooooo'
    xs = []
    ys = []
    expnames = []
    colors = []
    expmarkers = []
    for obs, name, marker in zip(obss, names, markers):
        obs = make_dists(obs)
        x0 = []
        y0 = []
        x1 = []
        y1 = []
        for dist in obs:
            # NOTE: Swapped x,y relative to the matlab code so it looks the
            # same...
            # empty vs multi:
            y0.append(np.true_divide(dist[0], sum(dist)))
            x0.append(np.true_divide(dist[2], sum(dist)))

            # single vs non-empty
            y1.append(np.true_divide(dist[1], sum(dist)))
            x1.append(np.true_divide(dist[1] + dist[2], sum(dist)))
        expnames.extend(['{:}_emptyVSmulti'.format(name),
                         '{:}_singleVSnonempty'.format(name)])
        xs.extend([np.array(x0), np.array(x1)])
        ys.extend([np.array(y0), np.array(y1)])
        colors.extend(['green', 'blue'])
        expmarkers.extend([marker, marker])
    scatters = scatter_set(expnames, xs, ys)
    fig, ax = gen_mh_multiscatter(scatters, colors=colors, markers=expmarkers)
    evmxs = []
    evmys = []
    svnxs = []
    svnys = []
    for lm in np.linspace(0, 6, 1000):
        pmf = list(poisson(lm).pmf(range(2)))
        pmf.append(1.0 - sum(pmf))
        evmxs.append(pmf[2])
        evmys.append(pmf[0])
        svnxs.append(pmf[1] + pmf[2])
        svnys.append(pmf[1])
    ax.plot(evmxs, evmys, 'green')
    ax.plot(svnxs, svnys, 'blue')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])

    ofn = outprefix + '_poisson_curves.png'
    savefig(fig, ofn)
    plt.close(fig)
    return [('_'.join(['poisson_curves',
                       os.path.basename(outprefix)]),
             ofn)]

def expec(samp, lm):
    numer = (lm * (1 - np.exp(-1 * lm)))
    denom = (1 - np.exp(-1 * lm) * (lm + 1))
    newval = 0.0
    if denom != 0.0:
        newval = np.true_divide(numer, denom)
    samp[samp >= 2] = newval
    return samp

def maxim(samp):
    return np.mean(samp)

def estlambda(samp):
    steps = 0
    eqc = 0
    lm = maxim(samp)
    lastlm = lm
    while eqc <= 2 and steps <= 100:
        samp = expec(samp, lm)
        lm = maxim(samp)
        if abs(lm - lastlm) < (0.0001 * lastlm):
            eqc += 1
        else:
            lastlm = lm
            eqc = 0
        steps += 1
        yield lm, samp

def findlambda(samp):
    steps = 0
    eqc = 0
    lm = maxim(samp)
    lastlm = lm
    while eqc <= 2 and steps <= 100:
        samp = expec(samp, lm)
        lm = maxim(samp)
        if abs(lm - lastlm) < (0.0001 * lastlm):
            eqc += 1
        else:
            lastlm = lm
            eqc = 0
        steps += 1
    return lm

def graph_lambdas(lambdas, outprefix, names):
    if not isinstance(lambdas[0], (list, np.ndarray)):
        lambdas = [lambdas]
        names = [names]
    if len(lambdas) == 0:
        return
    elif len(lambdas[0]) == 0:
        return
    elif len(lambdas[0]) == 1:
        fig, ax = gen_mh_boxplot(lambdas, names, ylabel='Lambda Estimates')
    else:
        fig, ax = gen_mh_violin(lambdas, labels=names,
                                ylabel='Lambda Estimates',
                                title=os.path.split(outprefix)[1])
    ofn = outprefix + '_EstLambdas.png'
    savefig(fig, ofn)
    plt.close(fig)
    return [('_'.join(['lambda_ests',
                       os.path.basename(outprefix)]), ofn)]

def csv_to_dists(csv):
    obsds = []
    expds = []
    for row in csv:
        row = list(row)
        obsd = np.array(row[1:4])
        obsd = np.true_divide(obsd, sum(obsd))
        expd = np.array(row[4:7])
        expd = np.true_divide(expd, sum(expd))
        obsds.append(obsd)
        expds.append(expd)
    return obsds, expds

def graph_dists(obsds, expds, outprefix, chi2=None, df=None):
    # Graph vals:
    fig, ax = None, None
    for obsd, expd in zip(obsds, expds):
        fig, ax = gen_mh_scatter(range(3), obsd, fig=fig, ax=ax,
                                 connect=True, color='green')
        fig, ax = gen_mh_scatter(range(3), expd, fig=fig, ax=ax,
                                 connect=True, color='red')
    ax.set_title('{n}: Observed = Green, Expected = Red'.format(
        n=os.path.split(outprefix)[1]))
    ax.set(xlabel='Loading (Empty, Single, Multi)',
           ylabel='Fraction of holes')
    if not chi2 is None and not df is None:
        ax.text(
            .95, .95, 'chi**2/df = {c}\np = {p}'.format(
                c=np.true_divide(chi2, df), p=chi2_p(chi2, df)),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ofn = outprefix + '_ObsVsExpDists.png'
    savefig(fig, ofn)
    plt.close(fig)
    return [('_'.join(['obs_vs_exp_dists',
                       os.path.basename(outprefix)]),
             ofn)]

    """
    # Graph errors:
    fig, ax = None, None
    error = defaultdict(list)
    for row in csv:
        row = list(row)
        obsd = row[1:4]
        expd = row[4:7]
        for i, (o, e) in enumerate(zip(obsd, expd)):
            error[i].append(abs(o - e)/sum(obsd))
    fig, ax = gen_mh_boxplot(error.values(), error.keys())
    ax.set_title("{n}: Per Category Error Fraction Distribution".format(
                 n=os.path.split(outprefix)[1]))
    ofn = outprefix + '_ObsVsExpError.png'
    savefig(fig, ofn)
    plt.close(fig)

    # Bands
    xs = [0, 1, 2]
    categories = zip(*obsds)
    obsbounds = np.argsort(categories[0])
    #obsmins = obsds[obsbounds[0]]
    #obsmaxs = obsds[obsbounds[-1]]
    #obsmedians = obsds[obsbounds[len(obsbounds)/2]]
    obsmins = [min(cat) for cat in categories]
    obsmedians = [np.median(cat) for cat in categories]
    obsmaxs = [max(cat) for cat in categories]
    #fig, ax = gen_mh_bandplot(xs, (obsmins, obsmedians, obsmaxs),
    fig, ax = gen_mh_bandplot(xs, zip(obsmins, obsmedians, obsmaxs),
                              color='green')
    categories = zip(*expds)
    expbounds = np.argsort(categories[0])
    #expmins = expds[expbounds[0]]
    #expmaxs = expds[expbounds[-1]]
    #expmedians = expds[expbounds[len(expbounds)/2]]
    expmins = [min(cat) for cat in categories]
    expmedians = [np.median(cat) for cat in categories]
    expmaxs = [max(cat) for cat in categories]
    #fig, ax = gen_mh_bandplot(xs, (expmins, expmedians, expmaxs), fig=fig,
    fig, ax = gen_mh_bandplot(xs, zip(expmins, expmedians, expmaxs), fig=fig,
                              ax=ax, color='red')
    ofn = outprefix + '_ObsVsExpBands.png'
    savefig(fig, ofn)
    plt.close(fig)
    """

#from profilehooks import profile
#@profile
def evalloading(obs, hns=None, error='variance'):
    """Evaluate the predictive ability of a method given the observed
    distribution of loading counts generated by that method, in the form of a
    list of lists of counts:

    [[Empty, Singleload, Multiload],...]

    This observation is likely a subset of a chip, and therefore probably a
    subset of an sts.h5 file.
    """
    assert error in ['variance', 'expected']
    if hns is None:
        hns = [None for _ in obs]

    orig_obs = copy.deepcopy(obs)

    # for each observation estimate a lambda value:
    lambdas = []
    for samp in obs:
        lambdas.append(findlambda(samp))

    # Given the observed and expected distrbutions, characterize the
    # goodness of fit with a chi2 statistic:
    chi2 = 0
    rows = []
    header = [('estlmbd', np.float),
              ('obsempt', np.int),
              ('obssing', np.int),
              ('obsmult', np.int),
              ('expempt', np.float),
              ('expsing', np.float),
              ('expmult', np.float),
              ('chi2', np.float),
              ('hns', object),
             ]
    for od, lm, hns in zip(make_dists(orig_obs), lambdas, hns):
        ed = []
        this_chi2 = 0
        #print ''
        for i, oc in enumerate(od):
            if i < 2:
                pi = ((lm ** i) * np.exp(-1 * lm)) / math.factorial(i)
            else:
                p0 = ((lm ** 0) * np.exp(-1 * lm)) / math.factorial(0)
                p1 = ((lm ** 1) * np.exp(-1 * lm)) / math.factorial(1)
                pi = 1 - (p0 + p1)
            totsamples = sum(od)
            exp = totsamples * pi
            ed.append(exp)
            if error == 'variance':
                # Binomial variance denominator formula
                var = pi * (1 - pi) * totsamples
                #print var
                if var > 0.0:
                    this_chi2 += np.true_divide((oc - exp)**2, var)
            else:
                # Expected count denominator formula:
                # TODO: If this is a crazy small number it can blow up your
                # chi2
                if exp > 0.0:
                    this_chi2 += np.true_divide((oc - exp)**2, exp)
        chi2 += this_chi2

        #print lm
        #print od
        #print ed
        #print this_chi2
        row = [lm]
        row.extend(od)
        row.extend(ed)
        row.append(this_chi2)
        if hns is None:
            row.append('')
        else:
            row.append(':'.join(hns))
        rows.append(tuple(row))
    csv = np.array(rows, dtype=header).view(np.recarray)
    # degrees of fredom here should be the number of chi2 values added
    # together.
    df = len(obs)
    #print 'chi2, p: {c}, {p}'.format(c=chi2, p=pval)
    return chi2, df, csv

def chi2_p(chi2, df):
    chi2d = chi2dist(df)
    pval = 1.0 - chi2d.cdf(chi2)
    return pval

def evalloadings(obss, names, hnss=None, outprefix=None, error='variance'):
    cchi2 = 0
    cdf = 0
    csvs = []
    figs = []
    if outprefix:
        figs.extend(poisson_curves(obss, names, outprefix=outprefix))
    if hnss is None:
        hnss = [None for _ in obss]
    for obs, hns, name in zip(obss, hnss, names):
        chi2, df, csv = evalloading(obs, hns=hns, error=error)
        cchi2 += chi2
        cdf += df
        csvs.append(csv)
        if outprefix:
            csvfn = outprefix + '_{n}_perblockdist.csv'.format(n=name)
            write_csv(csvfn, csv)
            obsds, expds = csv_to_dists(csv)
            figs.extend(graph_dists(obsds, expds,
                                    outprefix + '_{:}'.format(name),
                                    chi2=chi2, df=df))
            figs.extend(graph_lambdas(csv.estlmbd,
                                      outprefix + '_{:}'.format(name),
                                      name))

    obsds, expds = csv_to_dists(_stackRecArrays(csvs))
    if outprefix:
        figs.extend(graph_dists(obsds, expds, outprefix, chi2=chi2, df=df))
    if len(csvs) > 1:
        lconds = [cond.estlmbd for cond in csvs]
        lconds.append(np.array([l for cond in csvs for l in cond.estlmbd]))
        figs.extend(graph_lambdas(lconds, outprefix, names[:] + ['all']))
    pval = chi2_p(cchi2, cdf)
    return np.true_divide(cchi2, cdf), pval, figs, csvs

def prodVSpulserate(stsh5fhs, names, outprefix):
    figs = []
    fig, axs = plt.subplots(nrows=len(stsh5fhs), ncols=2,
                            sharex='col', sharey='col', figsize=(8.5, 11),
                            dpi=200)
    if len(stsh5fhs) == 1:
        axs = [axs]
    for h5, name, rax in zip(stsh5fhs, names, axs):
        hns = get_hns(h5)
        prod = get_prods(h5)
        pulserates = get_pulseratesVsT(h5)
        p0prs = pulserates[prod == 0].flatten()
        p1prs = pulserates[prod == 1].flatten()
        p2prs = pulserates[prod == 2].flatten()
        ax = rax[0]
        ax.set_title('{n} Pulserate Histograms'.format(n=name))
        fig, ax, handles0 = gen_histogram(p0prs, bins=np.linspace(0, 3, 10),
                                          step=True, color='blue',
                                          linewidth=1, fig=fig, ax=ax,
                                          rethandles=True)
        fig, ax, handles1 = gen_histogram(p1prs, bins=np.linspace(0, 3, 10),
                                          step=True, fig=fig, ax=ax,
                                          color='green', linewidth=1,
                                          rethandles=True)
        fig, ax, handles2 = gen_histogram(p2prs, bins=np.linspace(0, 3, 10),
                                          step=True, fig=fig, ax=ax,
                                          color='red', linewidth=1,
                                          rethandles=True)
        ax.legend([handles0[0], handles1[0], handles2[0]],
                  ['p0', 'p1', 'p2'])
        ax.grid(True)
        ax.set_ylim(ymax=900000)

        ax = rax[1]
        ax.set_title('{n} Pulserates Densities'.format(n=name))
        dists = [dist for dist in [p0prs, p1prs, p2prs] if len(dist)]
        labels = [label for label, dist in zip(
            ['p0', 'p1', 'p2'], [p0prs, p1prs, p2prs]) if len(dist)]
        fig, ax = stack_densities(dists, labels,
                                  fig=fig, ax=ax, bw=.1,
                                  colors=['blue', 'green', 'red'])
        ax.set_ylim(ymin=0, ymax=6)
        ax.set_xlim(xmin=0, xmax=4)
        ax.grid(True)
    fig.tight_layout()

    ofn = outprefix + "_pulserateVSproductivity.png"
    savefig(fig, ofn)
    figs.append(('_'.join([Constants.PROD_VS_PR_ID,
                           os.path.basename(outprefix)]),
                 ofn))
    plt.close(fig)
    return figs


def analyze_loading(stsfns, names, outprefix, dist):
    figs = []
    obss = []
    hnss = []
    for stsfn, name in zip(stsfns, names):
        this_outprefix = outprefix + '_' + name
        obs, hns, ofigs = stsh5_to_observations(stsfn,
                                                accessor=get_readtype_loading,
                                                outprefix=this_outprefix,
                                                prox=dist)
        obss.append(obs)
        hnss.append(hns)
        figs.extend(ofigs)
    chi2, pval, lfigs, tables = evalloadings(obss, names, hnss,
                                             outprefix=outprefix)
    figs.extend(lfigs)

    return chi2, pval, figs, tables

def loading_table(stsfhs, names, outprefix):
    header = [('name', object),
              ('empty', np.int),
              ('single', np.int),
              ('multi', np.int),
              ('ind', np.int),
             ]
    rows = []
    for stsfh, name in zip(stsfhs, names):
        row = [name]
        loads = get_loads(stsfh)
        row.append(np.count_nonzero(loads == Load.EMPTY))
        row.append(np.count_nonzero(loads == Load.SINGLE))
        row.append(np.count_nonzero(loads == Load.MULTI))
        row.append(np.count_nonzero(loads == Load.IND))
        rows.append(tuple(row))
    csv = np.array(rows, dtype=header).view(np.recarray)
    ofn = outprefix + "_loading_distributions.csv"
    write_csv(ofn, csv)
    return ofn

def productivity_table(stsfhs, names, outprefix):
    header = [('name', object),
              ('empty', np.int),
              ('prod', np.int),
              ('other', np.int),
             ]
    rows = []
    for stsfh, name in zip(stsfhs, names):
        row = [name]
        prods = get_prods(stsfh)
        row.append(np.count_nonzero(prods == 0))
        row.append(np.count_nonzero(prods == 1))
        row.append(np.count_nonzero(prods == 2))
        rows.append(tuple(row))
    csv = np.array(rows, dtype=header).view(np.recarray)
    ofn = outprefix + "_productivity_distributions.csv"
    write_csv(ofn, csv)
    return ofn

# This will need to be re-written for the conditional comparison internal
# report
def loading_report(stsfns, names, outdir, dist):
    # for lack of a singular name:
    name = os.path.join(outdir, 'loading')
    outprefix = os.path.join(outdir, name)

    chi2, pval, lfigs, tables = analyze_loading(stsfns, names, outprefix, dist)

    pfigs = []
    stsfhs = open_h5s(stsfns)
    pfigs.extend(prodVSpulserate(stsfhs, names,
                                 outprefix))
    ldist_tables = loading_table(stsfhs, names, outprefix)
    pdist_tables = productivity_table(stsfhs, names, outprefix)

    # Attributes
    attrs = [("loading_vs_poisson_version", __version__,
              "Loading Evaluation Version"),
             ("dataset_chi2_per_dof", chi2, "Chi2/df"),
             ("dataset_pval", pval, "p-value")]

    attributes = [Attribute(i, v, name=n) for i, v, n in attrs]

    # Plots
    lpg = PlotGroup(Constants.LOAD_PG)
    for plot_id, png in lfigs:
        plot = Plot(plot_id, os.path.relpath(png))
        lpg.add_plot(plot)

    ppg = PlotGroup(Constants.PROD_PG)
    for plot_id, png in pfigs:
        plot = Plot(plot_id, os.path.relpath(png))
        ppg.add_plot(plot)

    plotgroups = [lpg, ppg]

    # Tables
    pbctables = []
    for table in tables:
        pbccolumns = [Column(name, header=name, values=table[name])
                      for name in table.dtype.names]
        pbctables.append(Table('loading_table',
                               title='Loading Prediction Quality Metrics',
                               columns=pbccolumns))

    report = Report(Constants.REPORT_ID,
                    Constants.REPORT_TITLE,
                    attributes=attributes,
                    plotgroups=plotgroups,
                    tables=pbctables)
    return report

def loading_vs_poisson(stsh5_fn, report_fn, nproc, dist):
    outdir = os.path.dirname(report_fn)
    name = swap_extension(os.path.basename(stsh5_fn), '')
    report = loading_report([stsh5_fn], [name], outdir, dist)
    report.write_json(report_fn)
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pad = parser.add_argument
    pad('stsh5', type=str, nargs='+',
        help=('One or more sts.h5 files to analyze'))
    pad('-n', '--name', type=str, default=[], action='append',
        help=('Name for each stsh5'))
    pad('-o', '--outprefix', type=str,
        help=('A prefix to saved analysis information'))
    pad('-d', '--dataset', type=str, default='ReadType',
        help=('The name of the sts.h5 dataset to pull loading predictions '
              'from'))
    pad('--sim', default=False, action='store_true',
        help=('Simulate the loading observations'))
    pad('--dist', default=20, type=int,
        help=('Size of observation blocks'))
    pad('--prodonly', default=False, action='store_true')
    pad('--denom', default='variance', choices=['variance', 'expected'])
    args = parser.parse_args()
    if not args.prodonly:
        accessor = get_readtype_loading
        #accessor = get_hqreadtype_loading
        if args.sim:
            accessor = partial(get_variable_random_load, prox=args.dist)
            #accessor = get_fake_perfect_loading
        if len(args.name) == 0:
            args.name = [os.path.split(args.outprefix)[1] for _ in args.stsh5]
        elif len(args.name) == 1:
            args.name = [args.name[0] for _ in args.stsh5]
        else:
            assert len(args.name) == len(args.stsh5), (
                'Please provide 0, 1 or '
                'a number of names equal to the number of sts.h5 files')
        obss = [stsh5_to_observations(
                    stsh5,
                    outprefix=('{p}_{n}'.format(p=args.outprefix, n=name)),
                    accessor=accessor, prox=args.dist)
                for stsh5, name in zip(args.stsh5, args.name)]
        obss, hnss, ofigns = zip(*obss)
        chi2, pval, figns, csvs = evalloadings(obss, args.name, hnss=hnss,
                                               outprefix=args.outprefix,
                                               error=args.denom)
    stsfhs = open_h5s(args.stsh5)
    prodVSpulserate(stsfhs, args.name, args.outprefix)
    loading_table(stsfhs, args.name, args.outprefix)
    productivity_table(stsfhs, args.name, args.outprefix)
