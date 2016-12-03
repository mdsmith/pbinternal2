###############################################################################
# Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################

# Author: Martin D. Smith


"""
Utils that support fringe DataSet features.
"""
import os
import tempfile
import logging
import json
import shutil
import datetime
import pysam
import numpy as np
from pbcore.util.Process import backticks

log = logging.getLogger(__name__)

def which(exe):
    if os.path.exists(exe) and os.access(exe, os.X_OK):
        return exe
    path = os.getenv('PATH')
    for this_path in path.split(os.path.pathsep):
        this_path = os.path.join(this_path, exe)
        if os.path.exists(this_path) and os.access(this_path, os.X_OK):
            return this_path
    return None

RS = 65536

def xy_to_hn(x, y):
    return x * RS + y

def hn_to_xy(hn):
    x = hn/RS
    y = hn - (x * RS)
    return x, y

def shift(cx, cy, d):
    # 0 is up, 1 is right, 2 is down, 3 is left
    if d == 0:
        cy += 1
    elif d == 1:
        cx += 1
    elif d == 2:
        cy -= 1
    elif d == 3:
        cx -= 1
    return cx, cy, d

def change_d(d):
    d += 1
    d %= 4
    return d

def move(cx, cy, x, y, d):
    if cx == x and cy == y:
        return cx - 1, y, 0
    if abs(x - cx) == abs(y - cy):
        d = change_d(d)
        # expand the search
        if d == 0:
            cx -= 1
        return shift(cx, cy, d)
    else:
        return shift(cx, cy, d)

def find_closest(x, y, pos, limit=81):
    found = False
    cx = x
    cy = y
    d = None
    fails = 0
    while not found:
        hn = xy_to_hn(cx, cy)
        if hn in pos:
            return hn
        else:
            fails += 1
            cx, cy, d = move(cx, cy, x, y, d)
        if fails >= limit:
            return None

def quadratic_expand(lol):
    samples = [[p] for p in lol[0]]
    for ps in lol[1:]:
        newsamples = []
        for p in ps:
            for s in samples:
                newsamples.append(s[:] + [p])
        samples = newsamples
    return samples

def prodround(values, target):
    """Round the floats in values (whose product is <target>) to integers in a
    way that minimizes the absolute change in values
    Args:
        values: a list of numbers
        target: the product of values (perhaps approximate)
    Returns:
        The values array, rounded to integers
    """
    opts = [[np.floor(v), round(v), np.ceil(v)] for v in values]
    combos = quadratic_expand(opts)
    best = combos[0]
    for combo in combos[1:]:
        p = np.prod(combo)
        err = abs(target - p)
        berr = abs(target - np.prod(best))
        rnd = np.sum([abs(v-c) for v, c in zip(values, combo)])
        brnd = np.sum([abs(v-c) for v, c in zip(values, best)])
        if (err < berr) or ((err == berr) and (rnd < brnd)):
            best = combo
    return best

def sampleUniformly(nsamples, dimbounds):
    """dimbounds is list of tuples of range, inclusive"""
    volume = 1
    for dmin, dmax in dimbounds:
        volume *= dmax - dmin
    volume_per_sample = np.true_divide(volume, nsamples)
    sample_side_length = np.power(volume_per_sample,
                                  np.true_divide(1.0, len(dimbounds)))
    per_axis = [max(1.0, np.true_divide((dmax - dmin), sample_side_length))
                for dmin, dmax in dimbounds]
    per_axis = prodround(per_axis, nsamples)
    # Shrink the stride to account for end margins
    strides = [np.true_divide(dmax - dmin, nsam + 1)
               for (dmin, dmax), nsam in zip(dimbounds, per_axis)]
    # introduce a margin
    points = [np.linspace(dmin + dstride,
                          dmax - dstride,
                          round(nsamp))
              for (dmin, dmax), nsamp, dstride in zip(dimbounds, per_axis,
                                                      strides)]
    points = [map(round, ps) for ps in points]
    points = [map(int, ps) for ps in points]
    samples = quadratic_expand(points)
    return samples


def sampleHolesUniformly(nsamples, samplefrom, faillimit=25, rowstart=64,
                         colstart=64, nrows=1024, ncols=1144):
    xys = sampleUniformly(nsamples, [(colstart, ncols), (rowstart, nrows)])
    hns = [find_closest(x, y, samplefrom, limit=faillimit) for x, y in xys]
    return [hn for hn in hns if not hn is None]

