#!/usr/bin/env python
"""Simple storage classes for performing one-dimensional set operations on
integer ranges."""

import bisect


class Range(object):
    """Simple storage class for performing one-dimensional set operations on
    integer ranges.

    Ranges are stored as [start,end) with start<end (i.e. end is exclusive)
    """

    def __init__(self, start=0, end=0):
        start = int(start)
        end = int(end)
        if end < start:
            self.start = end
            self.end = start
        else:
            self.start = start
            self.end = end

    def get_length(self):
        return len(self)

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def copy(self):
        return Range(self.start, self.end)

    def __getitem__(self, i):
        if i == 0:
            return self.start
        if i == 1:
            return self.end
        raise Exception("Can't access key=%s in Range" % str(i))

    def __iter__(self):
        return xrange(self.start, self.end).__iter__()

    def __len__(self):
        return self.end - self.start

    def __contains__(self, x):
        return x >= self.start and x < self.end

    def __str__(self):
        return '[%d,%d)' % (self.start, self.end)

    def __eq__(self, v):
        return self.start == v.start and self.end == v.end

    def add_delta(self, delta):
        self.start += delta
        self.end += delta

    def intersects(self, range2):
        return self.contains(range2.start) or range2.contains(self.start)

    def intersect(self, range2):
        "Returns Range corresponding to intersection."
        if not self.intersects(range2):
            return Range(0, 0)
        return Range(max(self.start, range2.start), min(self.end, range2.end))

    def contains(self, x):
        return x in self

    def contains_range(self, other_range):
        return (self.contains(other_range.start) and
                self.contains(other_range.end - 1))

    def union(self, range2):
        """Replaces self with a range which encompasses the current broadest
        span of self and range2's limits"""
        self.start = min(range2.start, self.start)
        self.end = max(range2.end, self.end)


class Ranges(object):
    """Represents an ordered and non-overlapping collection of ranges.
    Not an efficient implementation for large N.
    TODO: more efficient would be i.e. an interval tree"""

    def __init__(self, key=lambda x: x.start):
        self._ranges = []
        self._key = key

    def add_range(self, pbrange):
        """Add a new range or extend an existing range"""
        bisect.insort_right(self._ranges, (self._key(pbrange), pbrange.copy()))
        self.__normalize()

    def remove_range(self, pbrange):
        """Remove a range or shrink an existing range"""
        del_list = []
        new_list = []
        for i, r in enumerate(self._ranges):
            if r[1].intersects(pbrange):
                if pbrange.contains_range(r[1]):
                    del_list.append(i)
                elif r[1].contains_range(pbrange):
                    del_list.append(i)
                    r1 = Range(r[1].start, pbrange.start)
                    if len(r1) > 0:
                        new_list.append(r1)
                    r2 = Range(pbrange.end, r[1].end)
                    if len(r2) > 0:
                        new_list.append(r2)
                else:
                    if r[1].start < pbrange.start:
                        r[1].end = pbrange.start
                    else:
                        r[1].start = pbrange.end
        n_del = 0
        for i in del_list:
            del self._ranges[i - n_del]
            n_del += 1
        for r in new_list:
            self.add_range(r)

    def __normalize(self):
        """Ensure that the class invariant is maintained.
           Namely that the list of ranges is ordered and non-overlapping.
        """
        if len(self._ranges) < 2:
            return
        del_list = []
        j = 0
        for i in xrange(0, len(self._ranges)):
            if i == j:
                continue
            r1 = self._ranges[j][1]
            r2 = self._ranges[i][1]
            if r1.intersects(r2):
                r1.union(r2)
                del_list.append(i)
            else:
                j = i
        n_del = 0
        for i in del_list:
            del self._ranges[i - n_del]
            n_del += 1

    def __len__(self):
        if len(self._ranges) == 0:
            return 0
        return self._ranges[-1][1].end - self._ranges[0][1].start

    def __iter__(self):
        """Iterates over Range"""
        for pbrange in self._ranges:
            yield pbrange[1]

    def __contains__(self, point):
        """Returns True if integer x is contained in these ranges."""
        for pbrange in self._ranges:
            if point in pbrange:
                return True
        return False

    def span(self):
        """Returns a Range object representing the span of these Ranges"""
        if len(self._ranges) == 0:
            return Range(0, 0)
        return Range(self._ranges[0][1].start, self._ranges[-1][1].end)

    def intersects(self, pbrange2):
        """Returns True if range2 intersects any Range in self."""
        for pbrange in self:
            if pbrange.intersects(pbrange2):
                return True
        return False

    def get_intersecting_range(self, pbrange2):
        """Returns the Range in self that overlaps range2 or None if
        no range overlaps."""
        for pbrange in self:
            if pbrange.intersects(pbrange2):
                return pbrange
        return None

    def get_intersecting_ranges(self, pbrange2):
        """Returns the Ranges in self that overlap range2"""
        for pbrange in self:
            if pbrange.intersects(pbrange2):
                yield pbrange

    def __str__(self):
        """For debugging"""
        return ('Rs {' + ' '.join(
            [str(pbrange[1]) for pbrange in self._ranges]) + '}')

    def gaps(self):
        """Iterates over 'gaps', which are defined by pairs of flanking \
        ranges."""
        for i in range(1, len(self._ranges)):
            yield (self._ranges[i - 1][1], self._ranges[i][1])

    def merge(self, ranges):
        """Merges the ranges in ranges (class Ranges) with this object."""
        for pbrange in ranges:
            bisect.insort_right(self._ranges, (self._key(pbrange),
                                               pbrange.copy()))
        self.__normalize()


class OverlappingRanges(object):
    """Represents an ordered and potentially overlapping collection of \
    ranges."""

    def __init__(self, key=lambda x: x.start):
        self._ranges = []
        self._key = key

    def add_range(self, pbrange):
        bisect.insort_right(self._ranges, (self._key(pbrange), pbrange.copy()))

    def __len__(self):
        if len(self._ranges) == 0:
            return 0
        return self._ranges[-1][1].end - self._ranges[0][1].start

    def __iter__(self):
        for r in self._ranges:
            yield r[1]

    def span(self):
        if len(self._ranges) == 0:
            return Range(0, 0)
        return Range(self._ranges[0][1].start, self._ranges[-1][1].end)

    def __str__(self):
        return 'Rs {' + ' '.join([str(r[1]) for r in self._ranges]) + '}'

    def overlapping_ranges(self, query):
        # TODO improve speed?

        start_idx = bisect.bisect(self._ranges, query.get_start())
        for (end, target) in self._ranges[start_idx:]:
            if target.intersects(query):
                yield target


# unit tests
if __name__ == '__main__':

    #--- Range unit tests
    r1 = Range(5, 10)
    r2 = Range(10, 15)
    r3 = Range(7, 15)
    r4 = Range(5, 10)

    print 'r1 = %s' % str(r1)
    print 'r2 = %s' % str(r2)
    print 'r3 = %s' % str(r3)
    print 'r4 = %s' % str(r4)

    print 'r1.contains(5) = %d' % r1.contains(5)
    print 'r1.contains(10) = %d' % r1.contains(10)

    print 'r1.intersects(r2) = %d' % r1.intersects(r2)
    print 'r1.intersects(r3) = %d' % r1.intersects(r3)

    print 'r1.intersect(r2) = %s' % r1.intersect(r2)
    print 'r1.intersect(r3) = %s' % r1.intersect(r3)

    print 'r1==r4 = %s' % str(r1 == r4)

    #--- Ranges unit tests
    r = Ranges()
    r.add_range(Range(1, 3))
    r.add_range(Range(9, 12))
    print '------'
    print 'r = [1,3)+[9,12) = %s' % str(r)
    r.add_range(Range(2, 5))
    print 'r + [2,5) = %s' % str(r)
    r.add_range(Range(1, 12))
    print 'r + [1,12) = %s' % str(r)
    r.add_range(Range(14, 15))
    print 'r + [14,15) = %s' % str(r)
    r.add_range(Range(20, 25))
    print 'r + [20,25) = %s' % str(r)
    r.add_range(Range(11, 22))
    print 'r + [11,22) = %s' % str(r)

    r5 = Ranges()
    r5.add_range(Range(1, 25))
    r5.add_range(Range(27, 29))
    r5.add_range(Range(35, 40))
    r6 = Ranges()
    r6.add_range(Range(2, 5))
    r6.add_range(Range(20, 30))
    r6.add_range(Range(42, 45))
    print '------'
    print 'r5 = %s' % str(r5)
    print 'r6 = %s' % str(r6)
    r5.merge(r6)
    print 'r5.merge(r6) = %s' % str(r5)

    # now r is 1,25
    print '------'
    print 'r = %s' % str(r)
    print 'r1 = %s' % str(r1)
    r.remove_range(r1)
    print 'r - r1 = %s' % str(r)
    # now r is [1,5) [10,25)
    r5 = Range(3, 12)
    print 'r5 = %s' % str(r5)
    r.remove_range(r5)
    print 'r - r5 = %s' % str(r)
    r.remove_range(Range(-1, 3))
    print 'r - [-1,3) = %s' % str(r)
    r.remove_range(Range(1, 3))
    print 'r - [1,3) = %s' % str(r)
    r.remove_range(Range(1, 25))
    print 'r - [1,25) = %s' % str(r)

    #--- OverlappingRanges unit tests
    r = OverlappingRanges()
    r.add_range(Range(1, 3))
    r.add_range(Range(0, 15))
    r.add_range(Range(9, 12))
    print ('\nOverlappingRanges tests\n\nr = [1,3)+[0,15)+[9,12) = %s\n' %
           str(r))

    query_ranges = [Range(0, 1), Range(2, 4), Range(13, 15), Range(15, 17)]
    for q_range in query_ranges:
        print "query %s is overlapped by:" % str(q_range)
        overlapping_length_sum = 0
        for o_range in r.overlapping_ranges(q_range):
            print str(o_range)
