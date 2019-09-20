#!/usr/bin/env python

import numpy as np
from weakref import ref
from itertools import chain

class Parent:

    def __init__(self):
        self.children = []

    def __iter__(self):
        return iter(self.children)

    def get_all(self):
        return [c for c in chain(c for c in self if isinstance(c, Parent))]

    def attach(self, c):
        if not isinstance(c, Child):
            raise ValueError('Parent only accept Child objects')
        c.parent = self
        self.children.append(c)


class Child:

    def __init__(self, parent=None):
        self._parent = parent

    @property
    def parent(self):
        return self._parent and self._parent()

    @parent.setter
    def parent(self, p):
        self._parent = ref(p) if p else None


class Origin(Child):

    def __init__(self, rpos=[0, 0], parent=None):
        self.rpos = np.asarray(rpos)
        Child.__init__(self, parent)

    def transform(self, rpos):
        return self.pos + rpos

    @property
    def pos(self):
        return self.parent.transform(self.rpos) if self.parent else self.rpos


class Space:

    @property
    def extent(self):
        raise NotImplementedError()

    @property
    def size(self):
        xmin, xmax, ymin, ymax = self.extent
        return xmax - xmin, ymax - ymin

    @property
    def center(self):
        xmin, xmax, ymin, ymax = self.extent
        return (xmax + xmin) / 2, (ymax + ymin) / 2

    def compute_total_extent(self, extents):
        xmin, xmax, ymin, ymax = np.asarray(extents).T
        xmin, ymin = np.min(xmin), np.min(ymin)
        xmax, ymax = np.max(xmax), np.max(ymax)
        return [xmin, xmax, ymin, ymax]


def attach(*args):
    for s0, s1 in zip(args, args[1:]):
        s0.attach(s1)
    return args
