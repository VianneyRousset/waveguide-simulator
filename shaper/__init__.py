#!/usr/bin/env python

# requierements:
# - meep
# - shapely
# - numpy
# - scipy
# - matplotlib

from shapely.geometry import Point, MultiLineString, CAP_STYLE, JOIN_STYLE
import numpy as np
from weakref import ref

default_linewidth = 0.2


def attach(*args):
    for s0, s1 in zip(args, args[1:]):
        s0.attach(s1)
    return args


class Origin:

    def transform(self, pos):
        return pos


class Space(Origin):

    def __init__(self, rpos=[0, 0]):
        Origin.__init__(self)
        self._parent = None
        self.rpos = np.asarray(rpos)

    @property
    def parent(self):
        return self._parent and self._parent()

    @parent.setter
    def parent(self, p):
        self._parent = ref(p) if p else None

    @property
    def pos(self):
        return self.parent.transform(self.rpos) if self.parent else self.rpos

    @property
    def boundaries(self):
        return [0, 0, 0, 0]

    @property
    def size(self):
        xmin, xmax, ymin, ymax = self.boundaries
        return xmax - xmin, ymax - ymin

    def transform(self, pos):
        return self.pos + np.asarray(pos)

    def contains(self, p):
        return False

    def rasterize(self, boundaries='auto', res=[32, 32]):
        resx, resy = res
        b = boundaries if boundaries is not 'auto' else self.boundaries
        xmin, xmax, ymin, ymax = b
        w, h = xmax - xmin, ymax - ymin
        X = np.linspace(xmin, xmax, w * resx)
        Y = np.linspace(ymin, ymax, h * resy)
        return self.rasterize_chunk(X, Y)

    def rasterize_chunk(self, X, Y):
        v = [[self.contains([x, y]) for x in X] for y in Y]
        return np.array(v, dtype=bool)

    def plot_raster(self, axis, res=[64, 64]):
        img = self.rasterize()
        axis.imshow(img, extent=self.boundaries, origin='lower')

    def draw(self):
        raise NotImplementedError()

    def draw_self(self):
        raise NotImplementedError()

    def plot(self, axis):
        self.plot_self(axis)

    def plot_self(self, axis):
        axis.plot(*self.pos, 'x')
        self.plot_boundaries(axis)

    def plot_boundaries(self, axis):
        x0, x1, y0, y1 = self.boundaries
        p = [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]
        axis.plot(*np.asarray(p).T, '--')

    def compute_total_boundaries(self, boundaries):
        xmin, xmax, ymin, ymax = np.asarray(boundaries).T
        xmin, ymin = np.min(xmin), np.min(ymin)
        xmax, ymax = np.max(xmax), np.max(ymax)
        return [xmin, xmax, ymin, ymax]


class World(Space, dict):

    def __init__(self, rpos=[0, 0]):
        Space.__init__(self, rpos)
        dict.__init__(self)

    def __setitem__(self, key, value):
        if not isinstance(value, Space):
            raise ValueError('Wold only accept spaces')
        if key in self:
            self[key].parent = None
        value.parent = self
        dict.__setitem__(self, key, value)

    @property
    def boundaries(self):
        return self.compute_total_boundaries([s.boundaries
                                              for s in self.values()])

    def contains(self, pos):
        for s in self.values:
            if s.contains(pos):
                return True
        return False

    def plot(self, axis):
        self.plot_self(axis)
        for s in self.values():
            s.plot(axis)


class Rectangle(Space):

    def __init__(self, width, height, rpos=[0, 0]):
        Space.__init__(self, rpos)
        self.width, self.height = width, height

    def contains(self, p):
        x, y = np.asarray(p) - self.pos
        xmin, xmax, ymin, ymax = self.boundaries
        return xmin <= x <= xmax and ymin <= y <= ymax

    @property
    def boundaries(self):
        x, y = self.pos
        w, h = self.width, self.height
        return [x, x + w, y, y + h]

    def plot_self(self, axis):
        w, h = self.width, self.height
        x0, y0 = self.pos
        x1, y1 = x0 + w, y0 + h
        axis.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0])
        Space.plot_self(self, axis)


class Waveguide(Space):

    def __init__(self, points, linewidth=default_linewidth, rpos=[0, 0]):
        Space.__init__(self, rpos)
        self._points = np.asarray(points)
        self._linewidth = linewidth
        self._shape = self._create_shape()
        self._child = None

    def contains(self, pos):
        rpos = np.asarray(pos) - self.pos
        contains_self = self._shape.contains(Point(*rpos))
        return  contains_self or (self._child and self._child.contains(pos))

    @property
    def linewidth(self):
        return self._linewidth

    @linewidth.setter
    def linewidth(self, w):
        self._linewidth = w
        self._shape = self._create_shape()

    @property
    def child(self):
        return self._child

    @child.setter
    def child(self, c):
        if not isinstance(c, Space) and c is not None:
            raise ValueError('Child can only be a space or None')
        if self.child:
            self.child.parent = None
        if c:
            c.parent = self
        self._child = c

    @property
    def points(self):
        return self._points

    @points.setter
    def points(self, p):
        self._points = p
        self._shape = self._create_shape()

    @property
    def boundaries(self):
        b = np.asarray(self._shape.bounds).reshape(2, 2)
        b = (self.pos + b).T.reshape(-1)
        if self.child:
            return self.compute_total_boundaries([b, self.child.boundaries])
        return b

    def attach(self, c):
        self.child = c
        return self

    def detach(self):
        self.child = None

    def transform(self, pos):
        return self.pos + self.points[-1] + pos

    def plot(self, axis):
        self.plot_self(axis)
        if self.child:
            self.child.plot(axis)

    def plot_self(self, axis):
        p = self.points + self.pos
        axis.plot(*p.T, linewidth=2)
        Space.plot_self(self, axis)

    def _create_shape(self):
        line = MultiLineString([self.points])
        return line.buffer(self.linewidth / 2, cap_style=CAP_STYLE.flat,
                           join_style=JOIN_STYLE.bevel)


class LineWaveguide(Waveguide):

    def __init__(self, width, height=0, linewidth=default_linewidth,
                 rpos=[0, 0]):
        self._width, self._height = width, height
        points = self.compute_points(width, height)
        Waveguide.__init__(self, points, linewidth, rpos)

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, w):
        self._width = w
        self.points = self.compute_points(self.width, self.height)

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, w):
        self._height = w
        self.points = self.compute_points(self.width, self.height)

    def compute_points(self, w, h):
        return np.asarray([[0, 0], [w, h]])


class ArcWaveguide(Waveguide):

    def __init__(self, radius, a0=0, a1='full', linewidth=default_linewidth,
                 res=64, rpos=[0, 0]):
        points = self.compute_points(radius, a0, a1, res)
        Waveguide.__init__(self, points, linewidth, rpos)

    def compute_points(self, r, a0, a1, res):
        a1 = a0 + 2 * np.pi if a1 == 'full' else a1
        a0, a1 = a0 - np.pi / 2, a1 - np.pi / 2
        a = np.linspace(a0, a1, int(res * abs(a1 - a0)) + 1)
        x, y = r * np.cos(a), r * np.sin(a)
        points = np.asarray([x, y]).T
        return points - points[0]


class SlineWaveguide(Waveguide):

    def __init__(self, width, height, linewidth=default_linewidth, res=64,
                 rpos=[0, 0]):
        self._width, self._height, self._res = width, height, res
        points = self.compute_points(width, height, res)
        Waveguide.__init__(self, points, linewidth, rpos)

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, w):
        self._width = w
        self.points = self.compute_points(self.width, self.height, self.res)

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, w):
        self._height = w
        self.points = self.compute_points(self.width, self.height, self.res)

    @property
    def radius(self):
        return self.compute_radius(self.width, self.height)

    def compute_radius(self, w, h):
        return (w**2 + h**2) / (4 * h)

    def compute_points(self, w, h, res):
        r = SlineWaveguide.compute_radius(self, w, h)
        a = np.pi / 2 - np.arccos(w / (2 * r))
        p1 = ArcWaveguide.compute_points(self, r, 0, a, res)
        p2 = ArcWaveguide.compute_points(self, r, a + np.pi, np.pi, res)
        p2 = p2 + p1[-1]
        return np.append(p1, p2[1:]).reshape(-1, 2)


class MZIWaveguide(Waveguide):

    def __init__(self, width, height, n, linewidth=default_linewidth, res=64,
                 rpos=[0, 0]):
        self._width, self._height, self._n, self._res = width, height, n, res
        points = self.compute_points(width, height, n, res)
        Waveguide.__init__(self, points, linewidth, rpos)

    @property
    def width(self):
        return self._width

    @width.setter
    def width(self, w):
        self._width = w
        self.points = self.compute_points(self.width, self.height, self.res)

    @property
    def height(self):
        return self._height

    @height.setter
    def height(self, w):
        self._height = w
        self.points = self.compute_points(self.width, self.height, self.res)

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, n):
        self._n = n
        self.points = self.compute_points(self.width, self.height, self.n,
                                          self.res)

    @property
    def radius(self):
        return SlineWaveguide.compute_radius(self, self.width, self.height)

    def compute_points(self, w, h, n, res):
        heights = [-h if i % 2 else h for i in range(2 * n)]
        points = np.asarray([[0, 0]])
        for h in heights:
            pnew = SlineWaveguide.compute_points(self, w, h, res)[1:]
            points = np.append(points, pnew + points[-1], axis=0)
        return points
