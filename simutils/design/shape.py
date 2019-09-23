#!/usr/bin/env python

from shapely import geometry as geom
import numpy as np
from .space import Origin, Space, Parent, Child


class Polygons(geom.MultiPolygon):

    def __init__(self, polygons=None, context_type='polygons'):
        polygons = [polygons] if type(polygons) is geom.Polygon else polygons
        geom.MultiPolygon.__init__(self, polygons, context_type)

    def contains(self, p):
        return geom.MultiPolygon.contains(self, geom.Point(*p))


class Shape(Origin, Space):

    def __init__(self, rpos=[0, 0], parent=None):
        Origin.__init__(self, rpos, parent)

    @property
    def polygons(self):
        raise NotImplementedError()

    def plot(self, ax, **kwargs):
        for p in self.polygons:
            ax.plot(*p.exterior.xy, **kwargs)

    def plot_debug(self, ax, **kwargs):
        ax.plot(*self.pos, 'x', **kwargs)
        # self.plot_boundaries(ax)

    def unifie(self, *args):
        if not args:
            return Polygons()
        u = args[0]
        for p in args[1:]:
            u = u.union(p)
        return Polygons(u)

    def plot_debug_self(self, ax):
        pass


class World(Shape, Parent, dict):

    def __init__(self, margin=4, rpos=[0, 0]):
        Origin.__init__(self, rpos)
        Parent.__init__(self)
        dict.__init__(self)
        self.margin = margin

    def __iter__(self):
        return iter(self.values())

    def __setitem__(self, k, c):
        if not isinstance(c, Child):
            raise ValueError('World only accept Child objects')
        c.parent = self
        dict.__setitem__(self, k, c)

    @property
    def polygons(self):
        return self.unifie(*[c.polygons for c in self
                             if isinstance(c, Shape)])

    @property
    def extent(self):
        extent = self.compute_total_extent([s.extent for s in self])
        xmin, xmax, ymin, ymax = extent
        m = self.margin
        return [xmin - m, xmax + m, ymin - m, ymax + m]


class Waveguide(Shape, Parent):

    def __init__(self, linewidth, rpos=[0, 0]):
        Shape.__init__(self, rpos)
        Parent.__init__(self)
        self.linewidth = linewidth

    @property
    def points(self):
        raise NotImplementedError()

    @property
    def line(self):
        return geom.LineString(self.pos + self.points)

    @property
    def polygons(self):
        p = self.line.buffer(self.linewidth / 2,
                             cap_style=geom.CAP_STYLE.flat,
                             join_style=geom.JOIN_STYLE.bevel)
        return self.unifie(p, *[c.polygons for c in self.children
                                if isinstance(c, Shape)])

    @property
    def extent(self):
        b = np.asarray(self.polygons.bounds).reshape(2, 2)
        return b.T.reshape(-1)

    def transform(self, pos):
        return self.pos + self.points[-1] + pos

    def coordinates_at(self, percent):
        p = self.line.interpolate(percent, normalized=True)
        return np.asarray([i[0] for i in p.xy])

    def angle_at(self, percent):
        line = self.line
        segments = self.split_line(line)
        progress = np.cumsum([seg.length for seg in segments]) / line.length
        for p, seg in zip(progress, segments):
            if percent <= p:
                dx, dy = np.asarray(seg.coords[1]) - np.asarray(seg.coords[0])
                return np.arctan(dy / dx)

    def split_line(self, line):
        coords = line.coords
        return [geom.LineString(l) for l in zip(coords[0:], coords[1:])]


class LineWaveguide(Waveguide):

    def __init__(self, linewidth, dx, dy=0, rpos=[0, 0]):
        self.dx, self.dy = dx, dy
        Waveguide.__init__(self, linewidth, rpos)

    @property
    def points(self):
        return np.asarray([[0, 0], [self.dx, self.dy]])


class Arc:

    EXTENSION_LENGTH = 1e-2

    def compute_arc(self, r, a0, a1, res, extensions=None):
        a1 = a0 + 2 * np.pi if a1 == 'full' else a1
        a = np.linspace(a0, a1, int(self.res * abs(a1 - a0)) + 1)
        x, y = r * np.cos(a), r * np.sin(a)
        points = np.asarray([x, y]).T
        if extensions in {'start', 'both'}:
            e = self.compute_extension(points[0], a0, -self.EXTENSION_LENGTH)
            points = np.append([e], points, axis=0)
        if extensions in {'end', 'both'}:
            e = self.compute_extension(points[-1], a1, self.EXTENSION_LENGTH)
            points = np.append(points, [e], axis=0)
        return points - points[0]

    def compute_extension(self, p, a, l):
        a += 0.5 * np.pi
        return np.asarray(p) + np.asarray([np.cos(a), np.sin(a)]) * l


class ArcWaveguide(Waveguide, Arc):

    def __init__(self, linewidth, radius, angle0=0, angle1='full', res=64,
                 rpos=[0, 0]):
        self.radius, self.angle0, self.angle1 = radius, angle0, angle1
        self.res = res
        Waveguide.__init__(self, linewidth, rpos)

    @property
    def points(self):
        return self.compute_arc(self.radius, self.angle0, self.angle1,
                                self.res, extensions='both')


class SlineWaveguide(Waveguide, Arc):

    def __init__(self, linewidth, dx, dy, res=64, rpos=[0, 0]):
        self.dx, self.dy, self.res = dx, dy, res
        Waveguide.__init__(self, linewidth, rpos)

    def compute_radius(self, w, h):
        return (w**2 + h**2) / (4 * h)

    @property
    def points(self):
        dx, dy, res = self.dx, self.dy, self.res
        r = SlineWaveguide.compute_radius(self, dx, dy)
        a = np.arccos(dx / (2 * r)) - 0.5 * np.pi
        p1 = self.compute_arc(r, -0.5 * np.pi, a, res, extensions='start')
        p2 = self.compute_arc(r, a + np.pi, 0.5 * np.pi, res, extensions='end')
        p2 = p2 + p1[-1]
        return np.append(p1, p2[1:]).reshape(-1, 2)

'''
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


    def contains(self, pos):
        rpos = np.asarray(pos) - self.pos
        contains_self = self._shape.contains(Point(*rpos))
        return  contains_self or (self._child and self._child.contains(pos))


default_linewidth = 0.2
'''
