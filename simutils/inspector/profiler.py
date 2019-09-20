import numpy as np
from math import ceil
from shapely import geometry as geom
from .plot import Plot


class Profile(Plot):

    def __init__(self, sim_name, name, vx, vy):
        Plot.__init__(self, sim_name, name)
        self.vx, self.vy = vx, vy

    def plot(self, ax, *args, **kwargs):
        ax.plot(self.vx, self.vy, *args, **kwargs)

    @property
    def abs(self):
        return Profile(self.sim_name, f'{self.name}_abs', self.vx,
                       np.abs(self.vy))

    @property
    def angle(self):
        return Profile(self.sim_name, f'{self.name}_angle', self.vx,
                       np.angle(self.vy))

    @property
    def real(self):
        return Profile(self.sim_name, f'{self.name}_real', self.vx,
                       self.vy.real)

    @property
    def imag(self):
        return Profile(self.sim_name, f'{self.name}_imag', self.vx,
                       self.vy.imag)

    @property
    def n(self):
        return len(self.vx)

    @property
    def sum(self):
        return np.sum(self.vy)


class Profiler:

    def __init__(self, sim_name, name, linewidth, color):
        self.sim_name = sim_name
        self.name = name
        self.linewidth = linewidth
        self.color = color

    @property
    def line(self):
        raise NotImplementedError()

    def get_profile(self, field, res=64, kind='linear'):
        from scipy.interpolate import interp2d
        f = interp2d(field.x, field.y, field.data.T)
        line = self.line
        length = line.length
        progress = np.linspace(0, 1, ceil(res * length))
        coords = (line.interpolate(p, normalized=True) for p in progress)
        coords = [list(list(p.coords[0])) for p in coords]
        vx = np.asarray([p * length for p in progress])
        vx = np.asarray([y for x, y in coords])
        vy = np.asarray([f(x, y) for x, y in coords])
        return Profile(self.sim_name, self.name, vx, vy)

    def plot(self, ax, *args, **kwargs):
        ax.plot(*self.line.xy, *args, linewidth=self.linewidth,
                color=self.color, **kwargs)


class WaveguideProfiler(Profiler):

    def __init__(self, sim_name, name, waveguide, percent, margin=0,
                 linewidth=1, color='green'):
        Profiler.__init__(self, sim_name, name, linewidth, color)
        self.waveguide = waveguide
        self.percent = percent
        self.margin = margin

    @property
    def line(self):
        p = self.waveguide.coordinates_at(self.percent)
        a = self.waveguide.angle_at(self.percent) - 0.5 * np.pi
        w = self.waveguide.linewidth + 2 * self.margin
        d = np.asarray([np.cos(a), np.sin(a)]) * w / 2
        return geom.asLineString([p + d, p - d])
