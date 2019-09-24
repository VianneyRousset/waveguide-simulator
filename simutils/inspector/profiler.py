import numpy as np
from math import ceil
from shapely import geometry as geom
from .plot import Plot
import meep as mp


def si_label(name, unit, delimiters=(r'\left[',r'\right]')):
    dl, dr = ('', '') if delimiters is None else delimiters
    al, ar = '{', '}'
    return fr'{name}$\;{dl}\,\si{al}{unit}{ar}\,{dr}$'


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


class SpectrumProfile(Profile):

    def plot(self, ax, autoxlabel=True, *args, **kwargs):
        ax.plot(self.vx, self.vy, *args, **kwargs)
        if autoxlabel:
            ax.set_xlabel(si_label('Frequencies', r'\hertz'))


class Profiler:

    def __init__(self, sim_name, name, linewidth, color):
        self.sim_name = sim_name
        self.name = name
        self.linewidth = linewidth
        self.color = color
        self.flux = None

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

    @property
    def flux_spectrum(self):
        return SpectrumProfile(self.sim_name, f'{self.name}_flux',
                               mp.get_flux_freqs(self.flux),
                               mp.get_fluxes(self.flux))

    def plot(self, ax, *args, **kwargs):
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = self.linewidth
        if 'color' not in kwargs:
            kwargs['color'] = self.color
        ax.plot(*self.line.xy, *args, **kwargs)

    @property
    def get_flux_region(self):
        raise NotImplementedError()


class WaveguideProfiler(Profiler):

    def __init__(self, sim_name, name, waveguide, percent, margin=0,
                 linewidth=1, color='green'):
        Profiler.__init__(self, sim_name, name, linewidth, color)
        self.waveguide = waveguide
        self.percent = percent
        self.margin = margin

    @property
    def pos(self):
        return self.waveguide.coordinates_at(self.percent)

    @property
    def angle(self):
        return self.waveguide.angle_at(self.percent) - 0.5 * np.pi

    @property
    def width(self):
        return self.waveguide.linewidth + 2 * self.margin

    @property
    def line(self):
        p, a, w = self.pos, self.angle, self.width
        d = np.asarray([np.cos(a), np.sin(a)]) * w / 2
        return geom.asLineString([p + d, p - d])

    def get_flux_region(self):
        p, a, w = self.pos, int(np.round(self.angle)) % 360, self.width
        if a % 90 < 1e-3:
            raise ValueError('Flux can only be measured on orthogonal lines'
                             f'got: {a}')
        s = [w, 0] if a in {0, 180} else [0, w]
        return mp.FluxRegion(center=mp.Vector3(*p), size=mp.Vector3(*s))
