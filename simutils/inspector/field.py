from .plot import Plot
import numpy as np

class Field(Plot):

    def __init__(self, sim_name, name, data, extent, **kwargs):
        Plot.__init__(self, sim_name, name, **kwargs)
        self.data = data
        self.extent = extent

    @property
    def abs(self):
        return Field(self.sim_name, f'{self.name}_abs', np.abs(self.data),
                     self.extent)

    @property
    def angle(self):
        return Field(self.sim_name, f'{self.name}_angle', np.angle(self.data),
                     self.extent)

    @property
    def real(self):
        return Field(self.sim_name, f'{self.name}_real', self.data.real,
                     self.extent)

    @property
    def imag(self):
        return Field(self.sim_name, f'{self.name}_imag', self.data.imag,
                     self.extent)

    @property
    def square(self):
        return Field(self.sim_name, f'{self.name}_square',
                     np.square(self.data), self.extent)

    @property
    def sum(self):
        return np.sum(self.data)

    @property
    def density(self):
        return self.sum / self.area

    @property
    def Dx(self):
        xmin, xmax, _, _ = self.extent
        return xmax - xmin

    @property
    def Dy(self):
        _, _, ymin, ymax = self.extent
        return ymax - ymin

    @property
    def area(self):
        return self.Dx * self.Dy

    @property
    def nx(self):
        return self.data.shape[0]

    @property
    def ny(self):
        return self.data.shape[1]

    @property
    def x(self):
        return np.linspace(*self.extent[0:2], self.nx)

    @property
    def y(self):
        return np.linspace(*self.extent[2:4], self.ny)

    def get_plot_opts(self, mode='normal'):
        if mode == 'normal':
            return {'cmap': 'viridis'}
        elif mode == 'symetric':
            vmax = max(np.max(self.data), -np.max(self.data))
            return {'cmap': 'seismic', 'vmax': vmax, 'vmin': -vmax}
        elif mode == 'angle':
            return {'cmap': 'viridis', 'vmax': 0, 'vmin': 2 * np.pi}
        else:
            ValueError(f'no mode named: {mode}')

    def plot(self, ax, mode='normal', colorbar=True, **kwargs):
        plot_opts = self.get_plot_opts(mode)
        ax.imshow(self.data.T, origin='lower', extent=self.extent,
                  **plot_opts, **kwargs)
        # colorbar and plt.colorbar()


class Region:

    def __init__(self, sim_name, name):
        self.sim_name = sim_name
        self.name = name

    @property
    def extent(self):
        raise NotImplementedError()

    def cut(self, field):
        xmin, xmax, ymin, ymax = self.extent
        xmask = np.logical_and(xmin < field.x, field.x < xmax)
        ixmin, ixmax = np.argwhere(xmask).reshape(-1)[[0, -1]]
        ymask = np.logical_and(ymin < field.y, field.y < ymax)
        iymin, iymax = np.argwhere(ymask).reshape(-1)[[0, -1]]
        xmin, xmax = field.x[[ixmin, ixmax]]
        ymin, ymax = field.y[[iymin, iymax]]
        return Field(field.sim_name, f'{field.name}_region',
                     data=field.data[ixmin:ixmax, iymin:iymax],
                     extent=[xmin, xmax, ymin, ymax])

    def plot(self, ax, *args, **kwargs):
        xmin, xmax, ymin, ymax = self.extent
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 0.2
        if 'color' not in kwargs:
            kwargs['color'] = 'cyan'
        ax.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin],
                *args, **kwargs)


class WaveguideRegion(Region):

    def __init__(self, sim_name, name, waveguide, percent, margin, length):
        Region.__init__(self, sim_name, name)
        self.waveguide = waveguide
        self.percent = percent
        self.margin = margin
        self.length = length
        self.get_size()

    def get_size(self):
        w, l = self.width, self.length
        a = self.waveguide.angle_at(self.percent) - 0.5 * np.pi
        a = int(np.round(np.degrees(a))) % 360
        if a % 90 > 1e-3:
            raise ValueError('Region has to be on orthogonal waveguides, ',
                             f'got: {a}')
        return [w, l] if a in {0, 180} else [l, w]

    @property
    def pos(self):
        return self.waveguide.coordinates_at(self.percent)

    @property
    def size(self):
        return self.get_size()

    @property
    def width(self):
        return self.waveguide.linewidth + 2 * self.margin

    @property
    def extent(self):
        x, y = self.pos
        dx, dy = self.size
        return [x - dx, x + dx, y - dy, y + dy]

