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
