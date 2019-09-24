from .plot import Plot


class Region(Plot):

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


