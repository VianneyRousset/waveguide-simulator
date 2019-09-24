from .polygons import Polygons
from .field import Field

__all__ = ['profiler']


class Inspector:

    def __init__(self, sim):
        self.sim = sim

    @property
    def shape(self):
        return Polygons(self.sim.name, 'shape', self.sim.dev.shape.polygons)

    def __getitem__(self, k):
        data, extent = self.sim[k]
        return Field(self.sim.name, k, data, extent)
