#!/usr/bin/env python

import meep as mp
from meep import Vector3
from . import dataio
from pathlib import Path
from .base import SimulationBase, SimulationOutputsBase


class Device:

    MATERIALS = {
        'vacuum': mp.Medium(epsilon=1),
        'silica': mp.Medium(epsilon=3.9),
        'siliconNitride': mp.Medium(epsilon=7.5)
    }

    def __init__(self, shape, default_material, shape_material):
        self.shape = shape
        self.polygons = shape.polygons
        self.default_material = self.MATERIALS[default_material]
        self.shape_material = self.MATERIALS[shape_material]

    def at(self, p):
        if self.polygons.contains(p):
            return self.shape_material
        else:
            return self.default_material

    @property
    def geometry(self):
        return [mp.Prism([Vector3(*xy) for xy in p.exterior.coords[:-1]],
                         mp.inf, material=self.shape_material)
                for p in self.polygons]

    @property
    def extent(self):
        return self.shape.extent

    @property
    def size(self):
        return Vector3(*self.shape.size)

    @property
    def center(self):
        return Vector3(*self.shape.center)

    @property
    def volume(self):
        return mp.Volume(center=self.center, size=self.size)


class Simulation(SimulationBase):

    def __init__(self, name, dev, res=32, pml=1):
        SimulationBase.__init__(self, name, dev, pml)
        self.sources = []
        self.resolution = res

    def set_resolution(self, res):
        self.resolution = res

    def create_waveguide_source(self, wavelength, pos, width, comp,
                                direction=[1, 0]):
        from numpy import pi
        k = mp.Vector3(*direction).unit() * 2 * pi / wavelength
        s = mp.EigenModeSource(src=mp.ContinuousSource(wavelength=wavelength),
                               center=mp.Vector3(*pos),
                               size=mp.Vector3(y=width),
                               direction=mp.NO_DIRECTION,
                               eig_kpoint=mp.Vector3(0.4),
                               eig_band=1,
                               eig_parity=mp.EVEN_Y+mp.ODD_Z,
                               eig_match_freq=True,
                               component=self.get_meep_comp(comp))
        self.sources.append(s)

    def create_source(self, wavelength, pos, size, comp):
        src = mp.Source(mp.ContinuousSource(wavelength=wavelength),
                        center=mp.Vector3(*pos),
                        size=mp.Vector3(*size),
                        component=self.get_meep_comp(comp))
        self.sources.append(src)

    def run(self, duration, dt):
        boundary_layers = self._create_boundary_layers()
        self.sim = mp.Simulation(cell_size=self.cell_size,
                                 geometry_center=self.center,
                                 default_material=self.dev.default_material,
                                 geometry=self.dev.geometry,
                                 resolution=self.resolution,
                                 sources=self.sources,
                                 boundary_layers=boundary_layers,
                                 force_complex_fields=True)
        self.sim.run(mp.in_volume(self.dev.volume), until=duration)
        mp.all_wait()

    def __getitem__(self, k):
        return self.sim.get_array(center=self.center, size=self.size,
                                  component=self.get_meep_comp(k))

    def _create_boundary_layers(self):
        return [mp.PML(self.pml)]


class SimulationOutputs(SimulationOutputsBase):

    RESULTS_DIR = Path('results')

    def __init__(self, name):
        SimulationOutputsBase.__init__(self, name)

    def new_output(self, name, cmx, dynamic):
        filepath = self.directory / f'{name}.h5'
        v = f'{name}.r', f'{name}.i' if cmx else name
        if dynamic:
            self[name] = dataio.DynamicFieldData(filepath, *v)
        else:
            self[name] = dataio.StaticFieldData(filepath, *v)

    def get_step_funcs(self, dt):
        sf = (mp.at_beginning(self.get_meep_output_comp(f))
              for f in self.dynamic_fields)
        df = (mp.at_every(dt, self.get_meep_output_comp(f))
              for f in self.static_fields)
        return (*sf, *df)
