#!/usr/bin/env python

import meep as mp
from meep import Vector3
from . import dataio


class SimulationOutputsBase(dict):

    def __init__(self, name):
        dict.__init__(self)
        self.name = name
        self.directory = self.create_output_dir(self.RESULTS_DIR / name)

    @property
    def static_fields(self):
        return {k: v for k, v in self.items()
                if type(v) is dataio.StaticFieldData}

    @property
    def dynamic_fields(self):
        return {k: v for k, v in self.items()
                if type(v) is dataio.DynamicFieldData}

    def create_output_dir(self, path):
        path.mkdir(parents=True, exist_ok=True)
        return path

    def get_meep_output_comp(self, name):
        return {
            'b': mp.output_bfield, 'bx': mp.output_bfield_x,
            'by': mp.output_bfield_y, 'bz': mp.output_bfield_z,
            'e': mp.output_efield, 'ex': mp.output_efield_x,
            'ey': mp.output_efield_y, 'ez': mp.output_efield_z,
            'd': mp.output_dfield, 'dx': mp.output_dfield_x,
            'dy': mp.output_dfield_y, 'dz': mp.output_dfield_z,
            'h': mp.output_hfield, 'hx': mp.output_hfield_x,
            'hy': mp.output_hfield_y, 'hz': mp.output_hfield_z,
            's': mp.output_sfield, 'sx': mp.output_sfield_x,
            'sy': mp.output_sfield_y, 'sz': mp.output_sfield_z,
            'eps': mp.output_epsilon, 'mu': mp.output_mu,
            'power': mp.output_tot_pwr,
        }[name]


class SimulationBase:

    def __init__(self, name, dev, pml):
        self.name = name
        self.dev = dev
        self.pml = pml

    @property
    def pmlv(self):
        return Vector3(self.pml, self.pml)

    @property
    def extent(self):
        return self.dev.extent

    @property
    def size(self):
        return self.dev.size

    @property
    def center(self):
        return self.dev.center

    @property
    def cell_size(self):
        return self.size + 2 * self.pmlv

    def get_meep_comp(self, name):
        return {
            'ex': mp.Ex, 'ey': mp.Ey, 'ez': mp.Ez, 'ep': mp.Ep, 'er': mp.Er,
            'bx': mp.Bx, 'by': mp.By, 'bz': mp.Bz, 'bp': mp.Bp, 'br': mp.Br,
            'dx': mp.Dx, 'dy': mp.Dy, 'dz': mp.Dz, 'dp': mp.Dp, 'dr': mp.Dr,
            'hx': mp.Hx, 'hy': mp.Hy, 'hz': mp.Hz, 'hp': mp.Hp, 'hr': mp.Hr,
            'sx': mp.Sx, 'sy': mp.Sy, 'sz': mp.Sz, 'sp': mp.Sp, 'sr': mp.Sr,
            'p': mp.P, 'eps': mp.Dielectric, 'mu': mp.Permeability,
        }[name.lower()]
