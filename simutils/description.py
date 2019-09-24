#!/usr/bin/env python

import meep as mp
# from . import dataio
from pathlib import Path


class Summary(dict):

    def __str__(self):
        import json

        def convert(x):
            if isinstance(x, Description):
                return x.summary
            if isinstance(x, set) or isinstance(x, mp.Vector3):
                return list(x)
            raise ValueError(f'Unable to convert type: {type(x)}')

        return json.dumps(self, indent=4, sort_keys=True,
                          default=convert)

    def save(self, filepath):
        with open(filepath, mode='w') as f:
            f.write(str(self))


class Description:

    @property
    def summary(self):
        return Summary(self.__dict__)

    def save(self, filepath):
        import pickle
        with open(filepath, mode='wb') as f:
            pickle.dump(self, f)


class FieldDescription(Description):

    def __init__(self, filepath, real, imag=None):
        self.filepath = filepath
        self.real = real
        self.imag = imag

    def load(self):
        return dataio.load_field(self.filepath, self.real, self.imag)


class Material(mp.Medium, Description):

    def __init__(self, name, *args, **kwargs):
        self.name = name
        mp.Medium.__init__(self, *args, **kwargs)

    @property
    def summary(self):
        return Summary({
            'name': self.name,
            'eps': self.epsilon(0)[0, 0],
            'mu': self.mu(0)[0, 0]})


class Source(mp.Source, Description):

    def __init__(self, fq, pos, size, component):
        self.pos = pos
        self.size = size
        self.comp = component
        mp.Source.__init__(self, mp.ContinuousSource(fq),
                           component=self._get_meep_comp(self.comp),
                           center=mp.Vector3(*self.pos),
                           size=mp.Vector3(*self.size))

    @property
    def summary(self):
        return Summary({
            'frequency': self.src.frequency,
            'position': self.pos,
            'size': self.size,
            'component': self.comp,
        })

    def _get_meep_comp(self, name):
        return {'ex': mp.Ex, 'ey': mp.Ey, 'ez': mp.Ez}[name]


class SimulationDescription(Description):

    REQUIRED = {'name', 'pml_thickness', 'resolution', 'duration', 'dt',
                'default_material', 'shape_material', 'shape', 'done',
                'sources', 'static_fields', 'dynamic_fields',
                'shape_filepath'}

    def __init__(self, name):
        self.name = name

    @property
    def size(self):
        return list(self.shape.size)

    @property
    def sim_size(self):
        return [x + 2 * self.pml_thickness for x in self.size]

    def get_sim_args(self):
        return {
            'cell_size': mp.Vector3(*self.sim_size),
            'geometry_center': mp.Vector3(*self.size) / 2,
            'resolution': self.resolution,
            'sources': self._create_sources(),
            'boundary_layers': self._create_boundary_layers,
            'force_complex_fields': self._create_material_func(),
            'filename_prefix': f'{self.name}-',
        }

    def create_simulation(self):
        return mp.Simulation(**self.sim_args)

    def get_step_funcs(self):
        sf = {self._get_meep_output_comp(f) for f in self.static_fields}
        sf = (mp.at_beginning(f) for f in sf)
        df = {self._get_meep_output_comp(f) for f in self.dynamic_fields}
        df = (mp.at_every(self.dt, f) for f in df)
        volume = mp.in_volume(mp.Volume(center=self.size / 2, size=self.size))
        fields = mp.with_prefix(f'{results_dir}/', *sf, *df)
        return volume, fields

    def _create_sources(self):
        return [src.create_source() for src in self.sources]

    def _create_boundary_layers(self):
        return [mp.PML(self.pml_thickness)]

    def _create_material_func(self):
        mat1 = self._get_meep_material(self.default_material)
        mat2 = self._get_meep_material(self.shape_material)

        def material(p):
            x, y, z = p
            return mat2 if self.world.contains([x, y]) else mat1

        return material

    def _get_meep_material(self, name):
        return {'vacuum': mp.Medium(epsilon=1),
                'silica': mp.Medium(epsilon=3.9)}[name]

    def _get_meep_output_comp(self, name):
        return {
            'b': mp.output_bfield, 'd': mp.output_dfield,
            'e': mp.output_efield, 'eps': mp.output_epsilon,
            'h': mp.output_hfield, 'mu': mp.output_mu,
            'poynting': mp.output_poynting, 's': mp.output_sfield,
            'power': mp.output_tot_pwr,
        }[name]

    def _load_geom(self, **kwargs):
        if 'world' in kwargs:
            self.extent = kwargs['world'].extent
        self._ensure_set({'extent'}, **kwargs)

    def _ensure_set(self, keys, **kwargs):
        for k in keys:
            try:
                self.__setattr__(k, kwargs[k])
            except KeyError:
                pass
        available = set(self.__dict__).union({'shape_filepath', 'filepath'})
        missing = keys.difference(available)
        if missing:
            raise ValueError(f'Missing description element: "{missing}"')




    REQUIRED = {'name', 'pml_thickness', 'resolution', 'duration', 'dt',
                'default_material', 'shape_material', 'shape', 'done',
                'sources', 'static_fields', 'dynamic_fields',
                'shape_filepath'}



class Simulation:

    RESULTS_DIR = Path('results')

    def __init__(self, name):
        self.name = name

    @property
    def summary_filepath(self):
        return self.RESULTS_DIR / f'{self.name}_summary.json'

    @property
    def desc_filepath(self):
        return self.RESULTS_DIR / f'{self.name}_desc.pickle'

    def autoload(self, filepath=None):
        self.desc = pickle.load(filepath or self.filepath)

    def init(self, desc):
        self.desc = desc

    def run(self):
        from subprocess import run, PIPE
        self.desc.save()
        print(self.desc)
        sim = run(['mpirun', '-np', '4', './simulation.py',
                   self.desc.filepath],
                  encoding='utf-8')
        if sim.returncode != 0:
            raise SystemError('Simulation failed')






'''
    def load(self, filepath=None):
        with open(filepath or self.filepath, mode='r') as f:
            print(json.loads(f.read()), flush=True)
            self._ensure_set({}, **json.loads(f.read()))
        with open(self.shape_filepath, mode='rb') as f:
            self.shape = pickle.load(f)
        self._ensure_set(self.REQUIRED)
'''


