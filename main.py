#!/usr/bin/env python

import shaper
import numpy as np
from math import ceil, floor
import matplotlib.pyplot as plt
from description import Simulation, SourceDescription

#from scipy.optimize import minimize
#from numpy.polynomial.polynomial import Polynomial

if __name__ == '__main__':

    # geometry and material
    geom = shaper.attach(
        shaper.LineWaveguide(8),
        shaper.SlineWaveguide(2, 2),
        shaper.MZIWaveguide(2, 2, 1),
    )

    # source
    sources = [SourceDescription(wavelength=1.55,
                                 position=[0, 0],
                                 size=[0, shaper.default_linewidth, 1e20],
                                 component='ez')]

    # simulation
    sim = Simulation(name='test',
                     shape=geom[0]._shape,
                     pml_thickness=1,
                     resolution=64,
                     duration=50,
                     dt=0.5,
                     default_material='vacuum',
                     shape_material='silica',
                     sources=sources,
                     static_fields={'eps', 'mu'},
                     dynamic_fields={'e', 'b', 'power'})

    sim.run()
    #mp.quiet()
