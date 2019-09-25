#!/usr/bin/env python

from simutils import Simulation, Device
from simutils import inspector
from simutils.inspector.profiler import WaveguideProfiler
from simutils.inspector.field import WaveguideRegion
import simutils.constants as const
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from sys import exit
import json
from subprocess import run
import meep as mp


def prepare_geom(linewidth, length, pml):
    from simutils.design.shape import LineWaveguide, World

    a = 1
    l = length + pml + a
    wg = LineWaveguide(linewidth, 0, 2*l, rpos=[0, -l])

    world = World(margins=[5, 5, -pml, -pml])
    world['wg'] = wg
    return world, wg


def plot_field(f, mode='normal', shape=None, prefix='', regions=[],
               title=None):
    fig, ax = plt.subplots()
    if title:
        ax.set_title(title)
    f.plot(ax, mode)
    if shape:
        shape.plot(ax)
    if prefix is None:
        filepath = f.filepath
    for reg in regions:
        reg.plot(ax)
    else:
        filepath = f.directory / f'{prefix}{f.name}.pdf'
    fig.savefig(filepath)
    print(f'Image saved to {filepath}')
    plt.close()


def simulate(wavelength, linewidth, n_eff, length, resolution, pml):

    duration = 1.5 * length * n_eff

    # geometry and material
    world, wg = prepare_geom(linewidth, length, pml)
    dev = Device(world, const.n_SiO2, n_eff)

    # sources
    sim = Simulation('test', dev, res=resolution, pml=pml)
    sim.create_waveguide_source(wavelength=wavelength,
                                pos=[0, 0],
                                width=8 * linewidth,
                                direction='+x',
                                comp='y')

    # run
    sim.run(duration)

    return sim


def inspect(sim):

    ins = inspector.Inspector(sim)
    matplotlib.rc('text', usetex=False)

    # plot
    try:
        plot_field(ins['eps'], 'normal', sim.dev.shape)
        for i in 'ex ey ez'.split():
            plot_field(ins[i].abs, 'normal', sim.dev.shape)
    except BaseException as e:
        print(f'Plotting failed: {e}')

if __name__ == '__main__':

    wavelength = 1.55
    linewidth = 1
    n_eff = 1.5842
    length = 3
    resolution = 16
    pml = 4

    sim = simulate(wavelength, linewidth, n_eff, length,
                   resolution, pml)
    inspect(sim)
