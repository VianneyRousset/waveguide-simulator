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


def prepare_geom(linewidth, radius, pml):
    from simutils.design.shape import World, LineWaveguide
    line = LineWaveguide(linewidth, 40, rpos=[-0.5, 0])
    world = World(margins=[4, -pml, 4, 4])
    world['line'] = line
    return world, line


def plot_field(f, mode='normal', shape=None):
    fig, ax = plt.subplots()
    f.plot(ax, mode)
    if shape:
        shape.plot(ax)
    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')
    plt.close()


if __name__ == '__main__':

    resolution = 40
    wavelength = 1.55
    linewidth = 1.0
    radius = 2.5
    comp = 'ez'
    duration = 80
    pml = 8

    # SIMULATION

    # geometry and material
    world, line = prepare_geom(linewidth, radius, pml)
    dev = Device(world, const.n_SiO2, const.n_eff)

    # sources
    sim = Simulation('line', dev, res=resolution, pml=pml)
    sim.create_waveguide_source(wavelength=wavelength,
                                pos=[0, 0],
                                width=8)

    # regions
    regions = {}
    for n, p in enumerate(np.linspace(0.15, 0.6, 16)):
        regions[n] = WaveguideRegion(sim.name, f'region_{n}',
                                     world['line'], p,
                                     margin=1, length=0.4)
    sim.run(duration)

    # INSPECTION
    ins = inspector.Inspector(sim)
    matplotlib.rc('text', usetex=True)

    # fields
    plot_field(ins['eps'], 'normal', shape=dev.shape)
    plot_field(ins['ez'].abs, 'normal', shape=dev.shape)
    plot_field(ins['ez'].real, 'symetric', shape=dev.shape)
    plot_field(ins['sx'], 'symetric', shape=dev.shape)
    plot_field(ins['sy'], 'symetric', shape=dev.shape)

    # fluxes
    fig, axes = plt.subplots(1, 2)
    f = ins['pwr']
    ax0 = plt.subplot2grid((2, 2), (0, 0))
    ax1 = plt.subplot2grid((2, 2), (0, 1))
    ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1)
    sums = []
    f.plot(ax2, 'normal')
    dev.shape.plot(ax2)
    for n in sorted(regions):
        reg = regions[n]
        reg.plot(ax2)
        f_reg = reg.cut(f)
        sums = np.append(sums, f_reg.density)
    print(f'Image saved to {f.filepath}')
    ax0.plot(sums)
    fig.savefig(f.filepath)

    plt.close()
