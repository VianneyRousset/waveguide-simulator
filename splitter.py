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


def prepare_geom(linewidth, radius, length, spacing, pml):
    from simutils.design.shape import World, LineWaveguide, ArcWaveguide
    from simutils.design.space import attach

    spacing += linewidth
    a = 1
    top_y = a + radius + spacing + radius + a

    top = attach(
        LineWaveguide(linewidth, 0, -a, rpos=[0, top_y]),
        ArcWaveguide(linewidth, radius, -180, -90),
        LineWaveguide(linewidth, length, 0),
        ArcWaveguide(linewidth, radius, -90, 0),
        LineWaveguide(linewidth, 0, a + pml)
    )

    bottom = attach(
        LineWaveguide(linewidth, 0, 2 * a, rpos=[0, -a]),
        ArcWaveguide(linewidth, radius, 180, 90),
        LineWaveguide(linewidth, length, 0),
        ArcWaveguide(linewidth, radius, 90, 0),
        LineWaveguide(linewidth, 0, -(a + pml))
    )

    world = World(margins=[5, 5, 5, -pml])
    world['top'] = top[0]
    world['bottom'] = bottom[0]
    return world, top, bottom


def plot_field(f, mode='normal', shape=None, prefix=None, regions=[],
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


def simulate(wavelength, linewidth, radius, n_eff, length, spacing,
             resolution, pml):

    duration = 1.5 * (length + np.pi * radius) * n_eff

    # geometry and material
    world, top, bottom = prepare_geom(linewidth=linewidth,
                                      radius=radius,
                                      length=length,
                                      spacing=spacing,
                                      pml=pml)
    dev = Device(world, const.n_SiO2, n_eff)

    # sources
    sim = Simulation('splitter', dev, res=resolution, pml=pml)
    sim.create_waveguide_source(wavelength=wavelength,
                                pos=[0, 0],
                                width=8 * linewidth,
                                direction='+y',
                                comp='x')

    # regions
    regions = [
        WaveguideRegion(sim.name, f'input', bottom[0], 0.90,
                        margin=1, length=1),
        WaveguideRegion(sim.name, f'output1', bottom[-1], 0.10,
                        margin=1, length=1),
        WaveguideRegion(sim.name, f'output2', top[-1], 0.10,
                        margin=1, length=1),
    ]

    # run
    sim.run(duration)

    return sim, regions


def inspect(sim, regions):

    ins = inspector.Inspector(sim)
    matplotlib.rc('text', usetex=False)

    # plot
    try:
        plot_field(ins['eps'], 'normal', sim.dev.shape)
        for i in 'ex ey'.split():
            plot_field(ins[i].abs, 'normal', sim.dev.shape)
            plot_field(ins[i].real, 'symetric', sim.dev.shape)
        plot_field(ins['sx'], 'symetric', sim.dev.shape)
        plot_field(ins['sy'], 'symetric', sim.dev.shape)
    except BaseException as e:
        print(f'Plotting failed: {e}')

    # save
    try:
        for f in 'eps ex ey sx sy pwr'.split():
            np.save(f'results/splitter/{f}.npy', ins[f].data)
    except BaseException as e:
        print(f'Saving failed: {e}')

    # splitting
    try:
        f = ins['pwr']
        reg_in1, reg_out1, reg_out2 = regions
        pwr_in1 = reg_in1.cut(f).density
        pwr_out1 = reg_out1.cut(f).density
        pwr_out2 = reg_out2.cut(f).density
        plot_field(ins['pwr'], 'normal', sim.dev.shape, regions=regions)
    except BaseException as e:
        print(f'Plotting pwr failed: {e}')

    return pwr_in1, pwr_out1, pwr_out2


def compute_splitting(wavelength, linewidth, radius, n_eff, length, spacing,
                      resolution, pml):

    sim, regions = simulate(wavelength, linewidth, radius, n_eff, length,
                            spacing, resolution, pml)
    pwr_in1, pwr_out1, pwr_out2 = inspect(sim, regions)

    print(f'({pwr_in1:5.2f}, 0) -> ({pwr_out1:5.2f} {pwr_out2:5.2f})')
    return pwr_out2 / pwr_out1


if __name__ == '__main__':

    param = {
        'wavelength': 1.55,
        'linewidth': 1,
        'radius': 50,
        'n_eff': 1.5842,
        'length': 20,
        'spacing': 1,
        'resolution': 8,
        'pml': 4,
    }

    splitting = compute_splitting(**param)
    print(f'>>> {splitting:05.2f}')
