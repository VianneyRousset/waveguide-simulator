#!/usr/bin/env python

from simutils import Simulation, Device
from simutils import inspector
from simutils.inspector.profiler import WaveguideProfiler
from matplotlib import pyplot as plt
import numpy as np


def prepare_geom(linewidth, radius):
    from simutils.design.shape import World, LineWaveguide, ArcWaveguide
    from simutils.design.space import attach
    line1, arc, line2 = attach(
        LineWaveguide(linewidth, 15, rpos=[-5, 0]),
        ArcWaveguide(linewidth, radius, np.radians(-90), np.radians(90)),
        LineWaveguide(linewidth, -15)
    )
    world = World()
    world['waveguide'] = line1
    return world, line1, arc, line2


if __name__ == '__main__':

    wavelength = 1.55
    linewidth = 0.6
    radius = 4
    comp = 'ez'
    duration = 100

    # SIMULATION

    # geometry and material
    world, _, arc, _ = prepare_geom(linewidth, radius)
    dev = Device(world, 'vacuum', 'siliconNitride')

    # saving shape

    # sources
    sim = Simulation('propagation_losses', dev, res=20)

    sim.create_waveguide_source(wavelength=wavelength,
                                pos=[0, 0],
                                width=8,
                                comp=comp)
    sim.run(duration, dt=1)

    # INSPECTION
    ins = inspector.Inspector(sim)

    # profilers
    print('Computing intensities...')
    N_prof = 32
    wg_percents = np.linspace(0, 1, N_prof)
    wg_length = 2 * np.pi * radius
    profilers = [WaveguideProfiler(sim.name, f'wg_profile_{n}', arc, n,
                                   margin=4)
                 for n in wg_percents]

    f = ins[comp].abs
    comp_abs_sum = [p.get_profile(f).sum for p in profilers]
    f = ins['p'].abs
    p_abs_sum = [p.get_profile(f).sum for p in profilers]

    # epsilon
    f = ins['eps']
    fig, ax = plt.subplots()
    f.plot(ax, 'normal')
    ax.set_title(f.name.replace('_', '\_'))
    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')
    plt.close()

    # comp real
    f = ins[comp].real
    fig, ax = plt.subplots()
    f.plot(ax, 'symetric')
    ins.shape.plot(ax)
    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')

    # intensities
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    f = ins[comp].abs
    f.plot(ax1, 'normal')
    ins.shape.plot(ax1)
    for p in profilers:
        p.plot(ax1)

    f = ins['p'].abs
    f.plot(ax2, 'normal')
    ins.shape.plot(ax2)
    for p in profilers:
        p.plot(ax2)

    ax3.plot(wg_percents * wg_length, comp_abs_sum)
    ax4.plot(wg_percents * wg_length, p_abs_sum)

    filepath = str(f.directory / 'comp_p.pdf')
    fig.savefig(filepath)
    print(f'Image saved to {filepath}')

    plt.close()
    print('End')
