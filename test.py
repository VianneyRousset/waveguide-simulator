#!/usr/bin/env python

from simutils import Simulation, Device
from simutils import inspector
from simutils.inspector.profiler import WaveguideProfiler
from matplotlib import pyplot as plt


def prepare_geom():
    from simutils.design.shape import LineWaveguide, SlineWaveguide
    from simutils.design.space import attach
    linewidth = 1
    return attach(
        LineWaveguide(linewidth, 20, rpos=[-10, 0]),
        SlineWaveguide(linewidth, 30, 10),
        LineWaveguide(linewidth, 10))


if __name__ == '__main__':

    # geometry and material
    from simutils.design.shape import World
    geom = prepare_geom()
    world = World()
    world['waveguide'] = geom[0]
    dev = Device(world, 'vacuum', 'silica')

    sim = Simulation('test', dev, res=10)

    sim.create_waveguide_source(wavelength=1.55,
                                pos=[0, 0],
                                width=8,
                                comp='ez')
    sim.run(120, dt=1)

    # inspect
    ins = inspector.Inspector(sim)
    prof = WaveguideProfiler(sim.name, 'wg_profile', geom[0], 0.6, margin=2)

    # epsilon
    f = ins['eps']
    fig, ax = plt.subplots()
    f.plot(ax, 'normal')
    prof.plot(ax)
    ax.set_title(f.name.replace('_', '\_'))
    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')
    plt.close()

    # ez real
    f = ins['ez'].real
    fig, [ax1, ax2] = plt.subplots(1, 2)

    f.plot(ax1, 'symetric')
    prof.plot(ax1)
    ins.shape.plot(ax1)

    prof.get_profile(ins['eps']).plot(ax2.twinx(), color='orange', label='eps')
    prof.get_profile(f).plot(ax2, label=f.name.replace('_', '\_'))

    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')

    # ez abs
    f = ins['ez'].abs
    fig, [ax1, ax2] = plt.subplots(1, 2)

    f.plot(ax1, 'normal')
    prof.plot(ax1)
    ins.shape.plot(ax1)

    prof.get_profile(ins['eps']).plot(ax2.twinx(), color='orange', label='eps')
    prof.get_profile(f).plot(ax2, label=f.name.replace('_', '\_'))

    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')

    plt.close()
    print('End')
