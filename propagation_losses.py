#!/usr/bin/env python

from simutils import Simulation, Device
from simutils import inspector
from simutils.inspector.profiler import WaveguideProfiler
from matplotlib import pyplot as plt
import numpy as np
from sys import exit


def prepare_geom(linewidth, radius, pml):
    from simutils.design.shape import World, LineWaveguide, ArcWaveguide
    from simutils.design.space import attach
    line1, arc, line2 = attach(
        LineWaveguide(linewidth, 2, rpos=[-0.5, 0]),
        ArcWaveguide(linewidth, radius, np.radians(-90), np.radians(0)),
        LineWaveguide(linewidth, 0, 10)
    )
    world = World(margins=[4, 4, 4, -pml])
    world['waveguide'] = line1
    return world, line1, arc, line2


def plot_field(f, mode='normal'):
    fig, ax = plt.subplots()
    f.plot(ax, mode)
    title = f.name.replace('_', ' ')
    ax.set_title(title)
    fig.savefig(f.filepath)
    print(f'Image saved to {f.filepath}')
    plt.close()


if __name__ == '__main__':

    resolution = 50
    wavelength = 1.55
    linewidth = 1.0
    radius = 2.5
    comp = 'ez'
    duration = 40
    pml = 1

    # SIMULATION

    # geometry and material
    world, line1, arc, _ = prepare_geom(linewidth, radius, pml)
    dev = Device(world, 'vacuum', 'siliconNitride')

    # flux regions

    import meep as mp
    fluxes = [mp.FluxRegion(mp.Vector3(*line1.coordinates_at(p)),
                            mp.Vector3(0, 2, 0), direction=mp.X)
              for p in [0.2, 0.4, 0.6, 0.8, 1.0]]
    flux = fluxes[0]

    # sources
    sim = Simulation('propagation_losses', dev, res=resolution, pml=pml)
    sim.create_waveguide_source(wavelength=wavelength,
                                pos=[0, 0],
                                width=8)
    # sim.create_source(wavelength, [0, 0], [0, linewidth], 'ez')

    # profilers
    sim.profilers['test'] = WaveguideProfiler(sim.name, 'wg_prof_test', line1,
                                              0.5)
    sim.run(duration, fluxes)

    print('saving fluxes')
    fig, ax = plt.subplots()
    flux = sim.profilers['test'].flux_spectrum
    flux.plot(ax)
    fig.savefig(flux.filepath)

    # INSPECTION
    ins = inspector.Inspector(sim)

    plot_field(ins['eps'], 'normal')

    plot_field(ins['ez'].abs, 'normal')
    plot_field(ins['ez'].real, 'symetric')

    plot_field(ins['sx'], 'symetric')
    plot_field(ins['sy'], 'symetric')

    exit()

    # intensities
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

    f = ins['p'].abs
    f.plot(ax2, 'normal')
    ins.shape.plot(ax2)
    for p in profilers:
        p.plot(ax2)

    ax3.plot(wg_percents * wg_length, comp_abs_integral)
    ax3.plot(wg_percents * wg_length, p_abs_integral)

    # profiles
    p.get_profile(ins['eps']).plot(ax4.twinx(), color='black', linewidth=0.2)
    for n, p in enumerate(profilers):
        p.get_profile(f).plot(ax4, label=str(n))

    filepath = str(f.directory / 'comp_p_profiles.pdf')
    fig.savefig(filepath)
    print(f'Image saved to {filepath}')

    plt.close()
    print('End')
