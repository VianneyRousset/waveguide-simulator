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


def prepare_geom(linewidth, radius, pml):
    from simutils.design.shape import World, LineWaveguide, ArcWaveguide
    from simutils.design.space import attach
    line1, arc, line2 = attach(
        LineWaveguide(linewidth, 4.5, 0, rpos=[-0.5, 0]),
        ArcWaveguide(linewidth, radius, -90, 0),
        LineWaveguide(linewidth, 0, 6 + pml),
    )
    world = World(margins=[5, 5, 5, -pml])
    world['waveguide'] = line1
    return world, line1, arc, line2


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


def simulate(n_eff, width, radius, resolution, wavelength, pml):

    duration = 1.5 * (4.5 + 0.5 * np.pi * radius + 5) * n_eff

    # geometry and material
    world, line1, arc, line2 = prepare_geom(width, radius, pml)
    dev = Device(world, const.n_SiO2, n_eff)

    # sources
    sim = Simulation('curve', dev, res=resolution, pml=pml)
    sim.create_waveguide_source(wavelength=wavelength,
                                pos=[0, 0],
                                width=8)

    # regions
    regions = [
        WaveguideRegion(sim.name, f'region_start', line1, 0.80,
                        margin=1, length=1),
        WaveguideRegion(sim.name, f'region_end', line2, 0.25,
                        margin=1, length=1),
    ]

    # run
    sim.run(duration)

    return sim, regions


def inspect(sim, regions, name):

    ins = inspector.Inspector(sim)
    matplotlib.rc('text', usetex=False)

    # plot
    prefix = f'{name}_'
    plot_field(ins['eps'], 'normal', sim.dev.shape, prefix)
    for i in 'ex ey'.split():
        plot_field(ins[i].abs, 'normal', sim.dev.shape, prefix)
        plot_field(ins[i].real, 'symetric', sim.dev.shape, prefix)
    plot_field(ins['sx'], 'symetric', sim.dev.shape, prefix)
    plot_field(ins['sy'], 'symetric', sim.dev.shape, prefix)

    # losses
    f = ins['pwr']
    reg_start, reg_end = regions
    insertion_losses = reg_end.cut(f).density / reg_start.cut(f).density
    plot_field(ins['pwr'], 'symetric', sim.dev.shape, prefix, regions,
               title=fr'Insertion: ${insertion_losses*100:05.2f}\%$')

    return insertion_losses


def compute_insertion(thickness, n_eff, width, radius, resolution,
                      wavelength, pml):
    name = f'h{thickness:04.3f}_w{width:04.2f}_n{n_eff:04.3}_r{radius:05.2f}'
    print(f'Computing insertion for {name}')
    sim, regions = simulate(n_eff, width, radius, resolution, wavelength, pml)
    losses = inspect(sim, regions, name)
    print(f'r{name} >>> {losses*100:05.2f}%')
    return losses


if __name__ == '__main__':

    param = {
        'resolution': 2,
        'wavelength': 1.55,
        'pml': 4,
    }

    radii = [50, 100, 200, 500]

    #              thickness    n_eff       width
    waveguides = [
                  (0.200,       1.5842,     0.80),
                  (0.200,       1.5842,     1.00),
                  (0.200,       1.5842,     1.10),
                  (0.500,       1.7952,     0.50),
                  (0.500,       1.7952,     0.60),
                  (0.500,       1.7952,     0.75)]

    insertions = [[compute_insertion(h, n, w, r, **param)
                   for r in radii]
                  for h, n, w in waveguides]

    # save summary
    waveguides = np.asarray(waveguides).T
    summary = {
        'radii': radii,
        'waveguides_width': list(waveguides[0]),
        'waveguides_thickness': list(waveguides[1]),
        'waveguides_n_eff': list(waveguides[2]),
        'insertions': insertions,
    }
    summary.update(param)
    print(json.dumps(summary, indent=4, sort_keys=True))
    with open('results/curve/summary.json', mode='w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)
