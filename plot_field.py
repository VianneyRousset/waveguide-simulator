#!/usr/bin/env python

from optparse import OptionParser
from dataio import load_field, StaticFieldData, DynamicFieldData
import numpy as np
import matplotlib.pyplot as plt
import sys


def parse_args():
    parser = OptionParser(usage='%prog [-v VALUE] [-e EXTENT] -f FIELD INPUT '
                          'OUTPUT')

    parser.add_option('-v', '--value', metavar='VALUE',
                      action='store', dest='value', default='real',
                      help='plot real, imag, angle or abs')

    parser.add_option('-e', '--extent', metavar='EXTENT',
                      action='store', dest='extent', default=None,
                      help='set data coordinate, format: xmin,xmax,ymin,xmax')

    parser.add_option('-f', '--field', metavar='FIELD',
                      action='store', dest='field', default=None,
                      help='field to read, use "real,imag" for complex fields')
    (options, args) = parser.parse_args()

    if len(args) != 2:
        print('Please provide exactly one input path and one output',
              file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if options.field and len(options.value.split(',')) not in {1, 2}:
        print('Please provide one or two field(s) to read', file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    return args[0], args[1], options


def convert_field(field, value):
    if value == 'real':
        return field.real
    elif value == 'imag':
        return field.imag
    elif value == 'angle':
        return np.angle(field)
    elif value == 'abs':
        return np.abs(field)
    else:
        raise KeyError(value)


def plot(field, value, output_path, vmin=None, vmax=None, extent=None):
    fig, ax = plt.subplots()
    field = convert_field(field, value)
    ax.imshow(field.T, origin='lower', extent=extent, vmin=vmin, vmax=vmax)
    fig.savefig(output_path)
    plt.close(fig)


if __name__ == '__main__':

    input_path, output_path, args = parse_args()
    names = args.field.split(',')
    field = load_field(input_path, *names)

    if isinstance(field, StaticFieldData):
        plot(field.data, args.value, f'{output_path}.png', args.extent)
    elif isinstance(field, DynamicFieldData):
        N = len(field)
        print('Finding min and max...')
        if args.value == 'angle':
            vmin, vmax = 0, 2 * np.pi
        else:
            vmin, vmax = field.get_min(args.value), field.get_max(args.value)
        print('Plotting...')
        for n, f in enumerate(field):
            print(f'{round(n/N * 100)}%')
            plot(f, args.value, f'{output_path}_{n:05d}.png', vmin, vmax,
                 args.extent)
    else:
        raise SystemError()


'''
        class Inspector:

            def __init__(self, geometry, simulation, margins):
                self.geom = geometry
        self.sim = simulation
        self.margins = margins

    def multi_plots(self, keys, with_shape_lines=True):
        fig, axes = plt.subplots(ceil(len(keys)/2), min(len(keys), 2))
        for k, ax in zip(keys, axes.reshape(-1)):
            self.plot(k, ax)
        return fig

    def plot(self, key, ax):
        if not isinstance(key, str):
            for k in key:
                self.plot(self, k, ax)
            return

        available_plots = {'geometry', 'epsilon', 'amplitude', 'phase'}
        if key not in available_plots:
            raise KeyError(k)

        extent = self.get_extent()

        def plot(field):
            ax.imshow(field, origin='lower', extent=extent)

        if (key == 'geometry'):
            self.geom.plot(ax)
        elif (key == 'epsilon'):
            plot(self.get_epsilon().T)
        elif (key == 'amplitude'):
            plot(self.get_amplitude())
        elif (key == 'phase'):
            plot(self.get_phase())
        ax.set_title(key)

    def get_extent(self):
        xmin, xmax, ymin, ymax = self.geom.boundaries
        t = self.margins
        return [xmin - t, xmax + t, ymin - t, ymax + t]

    def get_epsilon(self):
        return self.sim.get_epsilon()

    def get_field(self):
        return self.sim.get_array(center=mp.Vector3(),
                           size=mp.Vector3(*self.geom.size),
                           component=mp.Ez, cmplx=True).T

    def get_amplitude(self):
        return np.abs(self.get_field())

    def get_phase(self):
        return np.angle(self.get_field())
'''
