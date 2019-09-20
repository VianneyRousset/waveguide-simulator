#!/usr/bin/env python

from argparse import ArgumentParser
from description import SimulationDescription


def parse_args():
    parser = ArgumentParser(description='Run meep simulation')
    parser.add_argument('inputfile', metavar='SIM_DESC', type=str, nargs=1,
                        help='simulation description file')

    return parser.parse_args().inputfile[0]


if __name__ == '__main__':
    desc = SimulationDescription(parse_args())
    print(desc)
    sim = desc.create_simulation()
    sim.run(*desc.get_step_funcs(), until=desc.duration)
    #mp.all_wait()
