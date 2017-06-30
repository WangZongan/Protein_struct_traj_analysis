#!/usr/bin/env python
import sys
import cPickle as cp
import numpy as np
import argparse

sys.path.insert(1, '/home/zonganw/upside_devel/upside-devel/py/')
import upside_engine as ue

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help='Config file to modify with target information')
    parser.add_argument('output', help='Output path')
    args = parser.parse_args()

    engine = ue.Upside(args.config)
    ini_pos = engine.initial_pos
    energy  = engine.energy(ini_pos)
    env_cov = engine.get_output('environment_coverage')

    with open(args.output, 'w') as f:
        cp.dump(env_cov, f, -1)


if __name__ == '__main__':
    main()
