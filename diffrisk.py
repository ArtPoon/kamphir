"""
Driver script for fitting a differential risk SI model to tree shape.
"""
import os
import scipy.stats as stats
import argparse
from kamphir import Kamphir

parser = argparse.ArgumentParser(description='Fit a differential risk SI model to '
                                             'the shape of a phylogenetic tree using'
                                             'kernel-ABC',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('nwkfile', help='File containing Newick tree string.')
parser.add_argument('logfile', help='File to log ABC-MCMC traces.')
parser.add_argument('-skip', default=10, help='Number of steps in ABC-MCMC to skip for log.')
parser.add_argument('-tol0', default=0.005, help='Initial tolerance for simulated annealing.')
parser.add_argument('-mintol', default=0.001, help='Minimum tolerance for simulated annealing.')
parser.add_argument('-decay', default=0.0025, help='Simulated annealing decay rate.')

args = parser.parse_args()

# TODO: allow user to import settings from JSON?
# initialize model parameters - note variable names must match R script
settings = {'c1': {'initial': 0.8,
                   'min': 0.1,
                   'max': 10.,
                   'sigma': 0.1,
                   'weight': 1.,
                   'prior': stats.lognorm(1)},
            'c2': {'initial': 1.0,
                   'min': 0.1,
                   'max': 10.,
                   'sigma': 0.1,
                   'weight': 1.,
                   'prior': stats.lognorm(0.5)},
            'p': {'initial': 0.5,
                  'min': 0.,
                  'max': 1.,
                  'sigma': 0.05,
                  'weight': 1.,
                  'prior': stats.uniform(loc=0, scale=1)},
            'rho': {'initial': 0.9,
                    'min': 0.,
                    'max': 1.,
                    'sigma': 0.05,
                    'weight': 1.,
                    'prior': stats.uniform(loc=0, scale=1)},
            'N': {'initial': 20000,
                  'min': 1000,
                  'max': 1e8,
                  'sigma': 2000,
                  'weight': 1.,
                  'prior': stats.norm(loc=8000, scale=3000)},
            't.end': {'initial': 30*52.,
                      'min': 520,
                      'max': 2600,
                      'sigma': 50,
                      'weight': 0.,  # fixed
                      'prior': stats.uniform(loc=520, scale=2600)}}


kam = Kamphir(settings, ncores=6, nthreads=6)
kam.set_target_tree(args.nwkfile)


# prevent previous log files from being overwritten
modifier = ''
tries = 0
while os.path.exists(args.logfile+modifier):
    tries += 1
    modifier = '.%d' % tries

logfile = open(args.logfile+modifier, 'w')
kam.abc_mcmc(logfile,
                skip=args.skip,
                tol0=args.tol0,
                mintol=args.mintol,
                decay=args.decay)
logfile.close()

