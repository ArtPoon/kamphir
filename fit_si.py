"""
Driver script for fitting a simple SI model to tree shape.
"""
import os
import argparse
from kamphir import Kamphir

parser = argparse.ArgumentParser(description='Fit a simple SI model to '
                                             'the shape of a phylogenetic tree using '
                                             'kernel-assisted ABC-MCMC',
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
handle = open('tests/settings_si.json', 'rU')
settings = json.loads(handle.read())
handle.close()

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

