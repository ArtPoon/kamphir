"""
Driver script for fitting a differential risk SI model to tree shape.
"""
import os
import argparse
import json
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
parser.add_argument('-toldecay', default=0.0025, help='Simulated annealing decay rate.')
parser.add_argument('-kdecay', default=0.2, help='Decay factor for tree shape kernel.  '
                                                 'Lower values penalize large subset trees more severely.')

parser.add_argument('-ncores', default=6, help='Number of processes for tree simulation (rcolgem).')
parser.add_argument('-nthreads', default=6, help='Number of processes for kernel computation.')


args = parser.parse_args()

# TODO: allow user to import settings from JSON?
# initialize model parameters - note variable names must match R script
handle = open('tests/settings.json', 'rU')
settings = json.loads(handle.read())
handle.close()

kam = Kamphir(settings, ncores=6, nthreads=6, decayFactor=args.kdecay, normalize='none')
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
                decay=args.toldecay)
logfile.close()

