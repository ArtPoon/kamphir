#!/usr/bin/env python
"""
Estimate epidemic model parameters by comparing simulations to "observed" phylogeny.
"""
import sys
import os
from phyloK2 import *
import random

from copy import deepcopy
import time
from cStringIO import StringIO
import math
from scipy import stats

import logging

logging.getLogger().setLevel(logging.DEBUG)
FNULL = open(os.devnull, 'w')

# see http://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error/24673524#24673524
import dill

def run_dill_encoded(what):
    fun, args = dill.loads(what)
    return fun(*args)

def apply_async(pool, fun, args):
    return pool.apply_async(run_dill_encoded, (dill.dumps((fun, args)),))

class Kamphir (PhyloKernel):
    """
    Derived class of PhyloKernel for estimating epidemic model parameters
    by simulating coalescent trees and comparing simulations to the reference
    tree by the tree shape kernel function.
    
    [target_tree] = Bio.Phylo tree object to fit model to.
    
    [params] = dictionary of parameter values, key = parameter name as
        recognized by colgem2 HIVmodel, value = dictionary with following:
        'value': parameter value,
        'sigma': sigma value of Gaussian proposal,
        'min': minimum parameter value (optional),
        'max': maximum parameter value (optional)
    """
    
    def __init__(self, settings, script, driver, simfunc,
                 ncores=1, nreps=10, nthreads=1, gibbs=False, use_priors=False,
                 **kwargs):
        # call base class constructor
        PhyloKernel.__init__(self, **kwargs)

        self.use_priors = use_priors
        self.settings = deepcopy(settings)
        self.target_trees = []

        self.current = {}
        self.proposed = {}
        self.priors = {}
        for k, v in self.settings.iteritems():
            self.current.update({k: v['initial']})
            self.proposed.update({k: v['initial']})
            if self.use_priors:
                frozen_dist = eval('stats.'+v['prior'])
                self.priors.update({k: frozen_dist})

        # locations of files
        self.pid = os.getpid()  # make paths unique to this process
        logging.info('initializing Kamphir with pid {}'.format(self.pid))

        self.path_to_tree = None
        self.path_to_input_csv = '/tmp/input_%d.csv' % self.pid
        self.path_to_label_csv = '/tmp/tips_%d.csv' % self.pid
        self.path_to_output_nwk = '/tmp/output_%d.nwk' % self.pid
        self.path_to_script = script
        self.driver = driver

        # rcolgem functions
        self.simfunc = simfunc

        self.ntips = []
        self.tip_heights = []  # list of lists
        self.tree_heights = []
        self.ref_denom = []  # kernel score of target tree to itself

        self.ncores = ncores  # number of processes for rcolgem simulation
        self.nreps = nreps
        self.nthreads = nthreads  # number of processes for PhyloKernel
        self.gibbs = gibbs

        if "odelog" in kwargs:
            self.odelog = kwargs["odelog"]
        else:
            self.odelog = None

    def set_target_trees(self, path, treenum, delimiter=None, position=None):
        """
        Assign a Bio.Phylo Tree object to fit a model to.
        Parse tip dates from tree string in BEAST style.

        :param path: location of file containing Newick tree string
        :param delimiter: if set, partition tip label into tokens
        :param position: indicates which token denotes tip date
        :return: None
        """
        # TODO: Read states in from file.

        self.path_to_tree = path

        # reset lists
        self.target_trees = []  # tuple (newick string, tree height, [tip heights], denom)

        for index, tree in enumerate(Phylo.parse(path, 'newick')):
            if treenum is not None and index != treenum:
                # user asked to process only one tree from this file
                continue

            # record this before normalizing
            tree_height = max(tree.depths().values())

            tree.ladderize()
            self.normalize_tree(tree, self.normalize)
            self.annotate_tree(tree)

            ref_denom = self.kernel(tree, tree)

            tips = tree.get_terminals()
            ntips = len(tips)

            # record tip heights (list of lists)
            if delimiter is None:
                tip_heights = [0.] * ntips
            else:
                maxdate = 0
                tipdates = []
                for tip in tips:
                    try:
                        items = tip.name.strip("'").split(delimiter)
                        tipdate = float(items[position])
                        if tipdate > maxdate:
                            maxdate = tipdate
                    except:
                        logging.warning('Failed to parse tipdate from label {}'.format(tip.name))
                        tipdate = None  # gets interpreted as 0
                        pass

                    tipdates.append(tipdate)

                tip_heights = [str(maxdate-t) if t else 0 for t in tipdates]

            self.target_trees.append((tree, tree_height, tip_heights, ref_denom))

        if len(self.target_trees) == 0:
            # we didn't read any of the trees from the file!
            logging.error('File did not contain any Newick tree strings, '
                          'or -treenum ({}) exceeds number of trees!'.format(treenum))
            sys.exit()


    def proposal (self, tuning=1.0, max_attempts=100):
        """
        Generate a deep copy of parameters and modify one
        parameter value, given constraints (if any).
        :param tuning = factor to adjust sigma
        """

        # make deep copy
        for key in self.current.iterkeys():
            self.proposed[key] = self.current[key]
        
        if self.gibbs:
            # which parameter to adjust in proposal (component-wise)?
            choices = []
            for parameter in self.settings.iterkeys():
                choices.extend([parameter] * int(self.settings[parameter]['weight']))
            to_modify = random.sample(choices, 1)[0] # weighted sampling
            #to_modify = random.sample(self.proposed.keys(), 1)[0] # uniform sampling
        else:
            # full dimensional update
            to_modify = self.settings.keys()

        for key in to_modify:
            sigma = self.settings[key]['sigma'] * tuning
            if sigma == 0:
                # no modification
                continue

            attempts = 0
            this_min = self.settings[key].get('min', None)
            this_max = self.settings[key].get('max', None)
            while True:
                attempts += 1
                if attempts > max_attempts:
                    logging.error('Failed to update proposal, check initial/min/max settings.')
                    sys.exit()
                if self.settings[key]['log'].upper()=='TRUE':
                    # log-normal proposal - NOTE mean and sigma are on natural log scale
                    proposal_value = random.lognormvariate(math.log(self.current[key]), sigma)
                else:
                    # Gaussian
                    proposal_value = random.normalvariate(self.current[key], sigma)

                if this_min is not None and proposal_value < this_min:
                    delta = this_min - proposal_value  # how far past the minimum are we?
                    proposal_value = this_min + delta  # reflect this amount up from minimum

                if this_max is not None and proposal_value > this_max:
                    delta = proposal_value - this_max
                    proposal_value = this_max - delta

                # one more time to check that we are within bounds
                if this_min is not None and proposal_value < this_min:
                    continue  # try again

                self.proposed[key] = proposal_value
                break

    
    def log_priors (self):
        """
        Calculate the natural log-transformed prior probabilities for current and proposed
        parameter values.
        """
        retval = {'proposal': 0., 'current': 0.}
        """
        for key in self.current.iterkeys():
            retval['proposal'] += math.log(self.priors[key].pdf(self.current[key]))
            retval['current'] += math.log(self.priors[key].pdf(self.proposed[key]))
        """
        return retval

    def compute(self, tree, target_tree, ref_denom, output=None):
        """
        Calculate kernel score.  Allow for MP execution.
        """
        try:
            tree.root.branch_length = 0.
            tree.ladderize()
            self.normalize_tree(tree, self.normalize)
            self.annotate_tree(tree)
        except:
            logging.error('Failed to prepare tree for kernel computation')
            logging.error(tree)
            raise

        try:
            k = self.kernel(target_tree, tree)
            tree_denom = self.kernel(tree, tree)
        except:
            logging.error('Failed to compute kernel score for tree')
            logging.error(tree)
            logging.error(target_tree)
            raise

        #knorm = k / math.sqrt(self.ref_denom * tree_denom)
        try:
            knorm = math.exp(math.log(k) - 0.5*(math.log(ref_denom) + math.log(tree_denom)))
        except:
            logging.error('Failed to normalize kernel score')
            logging.error("{} {} {}".format(k, ref_denom, tree_denom))
            raise

        if output is None:
            return knorm

        output.put(knorm)  # MP

    def simulate_internal(self, tree_height, tip_heights):
        """
        Simulate trees using class function simfunc.
        Convert resulting Newick tree strings into Phylo objects.
        :return: List of Phylo BaseTree objects.
        """
        result = self.simfunc(self.proposed, tree_height, tip_heights)
        if self.odelog is None and result:
            newicks, ode_solution = result
        else:
            newicks = result

        trees = []
        for newick in newicks:
            try:
                tree = Phylo.read(StringIO(newick), 'newick')
            except:
                continue
            trees.append(tree)

        if self.odelog is not None and result:
            if not os.path.exists(self.odelog):
                ode_solution.to_csvfile(self.odelog, col_names=False, row_names=False, sep='\t')
            else:
                ode_solution.to_csvfile(self.odelog, col_names=False, row_names=False, sep='\t', append=True)

        return trees

    def simulate_external(self, tree_height, tip_heights, prune=True):
        """
        Estimate the mean kernel distance between the reference tree and
        trees simulated under the given model parameters.
        :returns List of Phylo Tree objects
        """
        # TODO: allow user to set arbitrary driver Rscript
        # TODO: generalize tip label annotation

        # generate input control CSV file
        handle = open(self.path_to_input_csv, 'w')
        handle.write('n.cores,%d\n' % self.ncores)  # parallel or serial execution
        handle.write('nreps,%d\n' % self.nreps)  # number of replicates
        handle.write('t_end,%f\n' % tree_height)
        for item in self.proposed.iteritems():
            handle.write('%s,%f\n' % item)  # parameter name and value
        handle.close()

        # generate tip labels CSV file
        # TODO: take user-specified tip labels
        handle = open(self.path_to_label_csv, 'w')
        for tip_height in tip_heights:
            handle.write('%d,%s\n' % (
                1 + int(random.random() < self.proposed['p']),
                #1 if i < (self.ntips*self.proposed['p']) else 2,
                tip_height
            ))
        handle.close()

        # external call to tree simulator script
        os.system(' '.join([self.driver, self.path_to_script, self.path_to_input_csv,
                            self.path_to_label_csv, self.path_to_output_nwk]) + ' >/dev/null')

        # retrieve trees from output file
        trees = []
        try:
            handle = open(self.path_to_output_nwk, 'rU')
        except IOError:
            # file does not exist, simulation failed
            return []
        
        for line in handle:
            try:
                tree = Phylo.read(StringIO(line), 'newick')
            except:
                # NewickError: Number of open/close parentheses do not match
                #print 'WARNING: Discarding mangled tree.'
                continue

            trees.append(tree)
        handle.close()

        return trees

    def prune_tree(self, tree):
        """
        Sample a random number of tips in the tree and prune the rest.
        Only used for forward-time simulation.
        :param tree:
        :param target_size:
        :return:
        """
        tips = tree.get_terminals()
        try:
            tips2 = random.sample(tips, self.ntips)
        except ValueError:
            tips2 = tips

        for tip in tips:
            if tip in tips2:
                continue
            _ = tree.prune(tip)

        return tree


    def evaluate(self):
        """
        Wrapper to calculate mean kernel score for a simulated set
        of trees given proposed model parameters.
        :param trees = list of Phylo Tree objects from simulations
                        in case we want to re-evaluate mean score (debugging)
        :return [mean] mean kernel score
                [trees] simulated trees (for debugging)
        """
        retval = 0.
        total_ntips = 0

        # iterate over target trees
        for target_tree, tree_height, tip_heights, ref_denom in self.target_trees:
            ntips = len(tip_heights)
            total_ntips += ntips

            # simulate trees for this target tree
            logging.debug("Simulating trees")
            if self.simfunc is None:
                trees = self.simulate_external(tree_height, tip_heights)
            else:
                trees = self.simulate_internal(tree_height, tip_heights)
            logging.debug("Finished simulating trees")

            if len(trees) == 0:
                # failed simulation
                return None

            if self.nthreads > 1:
                try:
                    async_results = [apply_async(pool,
                                                 self.compute,
                                                 args=(tree, target_tree, ref_denom))
                                     for tree in trees]
                except:
                    # TODO: dump trees to file for debugging
                    raise

                map(mp.pool.ApplyResult.wait, async_results)
                results = [r.get() for r in async_results]
            else:
                # single-threaded mode
                results = [self.compute(tree, target_tree, ref_denom) for tree in trees]

            # sum weighted by size of tree
            retval += sum(results)/len(results) * ntips

        return retval / total_ntips


    def abc_mcmc(self, logfile, max_steps=1e5, tol0=0.01, mintol=0.0005, decay=0.0025, skip=1, first_step=0):
        """
        Use Approximate Bayesian Computation to sample from posterior
        density over model parameter space, given one or more observed
        trees.
        [sigma2] = variance parameter for Gaussian RBF
                   A higher value is more permissive.
        :param odelog = log file for output of numerical ODE solution
        """
        # record settings in logfile header

        # report variables in alphabetical order
        keys = self.current.keys()
        keys.sort()

        logfile.write('# Kamphir log\n')
        logfile.write('# start time: %s\n' % time.ctime())
        logfile.write('# input file: %s\n' % self.path_to_tree)
        logfile.write('# annealing settings: tol0=%f, mintol=%f, decay=%f\n' % (tol0, mintol, decay))
        logfile.write('# MCMC settings: %s\n' % json.dumps(self.settings))
        logfile.write('# kernel settings: decay=%f normalize=%s tau=%f %s\n' % (
                      self.decayFactor, self.normalize, self.gaussFactor,
                      'gibbs' if self.gibbs else ''))

        logging.info('Calculating initial kernel score')
        cur_score = self.evaluate()
        if cur_score is None:
            logging.error('Failed to simulate trees under initial parameter values')
            sys.exit()
        logging.info('Current score is {}'.format(cur_score))

        step = first_step  # in case of restarting chain
        logfile.write('\t'.join(['state', 'score', 'prior'] + keys))
        logfile.write('\n')
        logfile.flush()

        # TODO: generalize screen and file log parameters
        print('\t'.join(['state', 'score', 'prior'] + keys)) # header for screen log
        while step < max_steps:
            next_score = None
            while next_score is None:
                self.proposal()  # update proposed values
                next_score = self.evaluate()  # returns None if simulations fail
                
            if next_score > 1.0 or next_score < 0.0:
                logging.error('next_score ({}) outside interval [0,1], dumping proposal and EXIT'.format(next_score))
                logging.error(self.proposal())
                sys.exit()
            
            # adjust tolerance, simulated annealing
            tol = (tol0 - mintol) * math.exp(-1. * decay * step) + mintol
            log_prior = self.log_priors()

            ratio = math.exp(-2.*(1.-next_score)/tol) / math.exp(-2.*(1.-cur_score)/tol)
            if self.use_priors:
                ratio *= math.exp(log_prior['proposal'] - log_prior['current'])
            accept_prob = min(1., ratio)

            # screen log
            to_screen = '%d\t%1.5f\t%1.5f\t%1.5f\t' % (step, cur_score, log_prior['proposal'], accept_prob)
            to_screen += '\t'.join(map(lambda x: str(round(x, 5)), [self.current[k] for k in keys]))
            print to_screen
            
            if random.random() < accept_prob:
                # accept proposal
                for key in self.current:
                    self.current[key] = self.proposed[key]
                cur_score = next_score
            
            if step % skip == 0:
                logfile.write('\t'.join(map(str, [step, cur_score, log_prior['proposal']] + [self.current[k] for k in keys])))
                logfile.write('\n')
                logfile.flush()
            step += 1

if __name__ == '__main__':
    import argparse
    import json
    from multiprocessing import cpu_count

    parser = argparse.ArgumentParser(description='KAMPHIR - Kernel-assisted ABC-MCMC for '
                                                 'Phylodynamic Inference\n==============='
                                                 '======================================='
                                                 '=======',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog='KAMPHIR uses Approximate Bayesian Computation to fit any model that '
                                            'can be used to generate a tree.')

    # positional arguments (required)
    parser.add_argument('model', help='Model to simulate trees with Rcolgem.  Use "*" to fit '
                                      'a model using another program and driver script.',
                        choices=['*', 'SI', 'SI2', 'DiffRisk', 'Stages'])
    parser.add_argument('settings', help='JSON file containing model parameter settings.  Ignored if'
                                         'restarting from log file (-restart).')
    parser.add_argument('nwkfile', help='File containing Newick tree string.')
    parser.add_argument('logfile', help='File to log ABC-MCMC traces.')

    # non-Rcolgem methods
    parser.add_argument('-script', default=None,
                        help='Driver script implementing model.  See examples in /drivers folder.')
    parser.add_argument('-driver', choices=['Rscript', 'python'],
                        help='Driver for executing script.')

    # log settings
    parser.add_argument('-skip', default=1, help='Number of steps in ABC-MCMC to skip for log.')
    parser.add_argument('-overwrite', action='store_true', help='Allow overwrite of log file.')
    parser.add_argument('-restart', default=None, help='Restart chain from log file specified.')
    parser.add_argument('-odelog', default=None, help='Log file for numerical ODE solutions.')

    # tree input settings
    parser.add_argument('-delimiter', default=None,
                        help='Delimiter used in tip label to separate fields.')
    parser.add_argument('-datefield', default=-1,
                        help='Index (from 0) of field in tip label containing date.')
    parser.add_argument('-treenum', type=int, default=None,
                        help='Index of tree in file to process.  Defaults to all.')
    
    # annealing settings
    parser.add_argument('-tol0', type=float, default=0.01,
                        help='Initial tolerance for simulated annealing.')
    parser.add_argument('-mintol', type=float, default=0.0025,
                        help='Minimum tolerance for simulated annealing.')
    parser.add_argument('-toldecay', type=float, default=0.0025,
                        help='Simulated annealing decay rate.')

    # MCMC settings
    parser.add_argument('-nreps', default=10, type=int, help='Number of replicate trees to simulate.')
    parser.add_argument('-maxsteps', type=int, default=1e5,
                        help='Maximum number of steps to run chain sample.')
    parser.add_argument('-gibbs', action='store_true',
                        help='Perform component-wise update; otherwise full-dimensional '
                             'Metropolis is the default.')
    parser.add_argument('-prior', action='store_true', help='Use prior distributions.')

    # kernel settings
    parser.add_argument('-kdecay', default=0.2, type=float,
                        help='Decay factor for tree shape kernel. Lower values penalize large subset '
                             'trees more severely.')
    parser.add_argument('-tau', default=2.0, type=float,
                        help='Precision for Gaussian radial basis function penalizing branch length '
                             'discordance. Lower values penalize more severely.  CAVEAT: if normalize '
                             'is set to "none" then make sure this parameter is scaled to the '
                             'typical branch length of the target tree.')
    parser.add_argument('-normalize', default='mean', choices=['none', 'mean', 'median'],
                        help='Scale branch lengths so trees of different lengths can be compared.')

    # parallelization
    parser.add_argument('-ncores', type=int, default=cpu_count(),
                        help='Number of processes for tree simulation (rcolgem).')
    parser.add_argument('-nthreads', type=int, default=cpu_count(),
                        help='Number of processes for kernel computation.')

    # reproducibility
    parser.add_argument('-seed', type=int, default=None,
                        help='Random seed, to make runs reproducible')

    args = parser.parse_args()

    # set the random seed
    random.seed(args.seed)

    # initialize multiprocessing thread pool at global scope
    pool = mp.Pool(processes=args.nthreads)

    # recover from log file if requested
    if args.restart:
        logfile = open(args.restart, 'rU')
        header = None
        tol0 = args.tol0
        mintol = args.mintol
        decay = args.toldecay
        items = []  # will hold the last

        for line in logfile:
            if line.startswith('#'):
                if line.startswith('# MCMC settings:'):
                    settings = json.loads(line.split('settings: ')[-1])
                if line.startswith('# annealing settings:'):
                    tol0, mintol, decay = map(float, [x.split('=')[-1]
                                                      for x in line.strip('\n').split('settings: ')[-1].split(', ')])
                if line.startswith('# kernel settings:'):
                    items = line.strip('\n').split('settings: ')[-1].split()
                    for item in items:
                        key, value = item.split('=')
                        if key == 'decay':
                            args.kdecay = float(value)
                        elif key == 'normalize':
                            args.normalize = value
                        elif key == 'tau':
                            args.tau = float(value)
                        else:
                            logging.error('Unrecognized key {} when parsing log file for restart'.format(key))
                            sys.exit()
            else:
                if line.endswith('\n'):
                    # complete row
                    items = line.strip('\n').split()
                    if header is None:
                    # take the first non-commented line as header row
                        header = items

        logfile.close()

        # reset initial values in settings JSON
        state = 0
        for i, key in enumerate(header):
            value = items[i]
            if key in settings:
                settings[key]['initial'] = float(value)
            if key == 'state':
                state = int(value)

        # adjust tolerance parameters for state
        args.mintol = mintol
        args.tol0 = (tol0 - mintol) * math.exp(-1. * decay * state) + mintol
        args.toldecay = decay


    else:
        if args.settings is None:
            logging.error('Settings is required if not restarting from log file')
            sys.exit()

        # initialize model parameters - note variable names must match R script
        handle = open(args.settings, 'rU')
        settings = json.loads(handle.read())
        handle.close()

    # select model
    simfunc = None
    if args.model == '*':
        if args.script is None:
            logging.error('Must specify (-script) if (-model) is "*".')
            sys.exit()
        # simfunc remains set to None
    else:
        import rcolgem
        r = rcolgem.Rcolgem(ncores=args.ncores, nreps=args.nreps)
        if args.model == 'SI':
            r.init_SI_model()
            simfunc = r.simulate_SI_trees
        elif args.model == 'SI2':
            r.init_SI_model()
            simfunc = r.simulate_SI2_trees
        elif args.model == 'DiffRisk':
            r.init_DiffRisk_model()
            simfunc = r.simulate_DiffRisk_trees
        elif args.model == 'Stages':
            r.init_stages_model()
            simfunc = r.simulate_stages_trees
        else:
            logging.error('Unrecognized rcolgem model type {}\n'
                          'Currently only SI, SI2, DiffRisk, and Stages are supported'
                          .format(args.model))
            sys.exit()

    kam = Kamphir(settings=settings,
                  driver=args.driver,
                  simfunc=simfunc,
                  script=args.script,
                  ncores=args.ncores,
                  nthreads=args.nthreads,
                  decayFactor=args.kdecay,
                  normalize=args.normalize,
                  gaussFactor=args.tau,
                  gibbs=args.gibbs,
                  nreps=args.nreps,
                  use_priors=args.prior,
                  odelog=args.odelog)

    kam.set_target_trees(args.nwkfile, delimiter=args.delimiter, position=args.datefield,
                        treenum=args.treenum)

    # prevent previous log files from being overwritten
    modifier = ''
    tries = 0
    while os.path.exists(args.logfile+modifier) and not args.overwrite:
        tries += 1
        modifier = '.%d' % tries

    logfile = open(args.logfile+modifier, 'w')
    kam.abc_mcmc(logfile,
                    max_steps=args.maxsteps,
                    skip=args.skip,
                    tol0=args.tol0,
                    mintol=args.mintol,
                    decay=args.toldecay)
    logfile.close()
