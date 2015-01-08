"""
Estimate epidemic model parameters by comparing simulations to "observed" phylogeny.
"""
import os
from pylab import *

from phyloK2 import *
import random

from copy import deepcopy
import time
from cStringIO import StringIO
import json
import subprocess
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
    
    def __init__(self, settings, script, driver,
                 ncores=1, nreps=10, nthreads=1, gibbs=False,
                 **kwargs):
        # call base class constructor
        PhyloKernel.__init__(self, **kwargs)

        self.settings = deepcopy(settings)
        self.target_tree = None

        self.current = {}
        self.proposed = {}
        for k, v in self.settings.iteritems():
            self.current.update({k: v['initial']})
            self.proposed.update({k: v['initial']})

        # locations of files
        self.pid = os.getpid()  # make paths unique to this process
        print 'initializing Kamphir with pid', self.pid

        self.path_to_tree = None
        self.path_to_input_csv = '/tmp/input_%d.csv' % self.pid
        self.path_to_label_csv = '/tmp/tips_%d.csv' % self.pid
        self.path_to_output_nwk = '/tmp/output_%d.nwk' % self.pid
        self.path_to_script = script
        self.driver = driver

        self.ntips = None
        self.tip_heights = []
        self.ref_denom = None  # kernel score of target tree to itself
        self.ncores = ncores  # number of processes for rcolgem simulation
        self.nreps = nreps
        self.nthreads = nthreads  # number of processes for PhyloKernel
        self.gibbs = gibbs

    def set_target_tree(self, path, delimiter=None, position=None):
        """
        Assign a Bio.Phylo Tree object to fit a model to.
        Parse tip dates from tree string in BEAST style.

        :param path: location of file containing Newick tree string
        :param delimiter: if set, partition tip label into tokens
        :param position: indicates which token denotes tip date
        :return: None
        """
        # TODO: If file contains more than one tree, then assign multiple trees.
        # TODO: Read states in from file.

        self.path_to_tree = path
        print 'reading in target tree from', path
        self.target_tree = Phylo.read(path, 'newick')

        tips = self.target_tree.get_terminals()
        self.ntips = len(tips)
        print 'read in', self.ntips, 'leaves'

        # constrain duration of simulation to tree depth
        if self.normalize != 'none':
            self.settings['t_end']['min'] = max(self.target_tree.depths().values())
            self.settings['t_end']['initial'] = self.settings['t_end']['min']
            self.settings['t_end']['max'] = 1.1*self.settings['t_end']['min']

            self.current['t_end'] = self.settings['t_end']['initial']
            self.proposed['t_end'] = self.settings['t_end']['initial']

        # parse tip heights from labels
        if delimiter is None:
            self.tip_heights = [''] * self.ntips
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
                    print 'Warning: Failed to parse tipdate from label', tip.name
                    tipdate = None  # gets interpreted as 0
                    pass

                tipdates.append(tipdate)

            self.tip_heights = [str(maxdate-t) if t else 0 for t in tipdates]

        # analyze target tree
        self.target_tree.ladderize()
        self.normalize_tree(self.target_tree, self.normalize)
        self.annotate_tree(self.target_tree)
        self.ref_denom = self.kernel(self.target_tree, self.target_tree)
    
    
    def proposal (self, tuning=1.0):
        """
        Generate a deep copy of parameters and modify one
        parameter value, given constraints (if any).
        :param tuning = factor to adjust sigma
        """
        if self.gibbs:
            # make deep copy
            for key in self.current.iterkeys():
                self.proposed[key] = self.current[key]

            # which parameter to adjust in proposal (component-wise)?
            choices = []
            for parameter in self.settings.iterkeys():
                choices.extend([parameter] * int(self.settings[parameter]['weight']))
            to_modify = random.sample(choices, 1)[0] # weighted sampling
            #to_modify = random.sample(self.proposed.keys(), 1)[0] # uniform sampling

            proposal_value = None
            current_value = self.proposed[to_modify]
            sigma = self.settings[to_modify]['sigma'] * tuning
            while True:
                proposal_value = current_value + random.normalvariate(0, sigma)
                if self.settings[to_modify].has_key('min') and proposal_value < self.settings[to_modify]['min']:
                    continue
                if self.settings[to_modify].has_key('max') and proposal_value > self.settings[to_modify]['max']:
                    continue
                break
            self.proposed[to_modify] = proposal_value
        else:
            # full-dimensional update
            for key in self.current.iterkeys():
                sigma = self.settings[key]['sigma'] * tuning
                while True:
                    if self.settings[key]['log'].upper()=='TRUE':
                        # log-normal proposal - NOTE mean and sigma are on natural log scale
                        proposal_value = random.lognormvariate(math.log(self.current[key]), sigma)
                    else:
                        # Gaussian
                        proposal_value = random.normalvariate(self.current[key], sigma)

                    if self.settings[key].has_key('min') and proposal_value < self.settings[key]['min']:
                        continue
                    if self.settings[key].has_key('max') and proposal_value > self.settings[key]['max']:
                        continue
                    break
                self.proposed[key] = proposal_value

    
    def prior (self, params):
        """
        Calculate the prior probability of a given parameter vector.
        """
        res = 1.
        for key in params.iterkeys():
            pass
            # work in progress

    def compute(self, tree, output=None):
        """
        Calculate kernel score.  Allow for MP execution.
        """
        try:
            tree.root.branch_length = 0.
            tree.ladderize()
            self.normalize_tree(tree, self.normalize)
            self.annotate_tree(tree)
        except:
            print 'ERROR: failed to prepare tree for kernel computation'
            print tree
            raise

        try:
            k = self.kernel(self.target_tree, tree)
            tree_denom = self.kernel(tree, tree)
        except:
            print 'ERROR: failed to compute kernel score for tree'
            print tree
            print self.target_tree
            raise

        #knorm = k / math.sqrt(self.ref_denom * tree_denom)
        try:
            knorm = math.exp(math.log(k) - 0.5*(math.log(self.ref_denom) + math.log(tree_denom)))
        except:
            print 'ERROR: failed to normalize kernel score'
            print k, self.ref_denom, tree_denom
            raise

        if output is None:
            return knorm

        output.put(knorm)  # MP


    def simulate(self):
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
        for item in self.proposed.iteritems():
            handle.write('%s,%f\n' % item)  # parameter name and value
        handle.close()

        # generate tip labels CSV file
        # TODO: take user-specified tip labels
        handle = open(self.path_to_label_csv, 'w')
        for i in range(self.ntips):
            handle.write('%d,%s\n' % (
                1 if i < (self.ntips*self.proposed['p']) else 2,
                self.tip_heights[i]
            ))
        handle.close()

        # external call to tree simulator script
        """
        p = subprocess.Popen([self.driver, self.path_to_script, self.path_to_input_csv,
                               self.path_to_label_csv, self.path_to_output_nwk],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        p.wait()
        """
        os.system(' '.join([self.driver, self.path_to_script, self.path_to_input_csv,
                            self.path_to_label_csv, self.path_to_output_nwk]))

        # retrieve trees from output file
        trees = []
        handle = open(self.path_to_output_nwk, 'rU')
        for line in handle:
            try:
                tree = Phylo.read(StringIO(line), 'newick')
            except:
                # NewickError: Number of open/close parentheses do not match
                print 'WARNING: Discarding mangled tree.'
                continue
            trees.append(tree)
        handle.close()

        #trees = Phylo.parse(self.path_to_output_nwk, 'newick')
        return trees

    def evaluate(self, trees=None):
        """
        Wrapper to calculate mean kernel score for a simulated set
        of trees given proposed model parameters.
        :param trees = list of Phylo Tree objects from simulations
                        in case we want to re-evaluate mean score (debugging)
        :return [mean] mean kernel score
                [trees] simulated trees (for debugging)
        """
        if trees is None:
            trees = self.simulate()

            if len(trees) < self.nreps:
                #print 'WARNING: tree sample size reduced to', len(trees)
                if len(trees) == 0:
                    print 'WARNING: none of the trees managed to coalesce in simulation time - returning 0.'
                    return 0.

        if self.nthreads > 1:
            # output = mp.Queue()
            #processes = [mp.Process(target=self.compute,
            #                        args=(trees[i], output)) for i in range(self.nthreads)]
            #map(lambda p: p.start(), processes)
            #map(lambda p: p.join(), processes)
            ## collect results and calculate mean
            #res = [output.get() for p in processes]
            try:
                async_results = [apply_async(pool, self.compute, args=(tree,)) for tree in trees]
            except:
                # TODO: dump trees to file for debugging
                raise

            #pool.close()  # prevent any more tasks from being added - once completed, workers exit
            map(mp.pool.ApplyResult.wait, async_results)
            results = [r.get() for r in async_results]

        else:
            # single-threaded mode
            results = [self.compute(tree) for tree in trees]

        try:
            mean = sum(results)/len(results)
        except:
            print res
            raise

        return mean


    def abc_mcmc(self, logfile, max_steps=1e5, tol0=0.01, mintol=0.0005, decay=0.0025, skip=1, first_step=0):
        """
        Use Approximate Bayesian Computation to sample from posterior
        density over model parameter space, given one or more observed
        trees.
        [sigma2] = variance parameter for Gaussian RBF
                   A higher value is more permissive.
        """
        # record settings in logfile header

        # report variables in alphabetical order
        keys = self.current.keys()
        keys.sort()

        logfile.write('# colgem_fitter.py log\n')
        logfile.write('# start time: %s\n' % time.ctime())
        logfile.write('# input file: %s\n' % self.path_to_tree)
        logfile.write('# annealing settings: tol0=%f, mintol=%f, decay=%f\n' % (tol0, mintol, decay))
        logfile.write('# MCMC settings: %s\n' % json.dumps(self.settings))
        logfile.write('# kernel settings: decay=%f normalize=%s tau=%f %s\n' % (
                      self.decayFactor, self.normalize, self.gaussFactor,
                      'gibbs' if self.gibbs else ''))
        
        cur_score = self.evaluate()
        step = first_step  # in case of restarting chain
        logfile.write('\t'.join(['state', 'score'] + keys))
        logfile.write('\n')

        # TODO: generalize screen and file log parameters
        while step < max_steps:
            self.proposal()  # update proposed values
            next_score = self.evaluate()
            if next_score > 1.0:
                print 'ERROR: next_score (%f) greater than 1.0, dumping proposal and EXIT' % next_score
                print self.proposal()
                sys.exit()
            
            # adjust tolerance, simulated annealing
            tol = (tol0 - mintol) * exp(-1. * decay * step) + mintol
            
            ratio = exp(-2.*(1.-next_score)/tol) / exp(-2.*(1.-cur_score)/tol)
            accept_prob = min(1., ratio)
            #step_down_prob = exp(-200.*(cur_score - next_score))
            #if next_score > cur_score or random.random() < step_down_prob:
            #rbf = math.exp(-(1-next_score)**2 / sigma2)  # Gaussian radial basis function
            
            # screen log
            #print '%d\t%1.5f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.5f\t%s' % (step, cur_score, next_score, 
            #    accept_prob, self.current['c1'], self.proposed['c1'], tol, time.ctime())
            # TODO: generalize log outputs for varying model parameters
            to_screen = '%d\t%1.5f\t%1.5f\t' % (step, cur_score, accept_prob)
            to_screen += '\t'.join(map(lambda x: str(round(x, 5)), [self.current[k] for k in keys]))
            print to_screen
            
            if random.random() < accept_prob:
                # accept proposal
                for key in self.current:
                    self.current[key] = self.proposed[key]
                cur_score = next_score
            
            if step % skip == 0:
                logfile.write('\t'.join(map(str, [step, cur_score] + [self.current[k] for k in keys])))
                logfile.write('\n')
            step += 1

if __name__ == '__main__':
    import os
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

    # first required arg is simulation Rscript
    parser.add_argument('script', help='Script for simulating trees.')
    parser.add_argument('settings',  help='JSON file containing model parameter settings.')
    parser.add_argument('-driver', default='Rscript', choices=['Rscript', 'python'],
                        help='Driver for executing script.')
    parser.add_argument('-nreps', default=10, type=int, help='Number of replicate trees to simulate.')

    # log settings
    parser.add_argument('nwkfile', help='File containing Newick tree string.')
    parser.add_argument('logfile', help='File to log ABC-MCMC traces.')
    parser.add_argument('-skip', default=1, help='Number of steps in ABC-MCMC to skip for log.')
    parser.add_argument('-overwrite', action='store_true', help='Allow overwrite of log file.')
    parser.add_argument('-restart', action='store_true', help='Restart chain from a log file.')

    # tree input settings
    parser.add_argument('-delimiter', default=None,
                        help='Delimiter used in tip label to separate fields.')
    parser.add_argument('-datefield', default=-1,
                        help='Index (from 0) of field in tip label containing date.')
    # annealing settings
    parser.add_argument('-tol0', type=float, default=0.01,
                        help='Initial tolerance for simulated annealing.')
    parser.add_argument('-mintol', type=float, default=0.005,
                        help='Minimum tolerance for simulated annealing.')
    parser.add_argument('-toldecay', type=float, default=0.0025,
                        help='Simulated annealing decay rate.')

    # MCMC settings
    parser.add_argument('-maxsteps', type=int, default=1e5,
                        help='Maximum number of steps to run chain sample.')
    parser.add_argument('-gibbs', action='store_true',
                        help='Perform component-wise update; otherwise full-dimensional '
                             'Metropolis is the default.')

    # kernel settings
    parser.add_argument('-kdecay', default=0.35, type=float,
                        help='Decay factor for tree shape kernel. Lower values penalize large subset '
                             'trees more severely.')
    parser.add_argument('-tau', default=0.5, type=float,
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

    args = parser.parse_args()

    # initialize multiprocessing thread pool at global scope
    pool = mp.Pool(processes=args.nthreads)

    # start analysis
    if args.restart:
        # TODO: work in progress
        logfile = open(args.logfile, 'rU')
        for line in logfile:
            if line.startswith('#'):
                if line.startswith('# MCMC settings:'):
                    settings = json.loads(line.split('settings: ')[-1])
                if line.startswith('# annealing settings:'):
                    tol0, mintol, decay = map(float, line.strip('\n').split('settings: ')[-1].split(', '))
                if line.startswith('# kernel settings:'):
                    kdecay = float(line.strip('\n').split('=')[-1])
            else:
                if line.endswith('\n'):
                    # complete row
                    items = line.strip('\n').split()

        logfile.close()

    else:
        # initialize model parameters - note variable names must match R script
        handle = open(args.settings, 'rU')
        settings = json.loads(handle.read())
        handle.close()

        kam = Kamphir(settings=settings,
                      driver=args.driver,
                      script=args.script,
                      ncores=args.ncores,
                      nthreads=args.nthreads,
                      decayFactor=args.kdecay,
                      normalize=args.normalize,
                      gaussFactor=args.tau,
                      gibbs=args.gibbs,
                      nreps=args.nreps)

        kam.set_target_tree(args.nwkfile, delimiter=args.delimiter, position=args.datefield)

        # prevent previous log files from being overwritten
        modifier = ''
        tries = 0
        while os.path.exists(args.logfile+modifier) and not args.overwrite:
            tries += 1
            modifier = '.%d' % tries

        logfile = open(args.logfile+modifier, 'w')
        kam.abc_mcmc(logfile,
                        skip=args.skip,
                        tol0=args.tol0,
                        mintol=args.mintol,
                        decay=args.toldecay)
        logfile.close()

