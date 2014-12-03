"""
Estimate epidemic model parameters by comparing simulations to "observed" phylogeny.
"""

import sys
sys.path.reverse() # to resolve some path problem

import os

from pylab import *

import colgem3
from Bio import Phylo

from phyloK2 import *
import random

from copy import deepcopy
import time

import scipy.stats as stats


class Fitter (PhyloKernel):
    """
    Derived class of PhyloKernel for estimating epidemic model parameters
    by simulating coalescent trees and comparing simulations to the reference
    tree by the tree shape kernel function.
    
    [ref_tree] = Bio.Phylo tree object to fit model to.
    
    [params] = dictionary of parameter values, key = parameter name as
        recognized by colgem2 HIVmodel, value = dictionary with following:
        'value': parameter value,
        'sigma': sigma value of Gaussian proposal,
        'min': minimum parameter value (optional),
        'max': maximum parameter value (optional)
    """
    
    def __init__(self, settings, initial):
        # call base class constructor
        PhyloKernel.__init__(self)
        self.initial = deepcopy(initial)
        self.settings = deepcopy(settings)
        self.ref_tree = None
        self.current = dict([(k,v) for k,v in initial.iteritems()])
        self.proposed = dict([(k,v) for k,v in initial.iteritems()])
        self.path_to_tree = None

        self.path_to_input_csv = '/tmp/input.csv'
        self.path_to_label_csv = '/tmp/tips.csv'
        self.path_to_output_nwk = '/tmp/output.nwk'
    
    def set_target_tree(self, path):
        """
        Assign a Bio.Phylo Tree object to fit a model to.
        TODO: If file contains more than one tree, then assign
        multiple trees.
        """
        self.path_to_tree = path
        print 'reading in target tree from', path
        self.ref_tree = Phylo.read(path, 'newick')
        print 'read in', len(self.ref_tree.get_terminals()), 'leaves'
        self.ref_tree.ladderize()
        self.normalize_tree(self.ref_tree, 'mean')
        self.annotate_tree(self.ref_tree)
        self.ref_denom = self.kernel(self.ref_tree, self.ref_tree)
    
    
    def proposal (self):
        """
        Generate a deep copy of parameters and modify one
        parameter value, given constraints (if any).
        """
        for key in self.current.iterkeys():
            self.proposed[key] = self.current[key]
        
        # which parameter to adjust in proposal
        choices = []
        for parameter in self.settings.iterkeys():
            choices.extend([parameter] * int(self.settings[parameter]['weight']))
        to_modify = random.sample(choices, 1)[0] # weighted sampling
        #to_modify = random.sample(self.proposed.keys(), 1)[0] # uniform sampling
        
        current_value = self.proposed[to_modify]
        sigma = self.settings[to_modify]['sigma']
        while True:
            proposal_value = current_value + random.normalvariate(0, sigma)
            if self.settings[to_modify].has_key('min') and proposal_value < self.settings[to_modify]['min']:
                continue
            if self.settings[to_modify].has_key('max') and proposal_value > self.settings[to_modify]['max']:
                continue
            break
        self.proposed[to_modify] = proposal_value
    
    def prior (self, params):
        """
        Calculate the prior probability of a given parameter vector.
        """
        res = 1.
        for key in params.iterkeys():
            pass
            # work in progress
        
    def evaluate (self, ntrees=5):
        """
        Estimate the mean kernel distance between the reference tree and
        trees simulated under the given model parameters.
        TODO: make N a free parameter
        """
        
        # simulated trees should have same number of tips as reference
        n = len(self.ref_tree.get_terminals())

        # generate input control CSV file
        handle = open(self.path_to_input_csv, 'w')

        handle.close()

        # external call to Rscript
        os.system('Rscript simulate.DiffRisk.R %s %s %s' % (input_csv, labels_csv, output_nwk))

        #solve ODE system
        m = DifferentialRiskModel(c1 = self.proposed['c1'], 
                                  c2 = self.proposed['c2'],
                                  rho = self.proposed['rho'],
                                  #x0 = array([4999., 5000., 1., 0]))
                                  x0 = array([self.proposed['N']*self.proposed['p'] - 1, 
                                              self.proposed['N']*(1-self.proposed['p']), 
                                              1., 
                                              0]))
        
        sf = float(n) / m.N # sampling fraction
        try:
            n1 = int(round(sf * m.i1[-1]))
        except:
            print sf
            print m.i1
            raise
        
        n2 = n - n1
        # generate taxon labels
        taxa = []
        for i in range(n):
            taxa.append('_%i_%d' % (i, 0 if i < n1 else 1))
    
        sampleStates = dict(zip(taxa, [eye(2)[0]]*n1 + [eye(2)[1]]*n2) ) # observed i1 and i2
        sampleTimes = dict(zip( taxa, [m.t[-1]]*n)) # homochronous sample at last simulation time
    
        F, G, Y = m.FGY() #births, migration, and population size
        fgy = colgem3.FGY(m.t, F, G, Y)
        
        # simulate tree
        #trees, nwks, A, daf = colgem3.simulate_coalescent(sampleTimes, sampleStates, fgy, singleMRCA = True)
        try:
            trees = colgem3.simulate_coalescent(sampleTimes, sampleStates, fgy, singleMRCA = True, ntrees=ntrees)[0]
        except:
            return 0.  # force reject
        
        res = []
        for tree in trees:
            tree.root.branch_length = 0.
            tree.ladderize()
            self.normalize_tree(tree, 'mean')
            self.annotate_tree(tree)
            k = self.kernel(self.ref_tree, tree)
            #print k/self.ref_denom  # examine variation among replicates
            tree_denom = self.kernel(tree, tree)  # normalize
            res.append(k / math.sqrt(self.ref_denom * tree_denom))
            #res.append(k/self.ref_denom)
            
            #print k / math.sqrt(self.ref_denom * tree_denom)
        
        
        # return expected value across replicates
        return sum(res)/len(res)
    
    
    def abc_mcmc (self, logfile, max_steps = 1e5, tol0=0.01, mintol=0.0005, decay=0.0025, skip=1):
        """
        Use Approximate Bayesian Computation to sample from posterior
        density over model parameter space, given one or more observed
        trees.
        [sigma2] = variance parameter for Gaussian RBF
                   A higher value is more permissive.
        """
        # record settings in logfile header
        logfile.write('# colgem_fitter.py log\n')
        logfile.write('# start time: %s\n' % time.ctime())
        logfile.write('# input file: %s\n' % self.path_to_tree)
        logfile.write('# annealing settings: tol0=%f, mintol=%f, decay=%f\n' % (tol0, mintol, decay))
        
        logfile.write('# proposal settings: ')
        logfile.write('\t'.join(['%s:(%f,%f); s=%f; w=%1.2f' % (k, v['min'], v['max'], v['sigma'], v['weight']) 
                                for k, v in self.settings.iteritems()]))
        logfile.write('\n')
        
        logfile.write('# initial: ')
        logfile.write(' '.join(['%s=%f' % (k, v) for k, v in self.initial.iteritems()]))
        logfile.write('\n')
        
        cur_score = self.evaluate()
        step = 0
        logfile.write('\t'.join(['state', 'score', 'c1', 'c2', 'p', 'rho', 'N', 'T']))
        logfile.write('\n')
        
        while step < max_steps:
            self.proposal()  # update proposed values
            next_score = self.evaluate(ntrees=100)
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
            to_screen = '%d\t%1.5f\t%1.5f\t' % (step, cur_score, accept_prob)
            to_screen += '\t'.join(map(lambda x: str(round(x, 5)), [
                self.current['c1'], 
                self.current['c2'], 
                self.current['p'], 
                self.current['rho'], 
                self.current['N'],
                self.current['T']]))
            print to_screen
            
            if random.random() < accept_prob:
                # accept proposal
                for key in self.current:
                    self.current[key] = self.proposed[key]
                cur_score = next_score
            
            if step % skip == 0:
                logfile.write('\t'.join(map(str, [step, cur_score,
                                              self.current['c1'], 
                                              self.current['c2'],
                                              self.current['p'],
                                              self.current['rho'],
                                              self.current['N'],
                                              self.current['T']])))
                logfile.write('\n')
            step += 1


# initialize model parameters - note variable names must match R script
# TODO: add prior distribution functions from scipy.stats
settings = {'c1': {'min': 0.1, 'max': 10., 'sigma': 0.05, 'weight': 1., 'prior': stats.lognorm(1)},
          'c2': {'min': 0.1, 'max': 10., 'sigma': 0.25, 'weight': 0., 'prior': stats.lognorm(0.5)},
          'p': {'min': 0., 'max': 1., 'sigma': 0.05, 'weight': 0., 'prior': stats.uniform(loc=0, scale=1)},
          'rho': {'min': 0., 'max': 1.1, 'sigma': 0.05, 'weight': 0., 'prior': stats.uniform(loc=0, scale=1)},
          'N': {'min': 1000, 'max': 1e8, 'sigma': 5000, 'weight': 2., 'prior': stats.norm(loc=8000, scale=3000)},
          't.end': {'min': 520, 'max': 2600, 'sigma': 50, 'weight': 0., 'prior': stats.uniform(loc=520, scale=2600)}}

initial = {'c1': 0.8, 'c2': 1.0, 'p': 0.5, 'rho': 1.0, 'N': 50000, 't.end': 30.*52}

fitter = Fitter(settings, initial)


# read in target tree
try:
    fn = sys.argv[1]
    fitter.set_target_tree(fn)
except:
    print '\nUsage: python colgem_fitter.py [newick]'
    sys.exit()


# prevent previous log files from being overwritten
logfn = fn.replace('.nwk', '.log')
modifier = ''
tries = 0
while os.path.exists(logfn+modifier):
    tries += 1
    modifier = '.%d' % tries

logfn += modifier

print 'writing to', logfn

logfile = open(logfn, 'w')
fitter.abc_mcmc(logfile, skip=10, tol0=0.005, mintol=0.001, decay=0.0025)
logfile.close()


