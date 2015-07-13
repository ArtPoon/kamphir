#!/usr/bin/env python

# For now, this has to be run from the kamphir root directory, not from
# the tests directory.

# TODO: this is an ugly hack, fix it
import sys
import os
sys.path.append(os.getcwd())

import unittest
import kamphir
import dill
import operator
import random
import math
import json
from cStringIO import StringIO
from Bio import Phylo

import phyloK2
import rcolgem
from rpy2 import robjects

# TODO: testing of multiprocessing should be optional
import multiprocessing

# make a phony class for testing apply on class methods
class Dummy:
    def multiply(self, x, y):
        """A deterministic function of x and y"""
        return x*y

    def random_multiply(self, x, y):
        """A non-deterministic function of x and y"""
        return x*y*random.random()


class ParallelTests(unittest.TestCase):
    """Basic tests for the parallelization functions"""

    def setUp(self):
        # set up a multiprocessing pool
        self.nthreads = 4
        self.pool = multiprocessing.Pool(processes=self.nthreads)

        # make an instance of our dummy class
        self.dummy = Dummy()

    def test_run_dill_encoded(self):
        """Test basic running of a dill-encoded function"""
        what = dill.dumps((operator.add, [2, 3]))
        self.assertEqual(kamphir.run_dill_encoded(what), 5)

    def test_apply_async(self):
        """Test dispatching of a task to the multiprocessing pool"""
        result = kamphir.apply_async(self.pool, operator.mul, [2, 3])
        multiprocessing.pool.ApplyResult.wait(result)
        self.assertEqual(result.get(), 6)

    def test_apply_async_classmethod(self):
        """Test apply_async on a class method"""
        args = [ (2, x) for x in range(10) ]
        result = [kamphir.apply_async(self.pool, self.dummy.multiply, x) for x in args]
        [multiprocessing.pool.ApplyResult.wait(x) for x in result]
        result = [x.get() for x in result]
        expected = [ 2*x for x in range(10) ]
        self.assertListEqual(result, expected)

    def test_apply_async_random(self):
        """Test apply_async on a non-deterministic class method"""
        args = [ (3, x) for x in range(10) ]
        result = [kamphir.apply_async(self.pool, self.dummy.random_multiply, x) for x in args]
        [multiprocessing.pool.ApplyResult.wait(x) for x in result]
        result = [x.get() for x in result]
        # make sure no two are the same
        self.assertListEqual(sorted(list(set(result))), sorted(result)) 

    def tearDown(self):
        self.pool.close()


class TreeKernelTests(unittest.TestCase):
    """Tests for the tree kernel"""

    def setUp(self):
        T1_str = "( ( A:0.5, B:0.25 )E:0.5, ( C:0.25, D:0.25 )F:0.5 )G;"
        T2_str = "( ( ( A:0.25, B:0.25 )E:0.5, C:0.25 )F:0.5, D:0.25 )G;"
        self.T1 = Phylo.read(StringIO(T1_str), "newick")
        self.T2 = Phylo.read(StringIO(T2_str), "newick")
        self.kernel = phyloK2.PhyloKernel(decayFactor=0.5, gaussFactor=1)

    def test_kernel(self):
        """Test that the tree kernel works on a simple example"""
        self.assertEqual(self.kernel.kernel(self.T1, self.T2), 1.125 * (1+math.exp(-0.0625)))

    def test_kernel_large(self):
        """Test that the tree kernel works on a larger example"""
        t = list(Phylo.parse("tests/test.tree", "newick"))
        k = self.kernel.kernel(t[0], t[1])
        self.assertAlmostEqual(k, 3042.58677467)

    def test_normalize_tree_median(self):
        expected_bl = [c.branch_length/0.375 for c in self.T1.find_clades(order="postorder")][:-1]
        self.kernel.normalize_tree(self.T1)
        bl = [c.branch_length for c in self.T1.find_clades(order="postorder")][:-1]
        self.assertListEqual(bl, expected_bl)

    def test_normalize_tree_mean(self):
        expected_bl = [c.branch_length/0.375 for c in self.T1.find_clades(order="postorder")][:-1]
        self.kernel.normalize_tree(self.T1, "mean")
        bl = [c.branch_length for c in self.T1.find_clades(order="postorder")][:-1]
        self.assertListEqual(bl, expected_bl)

class RcolgemTests(unittest.TestCase):
    """Tests for the Rcolgem simulation functions"""

    def setUp(self):
        self.r = rcolgem.Rcolgem(ncores=2, nreps=2, t0=0, fgy_resolution=500, integration_method='rk4', seed=0)
        self.params = { "N": 100, "beta": 0.05, "beta1": 0.05, "beta2": 0.05, "gamma": 0.005, "mu": 0.005, "lambd": 0.005 , "t_break": 50,
                        "c1": 1, "c2": 1, "rho": 0, "p": 0.5 }
        self.tree_height = 100
        self.tip_heights = [0] * 10

    def test_rcolgem_init(self):
        """Test that rcolgem is initialized as expected"""
        # make sure the needed variables were set when rcolgem was initialized in setUp
        self.assertEqual(robjects.r["n.cores"][0], 2)
        self.assertEqual(robjects.r["nreps"][0], 2)
        self.assertEqual(robjects.r["fgyResolution"][0], 500)
        self.assertEqual(robjects.r["integrationMethod"][0], "rk4")
        self.assertEqual(robjects.r["t0"][0], 0)

        # check that a cluster was initialized
        self.assertEqual(robjects.r('class(cl)')[1], "cluster")

    def test_simulate_SI_tree(self):
        """Test that an SI tree is correctly simulated"""
        self.r.init_SI_model()
        trees = self.r.simulate_SI_trees(self.params, self.tree_height, self.tip_heights)
        self.assertEqual(trees[0], "((3:30.77470915,5:30.77470915):44.54162798,(((7:6.293201,10:6.293201):27.90419663,(8:33.89050516,(1:25.09241256,4:25.09241256):8.798092599):0.3068924695):3.362019587,(2:34.72688185,(6:19.87778928,9:19.87778928):14.84909257):2.832535369):37.75691991);")
        self.assertEqual(trees[1], "(8:72.29052237,(5:65.24150433,((10:18.22316917,(1:10.60929307,(2:1.100461631,4:1.100461631):9.508831439):7.613876095):46.73482677,(3:56.22533304,(7:44.77487105,(6:33.05115908,9:33.05115908):11.72371197):11.450462):8.732662893):0.283508395):7.049018038);")

    def test_simulate_SI2_tree(self):
        """Test that SI2 reduces to SI when beta1 == beta2"""
        self.r.init_SI_model()
        trees = self.r.simulate_SI2_trees(self.params, self.tree_height, self.tip_heights)
        self.assertEqual(trees[0], "((3:30.77470915,5:30.77470915):44.54162798,(((7:6.293201,10:6.293201):27.90419663,(8:33.89050516,(1:25.09241256,4:25.09241256):8.798092599):0.3068924695):3.362019587,(2:34.72688185,(6:19.87778928,9:19.87778928):14.84909257):2.832535369):37.75691991);")
        self.assertEqual(trees[1], "(8:72.29052237,(5:65.24150433,((10:18.22316917,(1:10.60929307,(2:1.100461631,4:1.100461631):9.508831439):7.613876095):46.73482677,(3:56.22533304,(7:44.77487105,(6:33.05115908,9:33.05115908):11.72371197):11.450462):8.732662893):0.283508395):7.049018038);")

    def test_simulate_DiffRisk_tree(self):
        self.r.init_DiffRisk_model()
        robjects.r("set.seed(0)")
        robjects.r("clusterSetRNGStream(cl, 0)")
        trees = self.r.simulate_DiffRisk_trees(self.params, self.tree_height, self.tip_heights)
        self.assertEqual(trees[0], "((1:34.70552976,(5:30.75607927,7:30.75607927):3.949450494):40.52923452,((3:19.8658362,8:19.8658362):17.67013076,((2:6.28920109,10:6.28920109):27.8872289,(9:33.86974465,(4:25.07737022,6:25.07737022):8.792374435):0.3066853327):3.359536976):37.69879732);")
        self.assertEqual(trees[1], "((9:56.18259136,(4:44.74529922,6:44.74529922):11.43729214):16.03444887,((5:18.2121653,(8:10.60270091,10:10.60270091):7.609464388):46.97152823,(1:64.90076175,(7:33.030997,(2:1.099747134,3:1.099747134):31.93124987):31.86976475):0.282931775):7.0333467);")


class KamphirTests(unittest.TestCase):
    def setUp(self):
        with open("settings/example.DiffRisk.json") as f:
            self.settings = json.load(f)
        self.driver = None
        self.rcolgem = rcolgem.Rcolgem(ncores=1, nreps=1, t0=0, fgy_resolution=100, integration_method='rk4')
        self.simfunc = self.rcolgem.simulate_DiffRisk_trees
        self.kamphir = kamphir.Kamphir(self.settings, None, None, self.simfunc, ncores=2, nreps=2)

        T1_str = "( ( A:0.5, B:0.25 )E:0.5, ( C:0.25, D:0.25 )F:0.5 )G:0;"
        T2_str = "( ( ( A:0.25, B:0.25 )E:0.5, C:0.25 )F:0.5, D:0.25 )G:0;"
        self.T1 = Phylo.read(StringIO(T1_str), "newick")
        self.T2 = Phylo.read(StringIO(T2_str), "newick")

    def test_set_target_trees(self):
        self.kamphir.set_target_trees("tests/test.tree", 0)
        tree = next(Phylo.parse("tests/test.tree", "newick"))
        height = max(tree.depths().values())
        node_heights = [height-x[1] for x in tree.depths().items() if not x[0].is_terminal()]

        self.assertEqual(max(tree.depths().values()), self.kamphir.target_trees[0][1])
        self.assertListEqual([0.] * 100, self.kamphir.target_trees[0][2])
        self.assertListEqual(sorted(node_heights), self.kamphir.target_trees[0][3])

        tree.ladderize()
        tree.root.branch_length = 0
        self.kamphir.normalize_tree(tree, "mean")
        self.kamphir.annotate_tree(tree)
        self.assertEqual(tree.__format__("newick"), self.kamphir.target_trees[0][0].__format__("newick"))
        self.assertEqual(self.kamphir.kernel(tree, tree), self.kamphir.target_trees[0][4])

if __name__ == "__main__":
    unittest.main()
