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
from cStringIO import StringIO
from Bio import Phylo

import phyloK2

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


class TreeKernelTests(unittest.TestCase):
    """Tests for the tree kernel"""

    def setUp(self):
        T1_str = "( ( A:0.5, B:0.25 )E:0.5, ( C:0.25, D:0.25 )F:0.5 )G;"
        T2_str = "( ( ( A:0.25, B:0.25 )E:0.5, C:0.25 )F:0.5, D:0.25 )G;"
        self.T1 = Phylo.read(StringIO(T1_str), "newick")
        self.T2 = Phylo.read(StringIO(T2_str), "newick")

    def test_kernel(self):
        """Test that the tree kernel works on a simple example"""
        kernel = phyloK2.PhyloKernel(decayFactor=0.5, gaussFactor=1)
        self.assertEqual(kernel.kernel(self.T1, self.T2), 1.125 * (1+math.exp(-0.03125)))

if __name__ == "__main__":
    unittest.main()
