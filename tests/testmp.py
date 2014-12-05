import sys
sys.path.append('..')  # make modules in parent dir available
import time
import phyloK2
from Bio import Phylo

pk = phyloK2.PhyloKernel()

t1 = Phylo.read('large1.nwk', 'newick')
t2 = Phylo.read('large2.nwk', 'newick')
t1.ladderize()
t2.ladderize()
pk.normalize_tree(t1)
pk.normalize_tree(t2)
pk.annotate_tree(t1)
pk.annotate_tree(t2)

# single-threaded
t0 = time.time()
k = pk.kernel(t1, t2)
elapsed = time.time() - t0

print 'SP', k, elapsed, 'seconds'

t0 = time.time()
k = pk.kernel_parallel(t1, t2, 6)
elapsed = time.time() - t0

print 'MP', k, elapsed, 'seconds'
