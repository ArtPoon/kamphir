"""
Label tips with heights
"""
import sys
from Bio import Phylo

try:
    tree = Phylo.read(sys.argv[1], 'newick')
    outfile = open(sys.argv[2], 'w')
except:
    print 'python label_tips [input NWK] [output NWK]'
    sys.exit()

depths = tree.depths()
tips = tree.get_terminals()
for tip in tips:
    depth = depths[tip]
    tip.name = '%s_%1.6f' % (tip.confidence, depth)
    
Phylo.write(tree, outfile, 'newick')
