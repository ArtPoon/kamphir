"""
Populate BEAST2 XML template with new sequence data from FASTA
"""

from xml.etree.ElementTree import ElementTree as Tree
from xml.etree.ElementTree import Element as Node
from seqUtils import convert_fasta
import sys
import random

try:
    template_file = sys.argv[1]
    infile = sys.argv[2]
    outfile = sys.argv[3]
    nseqs = int(sys.argv[4])
except:
    print 'Usage: python populate_beast2_xml.py [template XML] [input FASTA] [output XML] [#seqs]'
    sys.exit()

if not outfile.endswith('.xml'):
    outfile += '.xml'

handle = open(infile, 'rU')
fasta = convert_fasta(handle)
handle.close()

if nseqs > len(fasta):
    print 'WARNING: requested #seqs exceeds length of FASTA (', len(fasta), '), truncating'
    nseqs = len(fasta)

# extract tip dates

template = Tree()
root = template.parse(template_file)

data = template.findall('data')[0]
data._children = []  # reset data block

sample = random.sample(fasta, nseqs)
tipdates = ''
for h, s in sample:
    tipname = h.strip()
    tipdate = tipname.split('_')[-1]
    tipdates += '%s = %s, ' % (tipname, tipdate)
    seq = Node('sequence', {'totalcount': '4', 'id': 'seq_seq_'+tipname, 'value': s,
                            'taxon': tipname})
    data._children.append(seq)


# replace tip date information
traits = template.find('trait')
traits.set('value', tipdates.strip(', '))

# update log file paths
run = template.find('run')
for logger in run.findall('logger'):
    if 'fileName' in logger.keys():
        if logger.get('mode', None) == 'tree':
            logger.set('fileName', outfile.replace('.xml', '.trees'))
        else:
            logger.set('fileName', outfile.replace('.xml', '.log'))

# output new XML
template.write(outfile)

