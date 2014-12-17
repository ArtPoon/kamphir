"""
Use XML template to generate control file for generating
a tree under an SI model with MASTER 2.0 (BEAST).
"""
import jinja2
import os
import sys
from Bio import Phylo
from random import sample
import subprocess

jarfile = '/Users/art/src/MASTER-2.0.0/dist/MASTER-2.0.0/MASTER-2.0.0.jar'
FNULL = open(os.devnull, 'w')

try:
    infile = sys.argv[1]
    tipfile = sys.argv[2]
    outfile = sys.argv[3]
except:
    print 'Usage: python MASTER.SIR.py [input CSV] [tip labels CSV] [output NWK]'
    sys.exit(1)

# get parent process ID from filename
pid = int(infile.split('.')[0].split('_')[-1])
tmpfile = '/tmp/MASTER.SIR.%d.xml' % pid

jenv = jinja2.Environment(
    block_start_string='{%',
    block_end_string='%}',
    variable_start_string='{{',
    variable_end_string='}}',
    loader=jinja2.FileSystemLoader(os.getcwd())
)

template = jenv.from_string(source=
"""
<beast version='2.0' namespace='master
                                :master.model
                                :master.conditions
                                :master.outputs
                                :master.postprocessors'>
    <run spec='InheritanceEnsemble'
         nTraj='{{ nreps|int }}'
         samplePopulationSizes="true"
         verbosity="1"
         simulationTime="{{ t_end }}">

        <model spec='Model' id='model'>
            <population spec='Population' id='S' populationName='S'/>
            <population spec='Population' id='I' populationName='I'/>
            <population spec='Population' id='R' populationName='R'/>
            <population spec='Population' id='I_sample' populationName='I_sample'/>

            <reaction spec='Reaction' reactionName="Infection" rate="{{ beta }}">
                S + I -> 2I
            </reaction>
            <reaction spec='Reaction' reactionName="Recovery" rate="{{ gamma }}">
                I -> R
            </reaction>
            <reaction spec='Reaction' reactionName="Sampling" rate="{{ phi }}">
                I -> I_sample
            </reaction>
        </model>

        <initialState spec='InitState'>
            <populationSize spec='PopulationSize' population='@S' size='{{ N - 1 }}'/>
            <lineageSeed spec='Individual' population='@I'/>
        </initialState>

        <inheritancePostProcessor spec="LineageFilter" reactionName="Sampling"/>
        <postSimCondition spec="LeafCountPostSimCondition" nLeaves="{{ ntips }}" exact="false"/>

        <output spec='NewickOutput' fileName='{{ outfile }}'/>
    </run>
</beast>
""")

# default parameter values
context = {
    'beta': 0.001,  # transmission rate
    'gamma': 0.3,  # mortality rate
    'N': 1000,      # total population size
    't_end': 30,    # length of simulation
    'phi': 0.15,    # sampling rate
    'ntips': 100,    # number of tips in tree
    'nreps': 10,     # number of trees to generate
    'outfile': outfile
}


# read parameter settings from input CSV
handle = open(infile, 'rU')
for line in handle:
    key, value = line.strip('\n').split(',')
    context.update({key: float(value)})
handle.close()

# FIXME: MASTER tends to generate larger trees than requested
# FIXME: setting post filter to "exact" is extremely inefficient
# FIXME: as a quick work-around, request same number of tips as simulation
# infer number of tips from tip label CSV
handle = open(tipfile, 'rU')
context['ntips'] = len(handle.readlines())
handle.close()

# populate template from context
handle = open(tmpfile, 'w')
handle.write(template.render(context))
handle.close()

# call MASTER
#os.system('master2 %s > /dev/null' % tmpfile)
subprocess.call(['java', '-jar', jarfile, tmpfile], shell=True, stdout=FNULL, stderr=FNULL)


# sample tips to enforce size of tree
trees = Phylo.parse(outfile, 'newick')
trees2 = []
for tree in trees:
    tips = tree.get_terminals()
    tips2 = sample(tips, context['ntips'])
    for tip in tips:
        tip.name = str(tip.confidence)
        if tip in tips2:
            continue
        _ = tree.prune(tip)
    trees2.append(tree)

Phylo.write(trees2, outfile, 'newick')
