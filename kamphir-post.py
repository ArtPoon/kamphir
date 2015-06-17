"""
Parse kamphir log to extract parameter estimates at regular intervals
and use these to output population trajectories and simulate trees.

Only rcolgem is supported for now.
"""
import sys
from rcolgem import Rcolgem
import argparse
from Bio import Phylo
from math import floor

def post_process(logfile, tree_height, tip_heights, model, ntrees, nrows, resol, burnin, nwkfile, csvfile):
    # select simulation function
    rcolgem = Rcolgem(ncores=1, nreps=1, fgy_resolution=resol)
    if model == 'SI':
        rcolgem.init_SI_model()
        simfunc = rcolgem.simulate_SI_trees
    elif model == 'SI2':
        rcolgem.init_SI_model()
        simfunc = rcolgem.simulate_SI2_trees
    elif model == 'DiffRisk':
        rcolgem.init_DiffRisk_model()
        simfunc = rcolgem.simulate_DiffRisk_trees
    elif model == 'Stages':
        rcolgem.init_stages_model()
        simfunc = rcolgem.simulate_stages_trees
    else:
        print 'ERROR: Unrecognized model', model
        sys.exit()

    # parse log data
    logdata = {}
    header = []
    ndata = 0
    for line in logfile:
        if line.startswith('#'):
            continue
        items = line.strip('\n').split('\t')

        # use header to prepare container
        if len(header) == 0:
            header = items
            for key in header:
                logdata.update({key: []})
            continue

        ndata += 1  # the total number of data rows
        if ndata <= burnin:
            continue  # discard step
        for i, key in enumerate(header):
            logdata[key].append(float(items[i]))


    maxrow = len(logdata.values()[0])

    # determine steps that we will output
    nwksteps = range(maxrow)  # return all available
    if ntrees < maxrow:
        step = float(maxrow)/ntrees
        nwksteps = [maxrow - int(round(i*step))-1 for i in range(ntrees)]

    csvsteps = range(maxrow)
    if nrows < maxrow:
        step = float(maxrow)/nrows
        csvsteps = [maxrow - int(round(i*step))-1 for i in range(nrows)]

    csvheader = False
    for step in range(maxrow):
        if step not in nwksteps and step not in csvsteps:
            continue
        print '(%d/%d)' % (step, maxrow)

        # extract parameter vector from log data
        params = {}
        for key, vals in logdata.iteritems():
            params.update({key: vals[step]})

        # solve ODE and simulate tree
        try:
            trees, trajectories = simfunc(params, tree_height, tip_heights, post=True)
        except:
            print 'Computation failed, skipping step', step
            sys.exit()

        if step in csvsteps:
            # output trajectories
            if csvheader is False:
                # output column names once only
                csvfile.write(','.join(map(str, ['step']+[x for x in trajectories.colnames])))
                csvfile.write('\n')
                csvheader = True
            for i in range(trajectories.nrow):
                row = trajectories.rx(i+1, True)  # R is 1-indexed
                csvfile.write(','.join(map(str, [step]+[x for x in row])))
                csvfile.write('\n')

        if step in nwksteps:
            # sometimes fgyResolution is set too low to simulate a tree
            nwkfile.write('' if len(trees) == 0 else trees[0])
            nwkfile.write('\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='KAMPHIR-post',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('log', help='<INPUT> Kamphir log for post-processing')
    parser.add_argument('tree', help='<INPUT> Newick tree string used to fit model')
    parser.add_argument('model', help='Rcolgem model used to generate log',
                        choices=['SI', 'SI2', 'DiffRisk', 'Stages'])
    parser.add_argument('nwk', help='<OUTPUT> file to write Newick tree strings')
    parser.add_argument('csv', help='<OUTPUT> file to write trajectories as CSV')

    parser.add_argument('-burnin', type=int, default=100, help='Number of steps to skip as burnin.')
    parser.add_argument('-ntrees', type=int, default=100, help='Number of trees to output.')
    parser.add_argument('-nrows', type=int, default=100, help='Number of trajectories to output.')
    parser.add_argument('-resol', type=int, default=100, help='Resolution for numerical solution of ODE.')
    parser.add_argument('-delimiter', default=None, help='Field separator for node names in tree.')
    parser.add_argument('-position', type=int, default=-1, help='Python index of field with tip date.')

    args = parser.parse_args()

    # open and parse tree
    try:
        tree = Phylo.read(args.tree, 'newick')
    except:
        print 'ERROR: Failed to parse tree from file', args.tree
        raise

    tree_height = max(tree.depths().values())
    tips = tree.get_terminals()
    ntips = len(tips)

    if args.delimiter is None:
        tip_heights = [0.] * ntips
    else:
        maxdate = 0
        tipdates = []
        for tip in tips:
            try:
                items = tip.name.strip("'").split(args.delimiter)
                tipdate = float(items[args.position])
                if tipdate > maxdate:
                    maxdate = tipdate
            except:
                print 'Warning: Failed to parse tipdate from label', tip.name
                tipdate = None  # gets interpreted as 0
                pass

            tipdates.append(tipdate)

        tip_heights = [str(maxdate-t) if t else 0 for t in tipdates]

    # prepare outputs
    nwkfile = open(args.nwk, 'w')
    csvfile = open(args.csv, 'w')

    # open handle to log file
    with open(args.log, 'rU') as handle:
        post_process(logfile=handle, tree_height=tree_height, tip_heights=tip_heights,
                     model=args.model, ntrees=args.ntrees, nrows=args.nrows, resol=args.resol, burnin=args.burnin,
                     nwkfile=nwkfile, csvfile=csvfile)

