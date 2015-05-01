"""
Parse kamphir log to extract parameter estimates at regular intervals
and use these to output population trajectories and simulate trees.

Only rcolgem is supported for now.
"""

import rcolgem
import argparse

class PostKamphir():
    def __init__ (self, logfile):
        self.logdata = {}
        header = []
        for line in logfile:
            if line.startswith('#'):
                continue
            items = line.strip('\n').split('\t')

            # use header to prepare container
            if len(header) == 0:
                header = items
                for key in header:
                    self.logdata.update({key: []})
                continue

            for i, key in enumerate(header):
                self.logdata[key].append(float(items[i]))

    def analyze(self):
        """
        Call Rcolgem
        :return:
        """






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='KAMPHIR-post',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('log', help='<INPUT> Kamphir log for post-processing')
    parser.add_argument('model', help='Rcolgem model used to generate log',
                        choices=['*', 'SI', 'SI2', 'DiffRisk', 'Stages'])
    parser.add_argument('nwk', help='<OUTPUT> file to write Newick tree strings')
    parser.add_argument('csv', help='<OUTPUT> file to write trajectories as CSV')

    parser.add_argument('-ntrees', default=100, help='number of trees to output')
    parser.add_argument('-nrows', default=100, help='number of trajectories to output')

    args = parser.parse_args()

    # open handle to log file
    with open(args.log, 'rU') as handle:
        pk = PostKamphir(logfile=handle)

