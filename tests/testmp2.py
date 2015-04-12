import sys
import scipy.stats as stats
sys.path.append('..')  # make modules in parent dir available
import time
from kamphir import Kamphir

settings = {'c1': {'initial': 0.8,
                   'min': 0.1,
                   'max': 10.,
                   'sigma': 0.05,
                   'weight': 1.,
                   'prior': stats.lognorm(1)},
            'c2': {'initial': 1.0,
                   'min': 0.1,
                   'max': 10.,
                   'sigma': 0.25,
                   'weight': 1.,
                   'prior': stats.lognorm(0.5)},
            'p': {'initial': 0.5,
                  'min': 0.,
                  'max': 1.,
                  'sigma': 0.05,
                  'weight': 1.,
                  'prior': stats.uniform(loc=0, scale=1)},
            'rho': {'initial': 0.9,
                    'min': 0.,
                    'max': 1.,
                    'sigma': 0.05,
                    'weight': 1.,
                    'prior': stats.uniform(loc=0, scale=1)},
            'N': {'initial': 10000,
                  'min': 1000,
                  'max': 1e8,
                  'sigma': 5000,
                  'weight': 1.,
                  'prior': stats.norm(loc=8000, scale=3000)},
            't.end': {'initial': 30*52.,
                      'min': 520,
                      'max': 2600,
                      'sigma': 50,
                      'weight': 0.,  # fixed
                      'prior': stats.uniform(loc=520, scale=2600)}}

kam = Kamphir(settings=settings, rscript='../simulate.DiffRisk.R', ncores=6, nthreads=1)
kam.set_target_tree('large1.nwk')

t0 = time.time()
trees = list(kam.simulate_external())
print 'required', time.time() - t0, 'seconds to simulate trees'

# single-process
t0 = time.time()
res = kam.evaluate(trees=trees)
elapsed = time.time() - t0

print 'SP', res, elapsed, 'seconds'

# multi-process
t0 = time.time()
res = kam.evaluate(trees=trees, nthreads=12)
elapsed = time.time() - t0

print 'MP', res, elapsed, 'seconds'


