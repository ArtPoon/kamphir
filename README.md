This repository comprises scripts and data that were used in a recent paper:

Poon AF. Phylodynamic inference with kernel ABC and its application to HIV epidemiology. Molecular biology and evolution. 2015 May 25:msv123.

to apply approximate Bayesian computation (ABC) to phylodynamic inference using a kernel method to compare the shapes of trees.

We are presently working on making a stable implementation of these concepts [here](http://github.com/PoonLab/Kaphi).
This current work replaces the MCMC sampling method with sequential Monte Carlo (SMC), with which we have obtained more consistent results (see [netabc](http://github.com/rmcclosk/netabc)) and published here:

McCloskey RM, Liang RH, Poon AF. Reconstructing contact network parameters from viral phylogenies. Virus evolution. 2016 Jul 1;2(2):vew029.


##Dependencies
* [Python](https://www.python.org/) - Kamphir was developed with Python 2.7.  Several required modules are only available in distributions of Python since version 2.6, such as [json](https://docs.python.org/2/library/json.html).
* [Biopython](http://biopython.org/wiki/Main_Page) - A collection of tools for working with biological data.  Kamphir makes extensive use of the Phylo module in Biopython for handling tree objects.
* [NumPy](http://www.numpy.org/) - A package for scientific computing in Python.  Kamphir makes use of its Array objects for improved performance.
* [dill](https://pypi.python.org/pypi/dill) - Python module that extends pickle module for serializing Python objects
* [jinja2](http://jinja.pocoo.org/) - Python module for populating templates with Python objects.
* [R](http://www.r-project.org/) - Widely used language for statistical computing.
* [rpy2](http://rpy.sourceforge.net/) - Python interface to R.
 
##Requires at least one of:
* [MASTER](http://compevol.github.io/MASTER/) - Java module for reaction system-based simulation of epidemiological processes in the BEAST2 package
* [rcolgem](http://colgem.r-forge.r-project.org/) - A modified version of this R module is already incorporated into this repository.  This module performs coalescent simulation and inference for epidemiological models through numerical solution of ODEs.  Requires [R](http://cran.r-project.org/) and packages [ape](http://cran.r-project.org/web/packages/ape/index.html), [deSolve](http://cran.r-project.org/web/packages/deSolve/index.html), [bbmle](http://cran.r-project.org/web/packages/bbmle/index.html)

