kamphir
=======

KAMPHIR stands for Kernel-assisted ABC-MCMC for PHylodynamic InfeRence.
Approximate Bayesian computing (ABC) attempts to fit a model to an observed data set by simulating additional data sets under different parameter settings of the model until the simulations resemble the observations.
ABC-MCMC is an implementation of ABC using Markov chain Monte Carlo sampling where new parameter values are drawn from a proposal distribution centered at the current state, and used to update the model if the resulting simulations meet some goodness-of-fit criterion. 
A kernel method 

##Dependencies
* [Python](https://www.python.org/) - Kamphir was developed with Python 2.7.  Several required modules are only available in distributions of Python since version 2.6, such as [json](https://docs.python.org/2/library/json.html).
* [Biopython](http://biopython.org/wiki/Main_Page) - A collection of tools for working with biological data.  Kamphir makes extensive use of the Phylo module in Biopython for handling tree objects.
* [NumPy](http://www.numpy.org/) - A package for scientific computing in Python.  Kamphir makes use of its Array objects for improved performance.
* [dill](https://pypi.python.org/pypi/dill) - Python module that extends pickle module for serializing Python objects
* [jinja2](http://jinja.pocoo.org/) - Python module for populating templates with Python objects.
 
##Requires at least one of:
* [MASTER](http://compevol.github.io/MASTER/) - Java module for reaction system-based simulation of epidemiological processes in the BEAST2 package
* [rcolgem](http://colgem.r-forge.r-project.org/) - R module for coalescent simulation and inference for epidemiological models through numerical solution of ODEs

