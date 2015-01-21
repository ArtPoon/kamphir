kamphir
=======

**KAMPHIR** stands for Kernel-assisted ABC-MCMC for PHylodynamic InfeRence.

Approximate Bayesian computing (ABC) attempts to fit a model to an observed data set by simulating additional data sets under different parameter settings of the model until the simulations resemble the observations.
ABC-MCMC is an implementation of ABC using Markov chain Monte Carlo sampling where new parameter values are drawn from a proposal distribution centered at the current state, and used to update the model if the resulting simulations meet some goodness-of-fit criterion. 

A kernel method is a way of comparing complex objects.  Objects can be mapped to a feature space by breaking each object down to its constituent parts.  We could then calculate the inner product of two objects in this feature space by multiplying the counts of each feature in their respective vectors.  An inner product provides a similarity measure of the two objects - it takes its maximum value when the vectors are the same.  However, the number of possible features can easily become enormous (for example, the number of all possible words in the English language, including "jabberwocky" and "w00t").  The kernel method is essentially an efficient way of computing this inner product by summing only features that appear in at least one of the objects.  

Kamphir uses a kernel function that compares the shapes of phylogenetic trees.  In this case, the features are subset trees extracted from each phylogeny.  When computing the inner product, these features are weighted by the discordance in branch lengths, which represent the passage of evolutionary or chronological time (depending on how the phylogeny was reconstructed).  It turns out that this kernel function provides a rather effective similarity measure for use in the ABC framework.

The take home message is that Kamphir enables you to fit practically any model to a phylogeny, so long as that model can be used to simulate trees.  It requires that you have a program that will generate trees, and a driver script that will automate the process of calling this program under different model parameter values.  I have provided driver scripts for two such programs, [MASTER](http://compevol.github.io/MASTER/) and [rcolgem](http://colgem.r-forge.r-project.org/).

##Dependencies
* [Python](https://www.python.org/) - Kamphir was developed with Python 2.7.  Several required modules are only available in distributions of Python since version 2.6, such as [json](https://docs.python.org/2/library/json.html).
* [Biopython](http://biopython.org/wiki/Main_Page) - A collection of tools for working with biological data.  Kamphir makes extensive use of the Phylo module in Biopython for handling tree objects.
* [NumPy](http://www.numpy.org/) - A package for scientific computing in Python.  Kamphir makes use of its Array objects for improved performance.
* [dill](https://pypi.python.org/pypi/dill) - Python module that extends pickle module for serializing Python objects
* [jinja2](http://jinja.pocoo.org/) - Python module for populating templates with Python objects.
 
##Requires at least one of:
* [MASTER](http://compevol.github.io/MASTER/) - Java module for reaction system-based simulation of epidemiological processes in the BEAST2 package
* [rcolgem](http://colgem.r-forge.r-project.org/) - R module for coalescent simulation and inference for epidemiological models through numerical solution of ODEs

