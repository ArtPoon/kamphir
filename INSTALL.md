Installing rpy2 
---------------

Kamphir requires a newer version of `rpy2` than the one available
through `apt-get`. First try installing through `pip`:

    pip install rpy2

This command may fail with an error about being unable to find `R.h`. If
this happens, you can install `rpy2` from source. Download and extract
the source tarball. If you are on OS X, replace `wget` below with
`curl`.

    wget https://pypi.python.org/packages/source/r/rpy2/rpy2-2.5.6.tar.gz
    tar xf rpy2-2.5.6.tar.gz
    cd rpy2-2.5.6
    sudo python ./setup.py install

If this fails with the same error, you may have to reinstall R from
source. Remove all installed R packages (like `r-base`, `r-base-core`,
and `r-base-dev`), then download the R source tarball from your
preferred mirror. As above, extract the archive, `cd` into it, and do

    ./configure --enable-R-shlib
    make
    sudo make install

Make sure to include the `--enable-R-shlib` flag. After installation,
you sholud be able to install `rpy2` from source.

Setting environment variables
-----------------------------

R needs to find shared libraries in order to work properly. Find out
where R installed them with `locate`.

    locate libRlapack.so

Some of the paths produced here may be in the R build directory if you
installed R from source, but there should be at least one in the `/usr`
or `/opt` hierarchy. On Ubuntu, this was `/usr/local/lib/R/lib`. Put
a line in your `.bash_profile` indicating where the shared libraries are
to be found.

    export LD_LIBRARY_PATH=/usr/local/lib/R/lib:$LD_LIBRARY_PATH

Installing rcolgem
------------------

Enter the R console (as root, unless you have a local R library set up),
and install the `deSolve` package.

    install.packages("deSolve")
    install.packages("ape")

Exit R. Enter the `rcolgem` subdirectory, and type these commands.

    R CMD build pkg
    sudo R CMD INSTALL rcolgem_0.0.4.tar.gz

The second command doesn't need to be run with `sudo` if you are
installing locally, but you will need to have a local R library set up.
