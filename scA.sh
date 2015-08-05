#!/bin/bash
python kamphir.py Stages /home/rmcclosk/Documents/pangea/settings/regional.Stages.json /home/rmcclosk/Documents/pangea/data/February2015/Regional/hyphy/root2tip/150129_PANGEAsim_Regional_FirstObj_scA_SIMULATED_all.nwk scA-stages.log -delimiter _ -ncores 5 -nthreads 5 -nreps 5 -tscale 0.142857 -tol0 0.03 -mintol 0.01
