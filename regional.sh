#!/bin/bash

export R_HOME="/usr/local/lib/R"
export LD_LIBRARY_PATH=$R_HOME/lib:$LD_LIBRARY_PATH

SCENARIO=$1
TREE_FILE=/home/rmcclosk/Documents/pangea/data/February2015/Regional/hyphy/root2tip/150129_PANGEAsim_Regional_FirstObj_sc${SCENARIO}_SIMULATED_all.nwk
SETTINGS_FILE=/home/rmcclosk/Documents/kamphir/settings/pangea_sc${SCENARIO}.json
LOG_FILE=sc${SCENARIO}-final.log
KAMPHIR_ARGS="-delimiter _ -ncores 5 -nthreads 5 -nreps 5 -tscale 0.142857 -tol0 0.03 -mintol 0.01 -prior"

while ((1)); do
    if [[ -f ${LOG_FILE} ]]; then
        echo "Restarting ${LOG_FILE} at `date`"
        LAST_LOG=`ls -1 ${LOG_FILE}* | sort -n -k 3 -t '.' | tail -n 1`
        python kamphir.py PANGEA ${SETTINGS_FILE} ${TREE_FILE} ${LOG_FILE} ${KAMPHIR_ARGS} -restart ${LAST_LOG}
    else
        python kamphir.py PANGEA ${SETTINGS_FILE} ${TREE_FILE} ${LOG_FILE} ${KAMPHIR_ARGS}
    fi
    
    LAST_LOG=`ls -1 ${LOG_FILE}* | sort -n -k 3 -t '.' | tail -n 1`
    if [[ `wc -l < ${LAST_LOG} ` -le 7 ]]; then
        echo "Removing unsuccesful log file ${LAST_LOG} at `date`"
        rm ${LAST_LOG}
        rm $(echo ${LAST_LOG} | sed s/log/trees/)
    fi
done
