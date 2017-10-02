#! /usr/bin/env bash

#set -x

export LD_LIBRARY_PATH=$LIBLOCATION
export DYLD_LIBRARY_PATH=$LIBLOCATION

echo PYTHONPATH=$PYTHONPATH
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
echo PATH=$PATH
echo
echo RIVET_ANALYSIS_PATH=$RIVET_ANALYSIS_PATH
echo RIVET_DATA_PATH=$RIVET_DATA_PATH
echo RIVET_REF_PATH=$RIVET_REF_PATH
echo RIVET_INFO_PATH=$RIVET_INFO_PATH
echo
echo PYTHON=$PYTHON


echo "trying to load rivet python module"
$PYTHON -c 'import rivet'  || exit $?
echo "Success"
