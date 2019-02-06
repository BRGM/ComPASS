#!/bin/bash
#
# WARNING check that this file hase Unix syle EOL
#

set -e

for script in *.py
do
    echo $script
    # redirect standard flux to file
    time (python3 $script > out-$script)
    # redirect standard and error fluxes to file
    # time (python3 $script > out-$script 2>&1)
    echo
done
