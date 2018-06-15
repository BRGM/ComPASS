#!/bin/sh

#TODO perform control of Environment variables to prevent from python compute error
#input data folder ?
#output data folder ?
#name of dataset ?


#TODO maybe should be more simple
if [ -z "$inputdir" ] || [ -z "$outputdir" ]; then
    echo "input or output dir is unset or empty, please be sure to set environment variable -e inputdir= and -e outputdir correctly"
    #exit 1;#To decomment if we want lock instead of warning only
else
    if [ -d "$inputdir" ] && [ -d "$outputdir" ]; then 
        if [ -w "$outputdir" ] && [ -r "$inputdir" ]; then #Test to be sure have right to write to output folder
            echo "Configuration is valid, start processing... "
        else
            echo "input dir is not readable or output is not writable , please be sure to mount volume with write permissions"
            #exit 1;#To decomment if we want lock instead of warning only
        fi
    else
        echo "input or output dir doesn't exists, please be sure to mount volume"
        #exit 1;#To decomment if we want lock instead of warning only
    fi
fi

#TODO to replace with python call with environnment variable
result=$(python3 <<EOF
import os
for key in os.environ.keys():
    print ("%30s=%s \\n" % (key,os.environ[key]))
EOF
)

echo $result
exit 0;



