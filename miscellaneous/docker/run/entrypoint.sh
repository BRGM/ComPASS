#!/bin/bash

export COMPASS_USER=compass
export COMPASS_DEFAULT_UID=1000
mkdir -p /etc/compass
useradd --create-home --shell /bin/bash --uid ${COMPASS_DEFAULT_UID} ${COMPASS_USER}

cat >>/etc/bash.bashrc <<EOL
#
# ComPASS specific
#
# The following is to pass mpi/petsc flags easily to cmake
export CC=mpicc
# The following is due to an output error with openmpi3
# cf. https://github.com/open-mpi/ompi/issues/4948
export OMPI_MCA_btl_vader_single_copy_mechanism=none
EOL

CUSTSHELL=/etc/compass/customize-session
SIMSCRIPT=/etc/compass/simulation-process

# creates a default file in case nothing is to be done
echo "echo" > $SIMSCRIPT

python3 /etc/compass/entrypoint.py \
    --customize-session $CUSTSHELL \
    --simulation-process $SIMSCRIPT "$@"
source $CUSTSHELL

# /localfs should already be present
mkdir -p /localfs
chown ${COMPASS_USER}:${COMPASS_USER} /localfs

echo "You are running a ComPASS container as user" $COMPASS_USER "with UID" `id -u $COMPASS_USER` 

su-exec compass `cat $SIMSCRIPT`

