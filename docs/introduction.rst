About ComPASS
=============

.. include:: README.rst

Installation instructions
=========================

Install on `Ubuntu 20.04 focal <https://releases.ubuntu.com/20.04/>`_
---------------------------------------------------------------------

All python packages can obviously be installed using your favorite package manager.

Packages used to build and run simulations:

.. code-block:: shell

  sudo apt-get install --yes python3-dev python3-click python3-numpy python3-setuptools \
    python3-pytest-xdist python3-pip python3-wheel \
    wget build-essential gcc gfortran cmake libopenmpi-dev \
    libmetis-dev libpetsc-real-dev python3-mpi4py cmake-curses-gui

Manually instal petsc4py the version must match PETSc version (3.12 for Ubuntu 20.04)

.. code-block:: shell

  wget https://bitbucket.org/petsc/petsc4py/downloads/petsc4py-3.12.0.tar.gz
  tar xf petsc4py-3.12.0.tar.gz
  cd petsc4py-3.12.0
  sudo python3 setup.py install
  python3 -m pip install sortedcontainers


Add the following definitions at the end your ```.bashrc`` file.
``mpicc`` is used to find PETSC through linker flags.

.. code-block:: shell

  export CC=mpicc
  export PETSC_DIR=/usr/lib/petsc

Then clone `ComPASS repository <https://gitlab.inria.fr/charms/ComPASS>`_:

  - through ssh:

  .. code-block:: shell

    git clone git@gitlab.inria.fr:charms/ComPASS.git

  - or https:

  .. code-block:: shell

    git clone https://gitlab.inria.fr/charms/ComPASS.git

Then cd to the root directory and run the installation:

  .. code-block:: shell

    pip3 install .

or alternatively in develop mode with verbose output:

  .. code-block:: shell

    pip3 -vvv install -e .

Other systems
-------------

.. include:: INSTALL.rst
