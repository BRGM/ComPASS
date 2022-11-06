About ComPASS
=============

.. include:: README.rst

Installation instructions
=========================

Following `The Zen of Python <https://www.python.org/dev/peps/pep-0020/>`_:

    There should be one -- and preferably only one -- obvious way to do it

Unfortunately subtle differences may exist from one system to another.
Please do not hesitate to submit issues/comments to improve this section.

.. _using conda environments:

Using conda environments
------------------------

As of today, using
`conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html>`_
is probably the fastest way to get ready to use ComPASS.
It will download a lot of dependencies: these may represent a subsequent payload
and will occupy a certain amount of space on your drive. Yet, you will
end up with an isolated environment to use and/or develop ComPASS.
Though ComPASS will still be compiled on your system,
if you are interested in maximising performance and/or exploiting
libraries that have been fine-tuned for your system (e.g. PETSc)
you may want to consider a native build (cf. the next section).

Prerequisites
^^^^^^^^^^^^^

Install conda
"""""""""""""

You will need to have
`conda <https://docs.conda.io/projects/conda/en/latest/index.html>`_
on your system.
Please consider reading the
`installation instructions
<https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.
The easiest (and lightest) way is probably to install
`Miniconda <https://docs.conda.io/projects/conda/en/latest/glossary.html#miniconda-glossary>`_.
Miniconda installers can be found on this
`download page <https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links>`_.

Once you have installed conda (either through Miniconda or Anaconda) you should have
the conda command available on your path. You may consider updating conda to have the
latest release.

.. code-block:: bash

  conda update -y conda


Replace conda with mamba (optional)
"""""""""""""""""""""""""""""""""""

`Mamba <https://github.com/mamba-org/mamba>`_
is an implementation of the conda package manager in C++
which is much faster than conda and will speed-up
many of the installation steps.

You can install it from `conda-forge <https://conda-forge.org/>`_.

.. code-block:: bash

  conda activate base
  conda install -c conda-forge mamba

or for one-liners:

.. code-block:: bash

  conda install -n base -c conda-forge mamba

Then, in the following you can replace all the occurences of
the :code:`conda` command with the :code:`mamba` command.


You're almost there
^^^^^^^^^^^^^^^^^^^

We provide `several conda environments <https://github.com/BRGM/ComPASS/tree/main/conda>`_ for ComPASS that only differ
in their final step.

In both following cases you will need to have conda (or mamba) on your path
and this may require to activate your conda :code:`base` environment
(:code:`conda activate base` or :code:`source activate base` depending
on your settings).

Case 1: you only want to use the latest version of ComPASS
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Just create the environment:

.. code-block:: bash

  conda env create -f https://raw.githubusercontent.com/BRGM/ComPASS/main/conda/compass-latest.yml

The whole process may be a bit long because MeshTools and ComPASS packages will be
downloaded and compiled under the hood, but that's it. Once it is finished,
you can activate the :code:`compass-latest` environment.

.. code-block:: bash

  conda activate compass-latest

And you can start :ref:`using simulation scripts <Setting-up a simulation>`.

.. note::
  Instead of *latest* you may also use the version tag *v4.4.1* that will
  skip the latest developments.

Case 2: you want to develop ComPASS
"""""""""""""""""""""""""""""""""""

First clone the ComPASS repository:

.. code-block:: bash

  git clone https://github.com/BRGM/ComPASS.git

Then go to ComPASS directory and create the *compass* conda environment
that will contain all ComPASS dependencies
(including the download, compilation and installation of the MeshTools package).

.. code-block:: bash

  cd ComPASS
  conda env create -f conda/compass.yml

Once it is finished, you can activate the :code:`compass` environment
and install ComPASS in development mode.

.. code-block:: bash

  conda activate compass
  pip install -e ComPASS

Then you can start :ref:`using simulation scripts <Setting-up a simulation>`.
Any modification in the ComPASS source python scripts will be reflected
immediately in the ComPASS package
(thanks to the pip development mode - cf. :code:`-e` install option).
If you modify Fortran or C++ source files you will have to re-run pip
that will trigger the compilation of modified source files.

Check dedicated sections of this documentation for more information
 (:ref:`developer's corner`).

.. note::
  *MeshTools* and *ComPASS* are developed on the
  `Inria gitlab server <https://gitlab.inria.fr/charms>`_.
  The main branches are mirrored to `github.com/BRGM` so that
  any github URL above can be replaced with `gitlab.inria.fr/charms`
  to have access to development branches.
  An access to `Inria gitlab server <https://gitlab.inria.fr/charms>`_
  can be provided upon `request <mailto:compass@brgm.fr>`_.


Useful `conda env create` options
""""""""""""""""""""""""""""""""""

:code:`--force`: force creation of environment (removing a previously existing environment of the same name)

:code:`-n new_name`: will rename the generated environment

.. _install with Ubuntu:

Native installation on Linux
----------------------------

Requirements
^^^^^^^^^^^^

All python packages can obviously be installed using your favorite package manager.

Example on `Ubuntu 20.04 focal <https://releases.ubuntu.com/20.04/>`_
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Packages used to build and run simulations:

.. code-block:: shell

  sudo apt-get install --yes python3-dev python3-click python3-numpy python3-setuptools \
    python3-pytest-xdist python3-pip python3-wheel \
    wget build-essential gcc gfortran cmake libopenmpi-dev \
    libmetis-dev libpetsc-real-dev python3-mpi4py python3-petsc4py-real cmake-curses-gui git

You may want to install the `python-is-python3 <https://packages.ubuntu.com/focal/python-is-python3>`_
package so that `python3` becomes the default python interpreter.

Add the following definitions at the end your ``.bashrc`` file.
``mpicc`` is used to find PETSC through linker flags.

.. code-block:: shell

  export CC=mpicc
  export PETSC_DIR=/usr/lib/petsc

You will also need a few additional python modules
(beware of the single quote to escape version constraints):

.. code-block:: shell

  python3 -m pip install scikit-build 'setuptools>=61' 'setuptools-scm>=6.2' sortedcontainers verstr vtkwriters 'numpy>=1.21' 'numba!=0.55.2'


MeshTools installation
^^^^^^^^^^^^^^^^^^^^^^

ComPASS v4 series rely on the MeshTools package
that provides various structures and utility functions to handle meshes.
The package is written in C++ and has python bindings
produced using `pybind11 <https://pybind11.readthedocs.io/en/stable>`_.

This dependency will disappear with the v5 series.

You can clone `MeshTools repository <https://github.com/BRGM/MeshTools>`_:

.. code-block:: shell

  git clone https://github.com/BRGM/MeshTools.git

As for ComPASS, the compilation relies on
`scikit-build <https://scikit-build.readthedocs.io/en/latest/index.html>`_
to run `CMake <https://cmake.org/>`_ through `setup.py`.
The following should do the job:

.. code-block:: shell

  cd MeshTools
  python3 setup.py install -DMESHTOOLS_TRIES_TO_USE_CGAL=OFF
  cd ..

The `-DMESHTOOLS_TRIES_TO_USE_CGAL=OFF` flag is optional. The default is that
MeshTools tries to wrap (a small) part of the CGAL library. Some (advanced)
example scripts use CGAL but now rely on the
`pyCGAL package <https://gitlab.brgm.fr/brgm/geomodelling/public/pycgal>`_
that provides bindings closer to the C++ CGAL code.


ComPASS installation
^^^^^^^^^^^^^^^^^^^^

Clone `ComPASS repository <https://github.com/BRGM/ComPASS>`_:

.. code-block:: shell

  git clone https://github.com/BRGM/ComPASS.git

Compilation relies on `scikit-build <https://scikit-build.readthedocs.io/en/latest/index.html>`_
to run `CMake <https://cmake.org/>`_ through `setup.py`.

cd to the root directory and run the installation (by default no physics is install so you will
get a warning, later is explained how to compile one or all physics):

  .. code-block:: shell

    cd ComPASS
    python3 setup.py install

or alternatively in develop mode with verbose output:

  .. code-block:: shell

    cd ComPASS
    python3 setup.py develop

and you may also use all `scikit-build <https://scikit-build.readthedocs.io/en/latest/index.html>`_
options to control the behavior of `CMake <https://cmake.org/>`_.

For example:

  .. code-block:: shell

    python3 setup.py develop --build-type Debug -j 4

will compile in `Debug` mode with 4 compilation threads.

.. note::
  *MeshTools* and *ComPASS* are developed on the
  `Inria gitlab server <https://gitlab.inria.fr/charms>`_.
  The main branches are mirrored to `github.com/BRGM` so that
  any github URL above can be replaced with `gitlab.inria.fr/charms`
  to have access to development branches.
  An access to `Inria gitlab server <https://gitlab.inria.fr/charms>`_
  can be provided upon `request <mailto:compass@brgm.fr>`_.


Troubleshooting
^^^^^^^^^^^^^^^

**PETSc**

The installation of PETSc is undoubtedly one of the trickiest part.
If you need to build PETSc from source please refer to
`PETSc installation instructions <https://www.mcs.anl.gov/petsc/documentation/installation.html>`_.

You will probably want to install `mpi4py` and `petsc4py` along with PETSc, using the
`--download-petsc4py=yes --download-mpi4py=yes --with-mpi4py=yes --with-petsc4py=yes` configuration flags.

It's likely that you need to define `PETSC_ARCH` environment variable along with `PETSC_DIR`.

Try to follow carefully hints given at each stage of the compilation / installation steps of PETSc.

If you need/want to manually instal petsc4py the version must match PETSc version (3.12 for Ubuntu 20.04)

.. code-block:: shell

  wget https://bitbucket.org/petsc/petsc4py/downloads/petsc4py-3.12.0.tar.gz
  tar xf petsc4py-3.12.0.tar.gz
  cd petsc4py-3.12.0
  sudo python3 setup.py install

**Permission related problems**

When running in a docker environment or on a distant machine you may face
permission related problems. Try to add the  setuptools :code:`--user` option at the
end of the compilation directive to use your local python site.


Using distributed docker environments (to be deprecated)
--------------------------------------------------------

If you want or need to use docker, you'd rather consider using
one of the two steps above in an appropriate docker container.

.. include:: using_docker.rst
