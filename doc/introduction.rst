About ComPASS
=============

.. include:: README.rst

Installation instructions
=========================

Following `The Zen of Python <https://www.python.org/dev/peps/pep-0020/>`_:

    There should be one -- and preferably only one -- obvious way to do it

Unfortunately subtle differences may exist from one system to another.
Please do not hesitate to submit issues/comments to improve this section.

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

You will also need a few additional python modules:

.. code-block:: shell

  python3 -m pip install scikit-build setuptools-scm sortedcontainers verstr vtkwriters numba


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

    python3 setup.py develop --build-type Debug -j 4 -DComPASS_WITH_water2ph_PHYSICS=ON

will compile in `Debug` mode with 4 compilation threads and will activate the *water2ph* physics,

  .. code-block:: shell

    python3 setup.py install -DComPASS_WITH_ALL_PHYSICS=ON

will compile and install all available physics.

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


Install in a conda environment
------------------------------

This relies on the conda-forge channel and will create an isolated `compass` environment.

Check the comment in the following bash script to adpat to your own needs.
Packages in the conda-forge channel are sometimes a bit too new...
You might get instabilities if you change version numbers.
The script has been successfully tested at the beginning of February 2022.


Linux script
^^^^^^^^^^^^

.. literalinclude:: ../miscellaneous/install-linux-with-conda
   :language: bash
   :linenos:

Download script:
:download:`install-on-linux-with-conda <../miscellaneous/install-linux-with-conda>`

MacOS script (with clang)
^^^^^^^^^^^^^^^^^^^^^^^^^

This script is not fully tested yet (feedback is wellcome).

.. literalinclude:: ../miscellaneous/install-mac-with-conda
   :language: bash
   :linenos:

Download script:
:download:`install-on-linux-with-conda <../miscellaneous/install-mac-with-conda>`


Using docker environments
-------------------------

.. include:: using_docker.rst
