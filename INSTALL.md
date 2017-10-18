# Log of previous installations

## From "scratch"

Scripts to download third party libraries and compile them are to be found in miscellaneous/install.


## Install on Centos 7:

https://gitlab.inria.fr/charms/ComPASS/wikis/Compass-install-on-centos-7

##Install on Mac OS X:

https://gitlab.inria.fr/charms/ComPASS/wikis/compass-install-on-macos-x

# General considerations

Cmake modules to find external packages :
* FindPETSC and others from Jed Brown's github repository https://github.com/jedbrown/cmake-modules
* FindMETIS from Dune setup

Usefull packages on ubuntu 16.10:
build-essential
(cmake) cmake-curses-gui or cmake-qt-gui
gfortran
liblapack-dev
libmetis-dev
mpi-default-dev
libvtk6-dev


## PetsC - PetsC version must NOT be greater than 3.5

We currently rely on ILU-0 from eculid package (Hypre) that has been removed from version 3.6 (cf. [v3.6 changes](https://www.mcs.anl.gov/petsc/documentation/changes/36.html)).
To download PETSC version <= 3.5.4 go [there](http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.5.tar.gz)

A possible configuration can be :
./configure PETSC_ARCH=linux-gnu --prefix=$PETSC_INSTALL_PREFIX --download-fblaslapack --with-cc=mpicc --with-cxx=mpicxx -with-fc=mpif90 --with-fortran=1 --with-clanguage=C --download-hypre

Where the PETSC_INSTALL_PREFIX stands for where you want to install petsC.

Useful environment variables (for CMake):

  * PETSC_DIR - directory in which PETSc resides
  * PETSC_ARCH - build architecture

FindPETSC is going to look in ${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h
PETSC_DIR can aslo be set to the compilation location


## HDF5

Useful environment variables (for CMake):

  * HDF5_ROOT

A few HDF5 calls rely on Fortran 2003  : h5offsetof [is available only in Fortran 2003 environments](https://support.hdfgroup.org/HDF5/doc/fortran/FortranFlags.html) (cf. CMAKE Flag HDF5_ENABLE_2003 if compiling HDF5 from source).

