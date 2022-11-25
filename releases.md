Version 4.5.1
=============

This is a bug fix release that aims at having a robust installation
procedure in a conda environment.
A nightly CI test has been set-up to test a full conda environment
creation from scratch and a simple ComPASS script is run in the
newly created environment.

We kept a single compass.yml conda environment file which is
easier to maintain and less prone to divergence than two separate
files. The documentation has been updated accordingly.

Pybind11 version 2.10.0 is skipped
as the make_iterator API changed introducing a regression.

A small bug has also been corrected in the selection of boundary fracture edges.

Version 4.5.0
=============

This new ComPASS release is the first of the v4.5.x series.
The minor version was increased to reflect the new ability to define
several physical properties directly in python.

Fluid molar and volumetric mass densities, molar enthalpy and dynamic viscosity
can be defined as standard python function. The functions are jit-compiled as
C functions using numba. No performance overhead has been observed (generated
code is even a bit faster). The documentation has been updated to introduce
this new API and CI tests were added or modified. The consistency
of properties' derivatives is also checked automatically using
finite differences.

The ComPASS.physics package was renamed into ComPASS.properties so as to
use *physics* to design what was previously and misleadlingly called *eos*
without being a true equation of state.
The documentation was updated to reflect these changes.

The *immiscible2ph* physics is now part of the default physics and has a CI test.
Some consistency checks were added for initial conditions if initial
context is not *diphasic*.

The simplification molar density = volumetric density for *water2ph* physics
(which has a single component) has been documented. This simplification is
bound to be removed in next releases.

This release includes several code cleaning, refactoring and simplification.
The timeloop analysis has been improved. The help of simulation *fake* methods
(i.e. function expecting a simulation object as their first argument) is now
correctly displayed in python shells (i.e. you can do `help(simulation.init)`).

Installation instructions have been updated relying on the conda environment
files and no longer using anaconda repositories.

This release also fixes several bugs including:
* the initialization of accumulation of absent components (*Ctilde*) that was skipped,
* a slight numerical error due to real kind of a fortran parameter
  in molar enthalpy,
* the fact that the fugacity function is now the true thermodynamic fugacity property
  (and not a fugacity coefficient).

An API change in pybind11 2.10.0 introduced a bug in ComPASS, this problem
was apparently solved with pybind11 2.10.1 (yet to be checked). The
conda environment files were configured to avoid the pybind11 releases
higher than v2.9.2 for now.

Several radioactive waste storage tests were added in the test/cases directory.

A multi-segmented well model coupled with the reservoir model was added.
This is still considered to be an experimental feature.

Release notes are now versioned in the releases.md file.

Version 4.4.1
=============

* Avoid duplicate pipelines
* Refactor pipeline management using workflow
* Default branch CI done only at night
* Documentation Mass balance equation in water2ph
* Documentation Glossary
* Standalone version of multi-segmented wells
* Allow mesh mapping when reloading a snapshot
* Simplify code
* Version specifications for numpy and numba
* Requires python>=3.7
* Clean atmospheric BC property
* BugFix Correct brine properties derivatives
* New ResiduMSWells and JacobianMSWells modules
* Simulation class becomes a true (singleton) class
* Make simulation.init a *fake* member method
* Documentation and comments
* Comments
* Use fugacity instead of fugacity coefficients
* Remove unused use clause
* Correct error message in linear water fugacity
* No timeloop logs with no_output flag set to false
* Improve utilies to analyse yaml log files
* Rename van Genuchten capillary pressure parameter
* Add comments for van Genuchten rel. permeability
* Document local generation of documentation
* Multi-segmented wells core modules
* New module to compute MSWells physical quantities
* Change subroutines interface in LoisThermoHdyro
* New IncPrimSecdMSWells module
* Improve conda recipe
* Refactor subroutine to avoid temporary array
* Add FIXME in flash about ideal gas hypothesis
* Use messages.warning instead of print
* Additional warning in standard_loop
* Set global flags after building connectivities
* New MSWellsData module
* New IncCVMSWells module
* MPI Barrier not useful in VAG volume distribution
* New data sctructures for multi-segmented wells
* Improve comment
* Simplify code
* Improve comments
* Improve subroutine interface
* Bugfix in reset_freeflow_nodes when empty array
* Comments in vanGenuchten pc
* Typo
* Doc Use of conda environments
* Doc Use sphinx.ext.autosectionlabel
* Doc Remove wrong info in tutorial
* Doc Add execution command for conda installation
* Doc Improve installation documentation
* Doc Improve rocktypes and add flags
* BugFix Update conda environments names
* Add conda environment files
* Switch to pyproject.toml and setup.cfg
* Additional cmake file to choose compiled physics
* Add test on interactions Dir and atm nodes
* Reset the FF info when modifying the Dir nodes
* Improve reset FF faces to remove all FF faces
* Dirichlet nodes (in FF face) are not flag as FF node
* No VAG volume at FF nodes when disrtibuting the cell vol
* Improve tests: use tick as argument in callbacks
* Reduce simulation domain size
* New test using atmospheric boundary conditions
* Fix test
* Change integer kind to avoid compilation warning
* Do not allow abitrary keywords arguments
* Use PetscCall for PETSc>3.17
* Fix cmake message
* Add kr in tutorial directory, change the plots
* Improve tutorial and documentation
* Upgrade black version
* Use github repository in documentation
* Small documentation reorganization
* Add a section about MeshTools installation
* Add missing dependencies in documentation
* Improve Newton behavior description
* Add comment
* Store last computed residual scales
* Doc: add an example of restriction condition in set()
* Typo in documentation
* In conda installation build all physics by default
* Starting script of the tutorial
* Tutorial instructions on the documentation
* Use python3 as the default python interpreter
* Typos
* Slightly refactor documentation structure
* Typo
* BugFix Definition of nb_output
* Change atmospheric boundary condition defaults
* Update scripts in cases directory
* Documentation: start of eos section
* Documentation: boundary conditions
* Documentation: Pc and kr with rocktypes
* Documentation: small introduction and reference to examples
* Documentation: time-step
* Documentation : fractures section
* Documentation: how to create or import a mesh
* Typo in documentation
* Add case test with atmospheric boundary conditions
* Update test case
* New CI test
* Change default pressure increment

Version 4.4.0
=============

This new ComPASS release is the first of the v4.4.x series.
We increased the minor version to reflect the changes
in installation that now relies on scikit-build.
Several installation scripts relying on conda
can be found in the miscellaneous directory.
Compilation with clang works though it
has not been thoroughly tested.
This release comes with several improvements:

* the possibility to dump mesh information before running the simulation
new API for the freeflow diphasic physics with
the possibility to have spatial variations for some of
its parameters
* the ability to import faulted Petrel grids that will
be remeshed using CGAL under the hood
* the export of meshes as generic polytope meshes
(fractures are exported as VTK PolyData)
* the Jacobian module has been slightly refactor to
prepare the transition to the v5.x series, more
functions are exposed on the python side
* the indexing of phases is now absolute throughout the code
* a tolerance has been added to create vertical wells
from mesh nodes
* the dependance on Metis is now optional and we can use
another partitioner
* van Genuchten capillary pressure and
relative permeability models have been introduced
* simulation.build_state API has been improved to init simulation
using rocktypes or diphasic equilibriums

This release includes several code cleaning, refactoring and simplification as well as performance enhancements:
* removing some of the temporary arrays in Fortran
* user defined functions such as capillary pressure functions
may rely on numba

We set up Architecture Decision Records in the doc/adr directory.

Version is now generated from git tags using setuptools_scm package.
Version numbers are wrapped in a verstr version object
to provide user friendly version checks.

The documentation has been updated to reflect code and API changes.
Several sections have been added.
The way the documentation is generated has been refactored.
C++-17 is now enforced which removes an external dependency on
the mapbox library.

The support for PETSc versions lower than 3.8 has been dropped.

The dependance on MeshTools is relaxed with the objective
to fully replace MeshTools by several independant packages in the v5.x series.

CI tests are now run in Debug mode which is more robust
and will allow to identify some potential memory bugs
(array bound checking).
New CI tests have been added and an
effort has been made to increase the performance of CI pipelines.
Pipelines can be skipped for branch names starting with NOCI.

This release also includes the following bug fixes:
* alignment matrix for immiscible2ph
* empty freeflow faces are correctly handled
* several bug fixes concern the freeflow diphasic module
(several CI tests have been added)
* PETSc failures are correctly documented
* a bug without incidence has been fixed in local Schur subroutine
