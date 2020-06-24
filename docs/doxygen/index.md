Main Page    {#mainpage}
=========

The objectives are to develop efficient parallel numerical methods
to simulate gas liquid compositional thermal flow in high energy
geothermal reservoirs. The spatial discretization is adapted to
unstructured polyhedral meshes and takes into account discrete
fracture or fault networks represented as 2D surfaces connected
to the surrounding 3D heterogeneous anisotropic porous medium (the matrix).
The difficulties result from the highly contrasted spatial and temporal scales
and petrophysical properties between the 3D matrix and the fault network.
The second difficulty is to account for the strong couplings and high
nonlinearities induced by the large range of pressure and temperature
and by the phase appearance and disapearance.  <br>

This project is a collaboration started in 2014 between:  <br>
<b>BRGM</B>: Simon Lopez and Farid Smai    <br>
<b>Joint INRIA-LJAD team Coffee</b>: Feng Xing, Laurence Beaude,
Konstantin Brenner and Roland Masson   <br>
An ANR project ?????????? has been submitted to strenghen the project
and enlarge it to other partners like Storengy, La Maison de la Simulation et
le Laboratoire Jacques Louis Lions.



<!--- One sentence on what is the PROJECT. -->
<!--- One sentence on who is the targeted audience. -->


## Documentation and Information

<a target="_blank" href="http://compass.gforge.inria.fr/">A website</a> is dedicated to
the code ComPASS with in particular illustrations of simulations run with ComPASS.      <br>

It is possible to download the code ComPASS from <a target="_blank"
href="https://gforge.inria.fr/">INRIA forge</a>.

## How to Compile the code ComPASS

To compile the code ComPASS you will need to fulfil the minimum requirements and compile
it through several stages.   ???

### Requirements

To be compiled the following tools are required :

        * cmake (version minimum required 2.8)
        * C and Fortran Compilers
         - gcc (minimum version 4.9)
         - Intel (minimum version 15.0)
        * MPI
         - OpenMPI (minimum version 1.8)
         - IntelMPI (minimum version ???)
        * METIS
        * PETSC (version 3.5.4 maximum)
        * LAPACK
        * VTK (version minimum 6.3)
        * HDF5 (version 1.8.16 recommended)


### Compilation

 - First, create a building directory, for example $My_ComPASS_DIRECTORY/build
  and generate native makefile using cmake:

\verbatim
mkdir $My_ComPASS_DIRECTORY/build
cd $My_ComPASS_DIRECTORY/build
cmake ..
\endverbatim

 - Once the Makefile has been generated, you can compile. It creates the executable NN:

\verbatim
make
\endverbatim

## How to execute the code ComPASS

\verbatim
mpirun -np $num_procs ./NN meshfile.msh report.txt outputvisu
\endverbatim



example of mesh file ????



# Current limitations
