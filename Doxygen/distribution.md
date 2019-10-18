Distribution    {#distributionpage} 
============

# METIS distribution

The ComPASS code is a MPI parallel code. 
The parallelism is performed by distributing a part of 
the global domain to each processor. 
We call "global" everything which concern the total 
domain (before the decomposition) and "local" 
what concerns only one processor in particular.
After reading or computing the global mesh, METIS 
is called to spread the cells over the
processors, it determine the set
of own cells for each processor.
The vector ProcbyCell of size NbCell is filled 
with the numero of the processor for which the cell 
is own. 

\image html global_local_distribution.png  
\image latex global_local_distribution.png "Illustration of the decomposition of the global domain." width=10cm

<br>

# Ghost cells

The processor is responsible for the computation 
in its own cells. To do so it needs 
to have access to one more layer of cells to
have the necessary data in its local memory. 
This additional elements are called ghost elements.
"Local" means every elements which are stored 
in one processor, then it includes the "own" plus 
the "ghost" elements.


\image html own_ghost_distribution2.png   
\image latex own_ghost_distribution2.png width=10cm


<br>


# What about the Nodes and Faces ?

METIS distributes only the cells, but nodes 
and faces also need to been spread. Indeed, 
once each cell is distributed, 
at the interface between two CPUs
one has to chose for which of the two processors 
the element will be own.
And for load balancing reason it is not a good 
choice that one of the two CPUs owns all of 
the interface nodes/faces.
This is why the determination is done in a 
random way. The ComPASS code use the following
discriminator : 
 - the node is own for the processor which has among 
its own cells the last cell of the list CellbyNode.
 - the face is own for the processor which has among 
its own cells the last cell of the list CellbyFace.



# Particularity of the wells

There is one distinction about the wells since
it is necessary to have the informations about 
all the nodes of the wells to be able to 
rebuild the informations at one well node. 
It means that every processor
concerned by a well (the well get through its local 
domain) will have every nodes of this well 
in his own or ghost elements. For the sake of simplicity
of the code, if one well node is own or ghost for the
processor, we suppose that the local domain is concerned
by the well. It means that even if a well cross only the ghost
domain it will be stored, which is not necessary for the 
computation but it was easier to implement, and it should
not slow down the performances.

Thus the presence of well does not change the own nodes of a 
processor but it can add ghost nodes.

## Own well

The well is own for the procesor to which the head node
is own. It implies that the concerned CPUs will be in charge
of the additionnal unknown Pw due to the well. 

# Code architecture
 
The connectivity of the global domain is done in 
the file GlobalMesh.F90. Then the determination of 
the own and ghost elements are performed in LocalMesh.F90 
<b>only by the master processor</b>. This informations are
finaly send in the file MeshSchema.F90 so every CPU gets its
local informations.




