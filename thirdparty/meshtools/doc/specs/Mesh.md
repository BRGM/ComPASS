# Mesh concepts

## Mesh elements


Mesh must define element types (*a priori* implemented by enums)

Each element type must define a compile time constant topology (`constexpr`):
* `nb_nodes` is the total number of nodes of the elements
* `nb_facets` is the total number of facets of the elements
* `facet_type` is the type of facet (which is defined if and only if `nb_facets` is non zero


* `dim` is the spatial dimension of the element 

## Connectivity

The basic access is :
ConnectivityTable which is a `CoC<NodeId>`
Vertex = Node ?
Cell Facet Edge Node (abstract topology)
Node id => vertices table

### Structured meshes

Special coc with iterator with node valarray

Full fledged mesh structure:
* nb nodes + id_node_min, id_node_max
* edges
  * nodes
* facets
  * edges
  * nodes
* cells
  * facets
  * edges
  * nodes

Lightweight mesh structure:
* cells
  * nodes

Adjacencies:
*edges
  * neighboring_facets
* facets
  * neighboring_cells
  * adjacent_cells (**fluxes**)
* cells
  * incident_facets
  * adjacent_cells
