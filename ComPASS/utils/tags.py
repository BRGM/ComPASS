import numpy as np

from .. import mpi


class EdgeTagger:
    """We use bitwise or to code different edge family flagging both nodes.
    """

    @staticmethod
    def family_id(i):
        # No more than 30 families because we encode family on a 4 bits integer
        assert i < 31, "too many families"
        return 2 ** i

    def __init__(self, edges_families, verbose=False):
        """:param edges_families: edges to be tagged with edges described by their nodes
           (i.e. edges_families is a sequence of sequence of tuple of size 2 or a sequence of 2 columns arrays) 
           :param verbose: if True will output a bit of information during the tagging operation 
        """
        assert len(edges_families) < 31, "no more than 30 families"
        self.edges_families = edges_families
        self.verbose = verbose

    def __call__(self, simulation):
        flags = simulation.global_nodeflags()
        assert np.all(flags == 0), "all nodes flags must be set to 0 to tag edges"
        for i, edges_nodes in enumerate(self.edges_families):
            if self.verbose:
                print("Family", i, "has", np.sum(edges_nodes), "nodes")
            flags[edges_nodes] = np.bitwise_or(
                flags[edges_nodes], EdgeTagger.family_id(i)
            )


def tag_edges_families(simulation, families, verbose=False):
    """Builds an EdgeTagger to tag all edges families.
       Parameters are the same as :py:class:`EdgeTagger`
    """
    EdgeTagger(families, verbose)(simulation)


def fracture_edges_and_node_flags(simulation, verbose=False):
    fracture_edges = simulation.find_fracture_edges(simulation.frac_face_id())
    if verbose:
        print(fracture_edges.shape[0], "fracture edges on proc", mpi.proc_rank)
    flags = simulation.nodeflags()[fracture_edges]
    return fracture_edges, flags


def retrieve_fracture_edges_families(simulation, verbose=False):
    """Retrieve all fracture edges families among fracture edges.
       :param verbose: if True will output a summary of what has been retrieved
       :return: a dicitonnary of edges families with the key being the family id
    """
    fracture_edges, flags = fracture_edges_and_node_flags(simulation, verbose)
    result = {}
    for fi in range(31):
        mask = np.bitwise_and(flags, EdgeTagger.family_id(fi))
        fi_edges = fracture_edges[np.nonzero(mask[:, 0] & mask[:, 1])]
        if len(fi_edges) > 0:
            result[fi] = fi_edges
            if verbose:
                print(
                    f"Family {fi} has {fi_edges.shape[0]} edges"
                    f" with {np.unique(fi_edges).shape[0]} nodes"
                    f" on proc {mpi.proc_rank}"
                )


def retrieve_fracture_edges_with_node_tag(simulation, node_tags, verbose=False):
    """Retrieve all fracture edges sharing a node_atgs.
       :param node_tags: a node tag or a list of node_tags
       :param verbose: if True will output a summary of what has been retrieved
       :return: a list of (possibly empty) list of edges in the same order as the node_tags parameter
    """
    try:
        node_tags = list(node_tags)
    except TypeError:
        node_tags = [int(node_tags)]
    fracture_edges, flags = fracture_edges_and_node_flags(simulation, verbose)
    result = [
        fracture_edges[(flags[:, 0] == tag) & (flags[:, 1] == tag)]
        for tag in node_tags
    ]
    if verbose:
        for tag, edges in zip(node_tags, result):
            print(
                f"{edges.shape[0]} edges found for node tag {tag}"
                f" with {np.unique(edges).shape[0]} nodes"
                f" on proc {mpi.proc_rank}"
            )
    if len(result) == 1:
        return result[0]
    return result
