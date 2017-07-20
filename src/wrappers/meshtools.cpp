#include "meshtools.h"

namespace MeshTools
{

	// FIXME: This is a bit redundant but gcc asks for both definition and declaration of static types
	// cf. https://stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
	constexpr decltype(Vertex_info::facets) Vertex_info::facets;
	constexpr decltype(Segment_info::facets) Segment_info::facets;
	constexpr decltype(Triangle_info::facets) Triangle_info::facets;
	constexpr decltype(Tetrahedron_info::facets) Tetrahedron_info::facets;

} // namespace MeshTools


