//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "meshtools.h"

namespace MeshTools
{

	// FIXME: This is a bit redundant but gcc asks for both definition and declaration of static types
	// cf. https://stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
	constexpr decltype(Vertex_info::facets) Vertex_info::facets;
	constexpr decltype(Segment_info::facets) Segment_info::facets;
	constexpr decltype(Triangle_info::facets) Triangle_info::facets;
	constexpr decltype(Tetrahedron_info::facets) Tetrahedron_info::facets;
	constexpr decltype(Hexahedron_info::facets) Hexahedron_info::facets;

} // namespace MeshTools


