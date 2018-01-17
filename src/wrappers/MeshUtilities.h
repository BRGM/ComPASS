//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <iosfwd>

#include "COC.h"

struct MeshConnectivity
{
	COC NodebyCell;
	COC NodebyFace;
	COC FacebyNode;
	COC FacebyCell;
	COC CellbyNode;
	COC CellbyFace;
	COC CellbyCell;
};

std::ostream& operator<<(std::ostream& os, const MeshConnectivity& connectivity)
{

	// Classical nested loops
	os << "Mesh cell nodes:" << std::endl;
	for (auto && cell : connectivity.NodebyCell) {
		for (auto && node : cell) {
			os << node << " ";
		}
		os << std::endl;
	}
	// Alternative way of dumping nodes using STL facilities
	os << "Mesh cell faces:" << std::endl;
	for (auto && face : connectivity.NodebyFace) {
		std::copy(std::begin(face), std::end(face), std::ostream_iterator<int>(os, " "));
		os << std::endl;
	}
	return os;

}
