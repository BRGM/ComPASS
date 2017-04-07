#pragma once

#include <iosfwd>

#include "COC.h"

struct Vertices
{
	double * p;
	std::size_t nb_points;
	Vertices() :
		p(NULL),
		nb_points(0)
	{}
};

struct MeshConnectivity
{
	COC NodebyCell;
	COC NodebyFace;
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
