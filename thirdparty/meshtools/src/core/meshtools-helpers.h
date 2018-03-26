#pragma once

#include <iostream>

#include "meshtools.h"

namespace MT = MeshTools;

template <typename T, std::size_t n>
std::ostream& operator<<(std::ostream& os, const std::array<T, n>& a)
{
	os << "(";
	for (std::size_t k = 0; k < n; ++k) {
		os << a[k];
		if (k < n - 1) os << ",";
	}
	os << ")";
	return os;
}

	template <typename Coordinate, std::size_t dimension>
std::ostream& operator<<(std::ostream& os, const MT::Point<Coordinate, dimension>& P) {
	os << "Point" << static_cast<const std::array<Coordinate,dimension>&>(P);
	return os;
}

template <std::size_t n>
std::ostream& operator<<(std::ostream& os, const MT::Element_by_nodes<n>& element) {
	os << "Element<" << n << ">" << element.nodes;
	return os;
}

std::ostream& operator<<(::std::ostream& os, const MT::FaceNeighbors& twins)
{
	os << twins.value().first;
	if (twins.is_inside()) {
		os << " " << twins.value().second;
	}
	return os;
}

template <typename Element_type>
std::ostream& dump_info(std::ostream& os)
{
	os << Element_type::name << ":"
		<< " " << Element_type::nbnodes() << " nodes"
		<< " " << Element_type::nbfacets() << " facets"
		<< " (VTK id " << Element_type::VTK_ID << ")"
		<< std::endl;
	return os;
}

template <typename Element_type>
std::ostream&  dump_nodes(const Element_type& element, std::ostream& os)
{
	os << Element_type::name << element.nodes;
	return os;
}

