#include <string>

#include "meshtools.h"
#include "meshtools-helpers.h"

namespace MT = MeshTools;

typedef MT::Point<double, 3> Point;
template <typename ... CellTypes>
using Mesh = MT::Mesh<Point, CellTypes...>;

int main()
{

	dump_info<MT::Tetrahedron>(std::cout);
	dump_info<MT::Wedge>(std::cout);
	dump_info<MT::Hexahedron>(std::cout);

	typedef MT::Tetrahedron Tet;
	auto mesh = Mesh<Tet>{};
	auto& connectivity = mesh.connectivity;
	auto& cells = connectivity.cells;
	cells.nodes.push_back(Tet{ 0, 1, 2, 3 });
	cells.nodes.push_back(Tet{ 4, 1, 2, 3 });
	connectivity.update_from_cellnodes();

	//auto selection = std::array<std::size_t, 2>{ { 0, 2 } };
	//auto selected = MT::extract(cells.nodes.back(), selection);
	//std::cout << "Selection:" << std::endl;
	//std::copy(selected.begin(), selected.end(), std::ostream_iterator<MT::NodeId>(std::cout, " "));
	//for (auto&& i : selected) {
	//	std::cout << i << " ";
	//}
	//std::cout << std::endl;
	
	for (Point::Coordinate_type i = 0; i < 10 * 3; i += 3) {
		mesh.vertices.emplace_back(Point{ i, i + 1, i + 2 });
	}
	auto cell_centers = mesh.cell_centers();
	std::cout << "Cell centers: " << std::endl;
	for (auto&& P : cell_centers)
		std::cout << P << std::endl;
	auto face_centers = mesh.face_centers();
	std::cout << "Face centers: " << std::endl;
	for (auto&& P : face_centers)
		std::cout << P << std::endl;

	const auto& faces = connectivity.faces.nodes;
	auto boundaries = connectivity.boundary_faces();
	std::cout << "Face nodes:" << std::endl;
	for (auto&& face : faces)
		dump_nodes(face, std::cout) << std::endl;
	std::cout << "Boundaries faces:" << std::endl;
	std::copy(boundaries.begin(), boundaries.end(), std::ostream_iterator<MT::FaceId>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "Face cells:" << std::endl;
	for (auto&& twins : connectivity.faces.cells)
		std::cout << twins << std::endl;
	std::cout << "Cell faces:" << std::endl;
	for (auto&& faces : connectivity.cells.faces)
		std::cout << faces << std::endl;

	return 0;

}

