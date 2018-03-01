#pragma once

#include <vector>
#include <limits>

#include <boost/optional.hpp>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

namespace ExperimentalBRGM
{

	template <typename BaseShape, typename Kernel = typename CGAL::Kernel_traits<BaseShape>::Kernel>
	class DTM 
	{
	public:
		typedef BaseShape Base_shape;
		typedef typename Kernel::Point_3 Point;
		typedef std::vector<Base_shape> Patches;
	protected:
		typedef typename Kernel::Vector_3 Vector;
		typedef typename Kernel::Line_3 Line;
		typedef typename Patches::const_iterator ShapeIterator;
		typedef typename CGAL::AABB_triangle_primitive<Kernel, ShapeIterator> Primitive;
		typedef typename CGAL::AABB_traits<Kernel, Primitive> PT; // Primitive Traits
		typedef CGAL::AABB_tree<PT> Tree;
		typedef typename Tree::template Intersection_and_primitive_id<Line>::Type Line_intersection_type;
		typedef boost::optional< Line_intersection_type > Line_intersection;
		typedef typename Kernel::FT FT;
		Patches patches;
		Tree tree;
		void rebuild_tree()
		{
			tree.rebuild(patches.begin(), patches.end());
			tree.accelerate_distance_queries();
		}
	public:
		DTM(DTM&& other) noexcept :
		patches{ std::move(other.patches) },
			tree{}
		{
			rebuild_tree();
		}
		DTM(const Patches& all_shapes):
		patches{ all_shapes },
			tree{}
		{
			rebuild_tree();
		}
		DTM(Patches&& all_shapes) noexcept:
		patches{ std::move(all_shapes) },
			tree{}
		{
			rebuild_tree();
		}
		FT depth(const Point& P) const
		{
			const Vector up(0, 0, 1);
			Line_intersection topo_intersection = tree.any_intersection(Line(P, up));
			assert(topo_intersection);
			if(topo_intersection) {
				if(const Point *I = boost::get<Point>(&(topo_intersection->first))) {
					return I->z() - P.z();
				}
			}
			assert(false); // should not reach this point
			return std::numeric_limits<FT>::max(); // may fails for some FT
		}

		/* return true if MNT is above point P, i.e. O is below surface */
		bool is_above(const Point& P) const { return depth(P)>0; }
		auto number_of_shapes() const { return patches.size(); }
		const auto& shapes() const { return shapes; }
	};

	template <typename Kernel, typename PointIterator>
	auto build_dtm_from_triangulation(PointIterator first, PointIterator last)
	{
		typedef typename CGAL::Projection_traits_xy_3<Kernel> Gt; // Geometric Traits
		typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
		auto triangulation = Delaunay{ first, last };
		typedef typename Kernel::Triangle_3 Triangle;
		auto triangles = typename DTM<Triangle>::Patches{};
		triangles.reserve(triangulation.number_of_faces());
		for (auto p = triangulation.finite_faces_begin(); p != triangulation.finite_faces_end(); ++p) {
			triangles.emplace_back(Triangle(p->vertex(0)->point(), p->vertex(1)->point(), p->vertex(2)->point()));
		}
		triangles.shrink_to_fit();
		auto result = DTM<Triangle>{ triangles };
		return result;
	}


} // namespace ExperimentalBRGM


