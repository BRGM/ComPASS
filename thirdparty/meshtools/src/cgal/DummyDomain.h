#pragma once

#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace ExperimentalOracle
{

	/** Implements the MeshDomain_3 concept cgal-4.4/doc_html/Mesh_3/classMeshDomain__3.html
	    but without anything done. */
	template <typename Kernel, typename SubdomainIndex, typename SurfacePatchIndex>
	class DummyDomain
	{
	public:
		// typedef from MeshDomain_3 concept
		typedef Kernel R;
		typedef typename Kernel::Point_3 Point_3;
		typedef typename Kernel::Segment_3 Segment_3;
		typedef typename Kernel::Ray_3 Ray_3;
		typedef typename Kernel::Line_3 Line_3;
		typedef CGAL::Tag_false Has_features;
		typedef SubdomainIndex Subdomain_index;
		typedef SurfacePatchIndex Surface_patch_index;
		typedef boost::variant<SubdomainIndex, SurfacePatchIndex> Index;
		typedef CGAL::cpp11::tuple<Point_3, Index, int> Intersection;

		struct Construct_initial_points {
			Construct_initial_points(const DummyDomain& domain) {}
			template<typename OutputIterator>
			OutputIterator operator()(OutputIterator out, int n = 10) const {
				assert(false);
				std::cerr << "DUMMY DOMAIN: Should not be called." << std::endl;
				return out;
			}
		};

		struct Is_in_domain {
			Is_in_domain(const DummyDomain& domain) {}
			boost::optional<Subdomain_index> operator()(const Point_3& P) const {
				assert(false);
				std::cerr << "DUMMY DOMAIN: Should not be called." << std::endl;
				return boost::optional<Subdomain_index>();
			}
		};

		struct Construct_intersection {
			Construct_intersection(const DummyDomain& domain) {}
			template <typename Query>
			Intersection operator()(const Query& query) const {
				assert(false);
				std::cerr << "DUMMY DOMAIN: Should not be called." << std::endl;
				return Intersection();
			}
		};

		// The following functions give access to the function objects:
		Construct_initial_points construct_initial_points_object() const { return Construct_initial_points(); }
		Is_in_domain is_in_domain_object() const { return Is_in_domain(); }
		Construct_intersection construct_intersection_object() const { return Construct_intersection(); }

		// Converters between index types
		Index index_from_subdomain_index(Subdomain_index subdomain_index) const { return Index(subdomain_index); }
		Index index_from_surface_patch_index(Surface_patch_index surface_patch_index) const { return Index(surface_patch_index); }
		Subdomain_index subdomain_index(Index index) const { return boost::get<Subdomain_index>(index); }
		Surface_patch_index surface_patch_index(Index index) const { return boost::get<Surface_patch_index>(index); }

	};


} // namespace ExperimentalOracle
