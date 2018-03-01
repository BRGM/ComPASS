#pragma once

#include <cassert>
#include <memory>
#include <vector>
#include <map>

#include <CGAL/Kernel_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

namespace ExperimentalBRGM
{

	template <typename Complex>
	class DistanceTrees {
	public:
		typedef typename Complex::Surface_patch_index Surface_patch_index;
		typedef typename Complex::Triangulation Triangulation;
		typedef typename Triangulation::Bare_point Point;
		typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
		typedef typename Kernel::Triangle_3 Triangle;
	protected:
		typedef typename Kernel::FT FT;
		typedef typename std::vector<Triangle>::const_iterator TriangleIterator;
		typedef CGAL::AABB_triangle_primitive<Kernel, TriangleIterator> Primitive;
		typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
		typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
		std::vector<Triangle> triangles;
		Tree bigtree;
		std::map<Surface_patch_index, std::shared_ptr<Tree> > patchtrees;
	public:
		inline DistanceTrees(const Complex& complex);
		/** Return the squared distance between point P and any complex facet. */
		inline FT squared_distance(const Point& P) const;
		/** Return the squared distance between point P and the complex facets of type spi. */
		template <typename SurfacePatchIndex>
		inline FT squared_distance(const Point& P, const SurfacePatchIndex& spi) const;
		template <typename OutputIterator>
		inline OutputIterator patch_indexes(OutputIterator out) const;
	};

	template <typename Complex>
	DistanceTrees<Complex>::DistanceTrees(const Complex& complex):
		bigtree(),
		patchtrees()
	{
		const Triangulation& triangulation = complex.triangulation();
		typedef std::list<Triangle> TriangleList;
		typedef std::map<Surface_patch_index, TriangleList> TriangleMap;
		TriangleMap patches;
		for(typename Complex::Facets_in_complex_iterator pf=complex.facets_in_complex_begin();
			pf!=complex.facets_in_complex_end(); ++pf)
		{
			const Surface_patch_index spi = complex.surface_patch_index(*pf);
			typename TriangleMap::iterator ptl = patches.find(spi);
			if(ptl==patches.end()) {
				TriangleList L;
				L.push_back(triangulation.triangle(*pf));
				patches[spi] = L; 
			} else {
				ptl->second.push_back(triangulation.triangle(*pf));
			}
		}
		size_t nbtriangles = 0;
		for(typename TriangleMap::const_iterator ptl=patches.begin(); ptl!=patches.end(); ++ptl) {
			nbtriangles+= ptl->second.size();
		}
		triangles.reserve(nbtriangles);
		for(typename TriangleMap::const_iterator ptl=patches.begin(); ptl!=patches.end(); ++ptl) {
			// WARNING all this is supposed to work because we have reserved the good amount of triangles
			//         and no reallocation of memory is supposed to occur
			TriangleIterator p = triangles.end();
			std::copy(ptl->second.begin(), ptl->second.end(), std::back_inserter(triangles));
			assert(std::distance(p, (TriangleIterator)triangles.end())==ptl->second.size());
			auto tree = std::make_shared<Tree>();
			tree->rebuild(p, (TriangleIterator) triangles.end());
			bool acceleration_ok = tree->accelerate_distance_queries();
			assert(acceleration_ok);
			patchtrees[ptl->first] = tree;
		}
		bigtree.rebuild(triangles.begin(), triangles.end());
		bool big_acceleration_ok = bigtree.accelerate_distance_queries();
		assert(big_acceleration_ok);
	}

		template <typename Complex>
	template <typename OutputIterator>
	OutputIterator DistanceTrees<Complex>::patch_indexes(OutputIterator out) const
	{

		for(typename std::map<Surface_patch_index, std::shared_ptr<Tree> >::const_iterator p=patchtrees.begin(); p!=patchtrees.end(); ++p) {
			*out = p->first;
			++out;
		}
		return out;
	
	}

		template <typename Complex>
	typename DistanceTrees<Complex>::FT DistanceTrees<Complex>::squared_distance(const Point& P) const
	{

		return bigtree.squared_distance(P);

	}

	template <typename Complex>
	template <typename SurfacePatchIndex>
	typename DistanceTrees<Complex>::FT DistanceTrees<Complex>::squared_distance(const Point& P, const SurfacePatchIndex& spi) const
	{

		assert(patchtrees.find(spi)!=patchtrees.end());
		return patchtrees.at(spi)->squared_distance(P);

	}

} // namespace ExperimentalBRGM
