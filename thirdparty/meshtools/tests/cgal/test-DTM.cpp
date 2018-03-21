#include <string>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "DTM.h"

namespace BRGM = ExperimentalBRGM;
typedef CGAL::Epick Kernel;

int main(int argc, const char* argv[])
{

	if(argc>1) {
		const std::string dtmfile(argv[1]);
		std::ifstream in;
		in.open(dtmfile);
		std::string buffer;
		std::getline(in, buffer); // skip header
		typedef typename Kernel::Point_3 Point;
		std::istream_iterator<Point> begin(in);
		std::istream_iterator<Point> end;
		auto dtm = BRGM::build_dtm_from_triangulation<CGAL::Epick>(begin, end);
		in.close();
		if(argc>2) {
			const std::string pointsfile(argv[2]);
			in.open(pointsfile);
			std::istream_iterator<Point> begin(in);
			std::istream_iterator<Point> end;
			for(std::istream_iterator<Point> p=begin; p!=end; ++p) {
				std::cout << dtm.depth(*p) << std::endl;
			}
		} else {
			std::cout << "Header:" << buffer << std::endl;
			std::cout << dtm.number_of_shapes() << " triangles built" << std::endl;
		}
	}

	return 0;

}
