#include <iostream>

#include "mapbox/optional.hpp"

int main()
{

	auto i = mapbox::util::optional<int>{};

	if (i) {
		std::cout << "Good" << std::endl;
	}
	i = 2;
	if (i) {
		std::cout << "Initialized: " << *i << std::endl;
	}
	return 0;
}