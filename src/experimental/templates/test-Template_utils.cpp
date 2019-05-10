#include <iostream>

#include "Template_utils.h"

int main(int argc, char * argv[])
{

    constexpr bool a = Is_set<int, double, char>::value;
    std::cout << "constexpr " << a << std::endl;
    constexpr bool b = Is_set<int, double, char, int>::value;
    std::cout << "constexpr " << b << std::endl;
    std::cout << "range 0 length " << Range<0>::list::length << std::endl;
    std::cout << "range 4 length " << Range<4>::list::length << std::endl;
    // std::cout << "rank " << Find<int, std::tuple<int, double, int>>::value << std::endl; // compilation fails <int, double, int> is not a set
    //std::cout << "rank " << Find<float, std::tuple<int, double, char>>::value << std::endl; // compilation fails float not in <int, double, int>
    std::cout << "rank " << Find<int, std::tuple<int, double, char>>::value << std::endl;
    std::cout << "rank " << Find<double, std::tuple<int, double, char>>::value << std::endl;
    std::cout << "rank " << Find<char, std::tuple<int, double, char>>::value << std::endl;

#if __cplusplus>=201703L // C++ 17 and over
    std::cout << "constexpr " << is_set<int, double, char> << std::endl;
    std::cout << "constexpr " << is_set<int, double, char, int> << std::endl;
#endif // C++ 17 and over

    return 0;
}
