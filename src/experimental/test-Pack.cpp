#include "Dof.h"
#include "Pack.h"

#include <iostream>

using namespace ComPASS;

template <typename View>
void dump(const View& v, std::string name = "")
{
    if(!name.empty()) {
        std::cout << "Number of " << name.c_str() << " values: " << v.size() << std::endl;
    }
    for (auto&& x : v) {
        std::cout << x << std::endl;
    }
}

int main(int argc, char * argv[])
{

    Pack<Node, Cell, Fracture> pack{ 5, 2, 3 };

    auto nodes = pack.view<Node>();
    std::cout << "Number of nodes: " << nodes.size() << std::endl;
    std::cout << "Nodes values: " << std::endl;
    for (auto&& v : nodes) {
        std::cout << v << std::endl;
    }
    // you can access directly a view element
    nodes[2] = 3.14;
    dump(nodes, "new nodes");
    // you can make a view of contiguous dof with the same type
    auto reservoir = view<Node, Cell, Fracture>(pack);
    dump(reservoir, "reservoir");
    //auto bad_view = view<Node, Fracture>(pack); // compilation fails not contiguous dofs

    // classical for loop with counter
    for (std::size_t k = 0; k != nodes.size(); ++k) {
        nodes[k] = static_cast<double>(k);
        nodes[k] *= 0.1;
    }
    dump(nodes, "new nodes again");

    return 0;

}
