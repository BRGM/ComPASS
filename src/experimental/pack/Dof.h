#pragma once

namespace ComPASS
{

    //// CHECKME: Is Undefined usefull?
    //enum class Dof_label : std::size_t {
    //    Undefined = 0, 
    //    Node, 
    //    Cell, 
    //    Fracture 
    //};

    //template <Dof_label dof_label, typename T>
    //struct Dof {
    //    static constexpr Dof_label label = dof_label;
    //    typedef typename T type;
    //};

    template <typename T> 
    struct Dof {
        static constexpr bool is_ComPASS_dof = true;
        typedef T type;
    };

    struct Node : Dof<double> {};
    struct Cell : Dof<double> {};
    struct Fracture : Dof<double> {};

} // namespace ComPASS