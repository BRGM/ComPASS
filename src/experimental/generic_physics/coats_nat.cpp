
/** Coats variable */

// Objectif: dissocier physique et choix des variables au niveau de
// l'implémentation

struct CoatsVariable {
   double value;
   operator const double&() const { return value; }
   operator double() const { return value; }
};

struct X {
   // std::variant<p, T, ...> ?
   std::array<double, n>;
};

// using Xnat = std::tuple<> //si tout C++ ??
enum struct variable_id {
   p = 0,
   T = 1,
};

// dans un premier temps au plus proche de l'existant

struct Xnat {  // Fortran compliant

   const std::size_t np;
   std::array<variable_id, nmax> ids;  // Constexpr by context
   double operator[](X& x, const std::size_t& i) { return x[ids[i]]; };
};
'

    // Xsol just reordering of Xnat... so primary / secondary

    struct Physical_property {};

// Physical property
struct rho : Physical_property {
   constexpr std::array<variable_id, n> depends_on;
   // double operator(const X&);
   double operator(const Xphase&) const;
};
