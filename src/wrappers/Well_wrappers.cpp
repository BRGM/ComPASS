//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>

#include "ArrayWrapper.h"
#include "COC.h"
#include "Well.h"

constexpr int NC = ComPASS_NUMBER_OF_COMPONENTS;
constexpr int NP = ComPASS_NUMBER_OF_PHASES;

namespace py = pybind11;

/** Common data structure shared by injectors and producers. */
struct Fortran_well_data {
   typedef std::array<double, NC> Component_vector;
   int id;
   double radius;
   double maximum_pressure;
   double minimum_pressure;
   double imposed_flowrate;
   Component_vector injection_composition;
   double injection_temperature;
   double actual_mass_flowrate;
   double actual_energy_flowrate;
   double actual_pressure;
   double actual_temperature;
   char operating_code;
};

/** Common data structure shared by injectors and producers. */
struct Perforation_state  // WellPerforationState_type in IncCVWell.F90
{
   double pressure;
   double temperature;
   double density;
   std::array<double, NP> saturation;
   double pressure_drop;
};

struct Perforation_data  // Fortran TYPE_DataNodeWell
{
   int parent_vertex_id;  // num of parent; -1 if head node
   int parent_offset;  // pt of parent; -1 if head node ! FIXME: improve doc !!!
   int parent_rank;    // pt of parent; -1 if head node ! FIXME: improve doc !!!
   double well_index_Darcy;
   double well_index_Fourier;
};

struct Well_perforations {
   Perforation_state* perforations_begin;
   std::size_t nb_perforations;
   auto size() const { return nb_perforations; }
};

struct PerforationDataCSRWrapper {
   const int nb_wells;            // Fortran TYPE_CSRDataNodeWell%Nb
   const int* well_offset;        // Fortran TYPE_CSRDataNodeWell%Pt
   const int* node_vertex;        // Fortran TYPE_CSRDataNodeWell%Num
   const Perforation_data* data;  // Fortran TYPE_CSRDataNodeWell%Val
   PerforationDataCSRWrapper()
       : nb_wells{0},
         well_offset{nullptr},
         node_vertex{nullptr},
         data{nullptr} {}
};

struct Well_information {
   std::size_t nb_perforations;
   const int* perforations_vertex;
   const Perforation_state* perforations_state;
   const Perforation_data* perforations_data;
};

// Fortran functions
#include "MeshSchema.fh"
extern "C" {
void Well_allocate_well_geometries(COC&, COC&);
void Well_set_wells_data(ArrayWrapper&, ArrayWrapper&);
Fortran_well_data* get_injectors_data();
std::size_t nb_injectors();
Fortran_well_data* get_producers_data();
std::size_t nb_producers();
Well_perforations get_producing_perforations(std::size_t);
Well_perforations get_injecting_perforations(std::size_t);
void retrieve_CSRData_producer(PerforationDataCSRWrapper&);
void retrieve_CSRData_injector(PerforationDataCSRWrapper&);
}

struct Well_group_information {
   std::vector<Well_information> info;
   ComPASS::Well::Well_type well_type;
};

auto unpack_well_information(const ComPASS::Well::Well_type well_type) {
   PerforationDataCSRWrapper wrapper;
   switch (well_type) {
      case ComPASS::Well::Well_type::injector:
         retrieve_CSRData_injector(wrapper);
         break;
      case ComPASS::Well::Well_type::producer:
         retrieve_CSRData_producer(wrapper);
         break;
      default:
         assert(false);  // Should never been reached
   }
   const auto nb_wells = static_cast<std::size_t>(wrapper.nb_wells);  // FIXME!
   std::vector<Well_information> info;
   info.reserve(nb_wells);
   Well_perforations perforations;
   for (std::size_t wk = 0; wk < nb_wells; ++wk) {
      switch (well_type) {
         case ComPASS::Well::Well_type::injector:
            perforations = get_injecting_perforations(wk);
            break;
         case ComPASS::Well::Well_type::producer:
            perforations = get_producing_perforations(wk);
            break;
         default:
            assert(false);  // Should never been reached
      }
      const auto nb_perforations = perforations.nb_perforations;
      assert(nb_perforations ==
             *(wrapper.well_offset + wk + 1) - *(wrapper.well_offset + wk));
      info.emplace_back(Well_information{
          nb_perforations, wrapper.node_vertex + *(wrapper.well_offset + wk),
          perforations.perforations_begin,
          wrapper.data + *(wrapper.well_offset + wk)});
   }
   return Well_group_information{info, well_type};
}

struct Well_perforations_iterator {
   ComPASS::Well::Well_type well_type;
   std::size_t well_id;
   Well_perforations operator*() {
      switch (well_type) {
         case ComPASS::Well::Well_type::producer:
            return get_producing_perforations(well_id);
            break;
         case ComPASS::Well::Well_type::injector:
            return get_injecting_perforations(well_id);
            break;
         default:
            assert(false);
      }
      return Well_perforations{nullptr, 0};
   }
   Well_perforations_iterator& operator++() {
      ++well_id;
      return *this;
   }
   bool operator<(const Well_perforations_iterator& other) const {
      assert(well_type == other.well_type);
      return well_id < other.well_id;
   }
   bool operator==(const Well_perforations_iterator& other) const {
      assert(well_type == other.well_type);
      return well_id == other.well_id;
   }
   bool operator!=(const Well_perforations_iterator& other) const {
      assert(well_type == other.well_type);
      return well_id != other.well_id;
   }
};

using namespace ComPASS::Well;

std::ostream& operator<<(std::ostream& os, const Well& well) {
   auto& geometry = well.geometry;
   os << "WELL" << std::endl;
   os << "    radius: " << geometry.radius << std::endl;
   os << "    segments: ";
   for (auto&& S : geometry.segments) {
      os << S[0] << "->" << S[1] << ", ";
   }
   return os;
}

struct WellGeometryInfo {
   std::vector<int> offsets;
   std::vector<int> edges;
   void shrink_to_fit() {
      offsets.shrink_to_fit();
      edges.shrink_to_fit();
   }
};

auto extract_geometry_info(py::list wells,
                           std::function<bool(const Well&)> select) {
   auto geometry_info = WellGeometryInfo{};
   auto& offsets = geometry_info.offsets;
   auto& edges = geometry_info.edges;
   assert(edges.empty());
   offsets.push_back(0);
   for (auto&& p : wells) {
      auto well = p.cast<const Well&>();
      if (select(well)) {
         for (auto&& edge : well.geometry.segments) {
            edges.push_back(edge[0]);
            edges.push_back(edge[1]);
         }
         offsets.push_back(edges.size());
      }
   }
   geometry_info.shrink_to_fit();
   return geometry_info;
}

void set_well_geometries(py::list wells) {
   auto producers_geometries = extract_geometry_info(
       wells, [](const Well& well) { return well.is_producing(); });
   auto injectors_geometries = extract_geometry_info(
       wells, [](const Well& well) { return well.is_injecting(); });
   auto producers =
       COC::wrap(producers_geometries.offsets, producers_geometries.edges);
   auto injectors =
       COC::wrap(injectors_geometries.offsets, injectors_geometries.edges);
   Well_allocate_well_geometries(producers, injectors);
}

struct Producer_data_info {
   int id;
   double radius;
   double minimum_pressure;
   double imposed_flowrate;
   char operating_code;  // 'p' for pressure mode; 'f' for flowrate mode; 'c'
                         // for closed
   Producer_data_info(const Well& well) {
      assert(well.is_producing());
      id = well.id;
      radius = well.geometry.radius;
      if (well.operates_on_pressure()) {
         operating_code = 'p';
         assert(well.control.operating_conditions
                    .is<Pressure_operating_conditions>());
         auto operating_conditions = well.control.operating_conditions
                                         .get<Pressure_operating_conditions>();
         minimum_pressure = operating_conditions.pressure;
         imposed_flowrate = operating_conditions.flowrate_limit;
      } else {
         if (well.operates_on_flowrate()) {
            operating_code = 'f';
            assert(well.control.operating_conditions
                       .is<Flowrate_operating_conditions>());
            auto operating_conditions =
                well.control.operating_conditions
                    .get<Flowrate_operating_conditions>();
            minimum_pressure = operating_conditions.pressure_limit;
            imposed_flowrate = operating_conditions.flowrate;
         } else {
            assert(well.is_closed());
            operating_code = 'c';
            minimum_pressure =
                std::numeric_limits<decltype(minimum_pressure)>::min();
            imposed_flowrate = 0;
         }
      }
      // std::cerr << "Create Producer_data_info: " << radius << " " <<
      // minimum_pressure << " " << imposed_flowrate << std::endl;
   }
};

struct Injector_data_info {
   int id;
   double radius;
   double temperature;
   double maximum_pressure;
   double imposed_flowrate;
   char operating_code;  // 'p' for pressure mode; 'f' for flowrate mode; 'c'
                         // for closed
   Injector_data_info(const Well& well) {
      assert(well.is_injecting());
      id = well.id;
      radius = well.geometry.radius;
      assert(well.control.status.is<Injection_well_status>());
      auto status = well.control.status.get<Injection_well_status>();
      temperature = status.temperature;
      if (well.operates_on_pressure()) {
         operating_code = 'p';
         assert(well.control.operating_conditions
                    .is<Pressure_operating_conditions>());
         auto operating_conditions = well.control.operating_conditions
                                         .get<Pressure_operating_conditions>();
         maximum_pressure = operating_conditions.pressure;
         imposed_flowrate = operating_conditions.flowrate_limit;
      } else {
         if (well.operates_on_flowrate()) {
            operating_code = 'f';
            assert(well.control.operating_conditions
                       .is<Flowrate_operating_conditions>());
            auto operating_conditions =
                well.control.operating_conditions
                    .get<Flowrate_operating_conditions>();
            maximum_pressure = operating_conditions.pressure_limit;
            imposed_flowrate = operating_conditions.flowrate;
         } else {
            assert(well.is_closed());
            operating_code = 'c';
            maximum_pressure =
                std::numeric_limits<decltype(maximum_pressure)>::max();
            imposed_flowrate = 0;
         }
      }
      if (imposed_flowrate > 0) {
         imposed_flowrate = -imposed_flowrate;
         std::cerr << "WARNING - setting negative injection flowrate"
                   << std::endl;
      }
      // std::cerr << "Create Injector_data_info: " << radius << " " <<
      // temperature << " " << maximum_pressure << " " << imposed_flowrate <<
      // std::endl;
   }
};

/** Change well ids so that each each has a unique id. */
void check_well_ids(py::list& wells) {
   std::set<std::size_t> ids;
   std::size_t new_id = 0;
   bool renumbered_wells = false;
   auto compute_new_id = [&]() {
      while (ids.count(new_id) == 1) {
         ++new_id;
      }
      return new_id;
   };
   for (auto& p : wells) {
      Well& well = p.cast<Well&>();
      if (ids.count(well.id) == 1) {
         well.id = compute_new_id();
         renumbered_wells = true;
      }
      ids.insert(well.id);
   }
   if (renumbered_wells) {
      py::print("WARNING - Some wells were renumbered to have a unique id.");
   }
}

void output_well_ids(py::list& wells) {
   std::ostringstream oss;
   oss << "Well ids:"
       << " ";
   for (auto&& w : wells) {
      oss << " " << w.cast<Well&>().id;
   }
   py::print(oss.str());
}

void set_well_data(py::list& wells, const bool display_well_ids) {
   check_well_ids(wells);
   if (display_well_ids) output_well_ids(wells);
   std::vector<Producer_data_info> producers_info;
   std::vector<Injector_data_info> injectors_info;
   for (auto& p : wells) {
      const Well& well = p.cast<const Well&>();
      if (well.is_producing()) producers_info.emplace_back(well);
      if (well.is_injecting()) injectors_info.emplace_back(well);
   }
   producers_info.shrink_to_fit();
   injectors_info.shrink_to_fit();
   auto wrapped_producers_info = ArrayWrapper::wrap(producers_info);
   auto wrapped_injectors_info = ArrayWrapper::wrap(injectors_info);
   Well_set_wells_data(wrapped_producers_info, wrapped_injectors_info);
}

void add_well_wrappers(py::module& module) {
   py::class_<Well_geometry>(module, "WellGeometry")
       .def(py::init<>())
       .def_readwrite("radius", &Well_geometry::radius)
       .def("add_segments",
            [](Well_geometry& self, py::iterable& segments) {
               auto to_segment = [](py::sequence& s) {
                  return Well_geometry::Segment_type{
                      {s[0].cast<Well_geometry::Node_id_type>(),
                       s[1].cast<Well_geometry::Node_id_type>()}};
               };
               for (auto&& S : segments) {
                  auto seq = S.cast<py::sequence>();
                  self.segments.push_back(to_segment(seq));
               }
            })
       .def("clear_segments",
            [](Well_geometry& self) { self.segments.clear(); })
       .def_property_readonly("segments", [](const Well_geometry& self) {
          const auto& segments = self.segments;
          return py::array_t<Well_geometry::Node_id_type, py::array::c_style>{
              {segments.size(), static_cast<std::size_t>(2)},
              reinterpret_cast<const Well_geometry::Node_id_type*>(
                  segments.data())};
       });

   py::class_<Well>(module, "Well")
       .def(py::init<>())
       .def_readwrite("id", &Well::id)
       .def_readwrite("geometry", &Well::geometry)
       .def_property_readonly("closed", &Well::is_closed)
       .def_property_readonly("injecting", &Well::is_injecting)
       .def_property_readonly("producing", &Well::is_producing)
       .def("close",
            [](Well& instance) {
               instance.control.status = Closed_well_status{};
            })
       .def("produce",
            [](Well& instance) {
               instance.control.status = Production_well_status{};
            })
       .def("inject",
            [](Well& instance, double T) {
               instance.control.status = Injection_well_status{T};
            })
       .def_property(
           "operate_on_pressure",
           [](Well& instance) -> py::object {
              py::object result = py::none{};
              if (instance.operates_on_pressure()) {
                 assert(instance.control.operating_conditions
                            .is<Pressure_operating_conditions>());
                 auto operating_conditions =
                     instance.control.operating_conditions
                         .get<Pressure_operating_conditions>();
                 return py::make_tuple(operating_conditions.pressure,
                                       operating_conditions.flowrate_limit);
              }
              return result;
           },
           [](Well& instance, py::tuple args) {
              instance.control.operating_conditions =
                  Pressure_operating_conditions{args[0].cast<double>(),
                                                args[1].cast<double>()};
           })
       .def_property(
           "operate_on_flowrate",
           [](Well& instance) -> py::object {
              py::object result = py::none{};
              if (instance.operates_on_flowrate()) {
                 assert(instance.control.operating_conditions
                            .is<Flowrate_operating_conditions>());
                 auto operating_conditions =
                     instance.control.operating_conditions
                         .get<Flowrate_operating_conditions>();
                 return py::make_tuple(operating_conditions.flowrate,
                                       operating_conditions.pressure_limit);
              }
              return result;
           },
           [](Well& instance, py::tuple args) {
              instance.control.operating_conditions =
                  Flowrate_operating_conditions{args[0].cast<double>(),
                                                args[1].cast<double>()};
           })
       .def("__repr__", [](Well& instance) {
          std::ostringstream os;
          os << instance;
          return os.str();
       });

   module.def("set_well_geometries", &set_well_geometries,
              "Set well geometries.");
   module.def("set_well_data", &set_well_data, "Set well data.",
              py::arg("wells"), py::arg("display_well_ids") = false);

   py::class_<Fortran_well_data>(module, "WellData")
       .def_readonly("id", &Fortran_well_data::id)
       .def_readwrite("operating_code", &Fortran_well_data::operating_code)
       .def_readwrite("radius", &Fortran_well_data::radius)
       .def_readwrite("maximum_pressure", &Fortran_well_data::maximum_pressure)
       .def_readwrite("minimum_pressure", &Fortran_well_data::minimum_pressure)
       .def_readwrite("imposed_flowrate", &Fortran_well_data::imposed_flowrate)
       .def_readwrite("injection_temperature",
                      &Fortran_well_data::injection_temperature)
       .def_readonly("actual_mass_flowrate",
                     &Fortran_well_data::actual_mass_flowrate)
       .def_readonly("actual_energy_flowrate",
                     &Fortran_well_data::actual_energy_flowrate)
       .def_readonly("actual_pressure", &Fortran_well_data::actual_pressure)
       .def_readonly("actual_temperature",
                     &Fortran_well_data::actual_temperature)
       .def_property_readonly("is_closed",
                              [](const Fortran_well_data& self) {
                                 return self.operating_code == 'c';
                              })
       .def("open", [](Fortran_well_data& self) { self.operating_code = 'f'; })
       .def("close",
            [](Fortran_well_data& self) { self.operating_code = 'c'; });

   py::class_<Perforation_state>(module, "PerforationState")
       .def_readonly("pressure", &Perforation_state::pressure)
       .def_readonly("temperature", &Perforation_state::temperature)
       .def_readonly("density", &Perforation_state::density)
       .def_property_readonly("saturation",
                              [](py::object& self) {
                                 auto state = self.cast<Perforation_state*>();
                                 return py::array_t<double, py::array::c_style>{
                                     NP, state->saturation.data(), self};
                              })
       .def_readonly("pressure_drop", &Perforation_state::pressure_drop);

   py::class_<Perforation_data>(module, "PerforationData")
       .def_readonly("parent_vertex", &Perforation_data::parent_vertex_id)
       .def_readonly("parent_offset", &Perforation_data::parent_offset)
       .def_readonly("parent_rank", &Perforation_data::parent_rank)
       .def_readonly("well_index_Darcy", &Perforation_data::well_index_Darcy)
       .def_readonly("well_index_Fourier",
                     &Perforation_data::well_index_Fourier);

   auto pyWellInfo =
       py::class_<Well_information>(module, "WellInformation")
           .def_readonly("nb_perforations", &Well_information::nb_perforations)
           .def_property_readonly(
               "vertices",
               [](const Well_information& self) {
                  // Convert to C rank
                  const auto n = self.nb_perforations;
                  auto result = py::array_t<int, py::array::c_style>{n};
                  const auto p = self.perforations_vertex;
                  std::transform(p, p + n, result.mutable_data(0),
                                 [](auto i) { return i - 1; });
                  return result;
               },
               py::return_value_policy::reference_internal);

   auto add_perforation_state_property = [&pyWellInfo](const char* name,
                                                       std::size_t offset) {
      pyWellInfo.def_property_readonly(
          name,
          [offset](py::object& self) {
             auto info = self.cast<Well_information*>();
             std::vector<std::size_t> shape;
             shape.push_back(info->nb_perforations);
             std::vector<std::size_t> strides;
             strides.push_back(sizeof(Perforation_state));
             auto ptr = reinterpret_cast<const double*>(
                 reinterpret_cast<const unsigned char*>(
                     info->perforations_state) +
                 offset);
             return py::array_t<double, py::array::c_style>{shape, strides, ptr,
                                                            self};
          },
          py::return_value_policy::reference_internal);
   };

   add_perforation_state_property("pressure",
                                  offsetof(Perforation_state, pressure));
   add_perforation_state_property("temperature",
                                  offsetof(Perforation_state, temperature));
   add_perforation_state_property("density",
                                  offsetof(Perforation_state, density));
   add_perforation_state_property("pressure_drop",
                                  offsetof(Perforation_state, pressure_drop));

   pyWellInfo.def_property_readonly(
       "saturation",
       [](py::object& self) {
          auto info = self.cast<Well_information*>();
          std::vector<std::size_t> shape;
          shape.push_back(info->nb_perforations);
          shape.push_back(NP);
          std::vector<std::size_t> strides;
          strides.push_back(sizeof(Perforation_state));
          strides.push_back(sizeof(double));
          auto ptr = reinterpret_cast<const double*>(
              reinterpret_cast<const unsigned char*>(info->perforations_state) +
              offsetof(Perforation_state, saturation));
          return py::array_t<double, py::array::c_style>{shape, strides, ptr,
                                                         self};
       },
       py::return_value_policy::reference_internal);

   auto add_perforation_data_property = [&pyWellInfo](const char* name,
                                                      std::size_t offset) {
      pyWellInfo.def_property_readonly(
          name,
          [offset](py::object& self) {
             auto info = self.cast<Well_information*>();
             auto ptr = reinterpret_cast<const double*>(
                 reinterpret_cast<const unsigned char*>(
                     info->perforations_data) +
                 offset);
             std::vector<std::size_t> shape;
             shape.push_back(info->nb_perforations);
             std::vector<std::size_t> strides;
             strides.push_back(sizeof(Perforation_data));
             return py::array_t<double, py::array::c_style>{shape, strides, ptr,
                                                            self};
          },
          py::return_value_policy::reference_internal);
   };

   add_perforation_data_property("well_index_Darcy",
                                 offsetof(Perforation_data, well_index_Darcy));
   add_perforation_data_property(
       "well_index_Fourier", offsetof(Perforation_data, well_index_Fourier));

   auto add_perforation_rank_property = [&pyWellInfo](const char* name,
                                                      std::size_t offset) {
      pyWellInfo.def_property_readonly(name, [offset](py::object& self) {
         auto info = self.cast<Well_information*>();
         auto ptr = reinterpret_cast<const int*>(
             reinterpret_cast<const unsigned char*>(info->perforations_data) +
             offset);
         const auto nbperfs = info->nb_perforations;
         assert(nbperfs > 1);
         auto result = py::array_t<int, py::array::c_style>{nbperfs - 1};
         auto p = result.mutable_data(0);
         for (std::size_t k = 0; k + 1 < nbperfs; ++k, ++p) {
            auto rank = *reinterpret_cast<const int*>(
                reinterpret_cast<const unsigned char*>(ptr) +
                k * sizeof(Perforation_data));
            assert(rank > 0);  // Fortran index > 0
            *p = rank - 1;
         }
         assert(*reinterpret_cast<const int*>(
                    reinterpret_cast<const unsigned char*>(ptr) +
                    (nbperfs - 1) * sizeof(Perforation_data)) ==
                -1);  // Well head
         return result;
      });
   };

   // parent reservoir node
   add_perforation_rank_property("parent_vertex",
                                 offsetof(Perforation_data, parent_vertex_id));
   // parent offset in *wells CSR*
   add_perforation_rank_property("parent_offset",
                                 offsetof(Perforation_data, parent_offset));
   // parent local node index
   add_perforation_rank_property("parent_rank",
                                 offsetof(Perforation_data, parent_rank));

   py::class_<Well_group_information>(module, "WellGroup")
       .def_property_readonly(
           "nb_wells",
           [](const Well_group_information& self) { return self.info.size(); })
       .def("__iter__", [](const Well_group_information& self) {
          return py::make_iterator(begin(self.info), end(self.info));
       });

   module.def("injectors_information", []() {
      return unpack_well_information(ComPASS::Well::Well_type::injector);
   });

   module.def("producers_information", []() {
      return unpack_well_information(ComPASS::Well::Well_type::producer);
   });

   module.def(
       "injectors_data",
       [](bool own_only) {
          auto p = get_injectors_data();
          const std::size_t n =
              own_only ? number_of_own_injectors() : nb_injectors();
          return py::make_iterator(p, p + n);
       },
       py::return_value_policy::reference, py::arg("own_only") = false);

   module.def(
       "producers_data",
       [](bool own_only) {
          auto p = get_producers_data();
          const std::size_t n =
              own_only ? number_of_own_producers() : nb_producers();
          return py::make_iterator(p, p + n);
       },
       py::return_value_policy::reference, py::arg("own_only") = false);

   py::class_<Well_perforations>(module, "WellPerforations")
       .def("__len__", &Well_perforations::size)
       .def(
           "__iter__",
           [](Well_perforations& self) {
              return py::make_iterator(
                  self.perforations_begin,
                  self.perforations_begin + self.nb_perforations);
           },
           py::return_value_policy::reference)
       .def(
           "__getitem__",
           [](Well_perforations& self, int i) {
              assert(self.nb_perforations < std::numeric_limits<int>::max());
              const auto n = static_cast<int>(self.nb_perforations);
              if (i >= n || i < -n) {
                 py::print("wrong perforation index (", i, "out of",
                           self.nb_perforations, ")");
                 throw pybind11::index_error{};
              }
              if (i < 0) return self.perforations_begin + (n + i);
              return self.perforations_begin + i;
           },
           py::return_value_policy::reference);

   module.def("producers_perforations", []() {
      constexpr auto well_type = ComPASS::Well::Well_type::producer;
      return py::make_iterator(
          Well_perforations_iterator{well_type, 0},
          Well_perforations_iterator{well_type, nb_producers()});
   });

   module.def("injectors_perforations", []() {
      constexpr auto well_type = ComPASS::Well::Well_type::injector;
      return py::make_iterator(
          Well_perforations_iterator{well_type, 0},
          Well_perforations_iterator{well_type, nb_injectors()});
   });

   module.def("nb_producers", &nb_producers);
   module.def("nb_injectors", &nb_injectors);
}
