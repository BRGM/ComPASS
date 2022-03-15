//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <array>
#include <iterator>
#include <list>
#include <variant>
#include <vector>

namespace ComPASS {

namespace Well {

enum struct Well_type { injector, producer, undefined };

struct Well_geometry {
   typedef int Node_id_type;
   typedef std::array<Node_id_type, 2> Segment_type;
   double radius;
   std::vector<Segment_type> segments;
};

struct Undefined_operating_conditions {};

struct Flowrate_operating_conditions {
   double flowrate;
   double pressure_limit;
   Flowrate_operating_conditions(double Q, double Plim)
       : flowrate{Q}, pressure_limit{Plim} {}
};

struct Pressure_operating_conditions {
   double pressure;
   double flowrate_limit;
   Pressure_operating_conditions(double P, double Qlim)
       : pressure{P}, flowrate_limit{Qlim} {}
};

using Operating_conditions =
    std::variant<Undefined_operating_conditions, Flowrate_operating_conditions,
                 Pressure_operating_conditions>;

struct Closed_well_status {};

struct Production_well_status {};

struct Injection_well_status {
   // FIXME: Add PhysicalState structures, including temperature
   double temperature;
   Injection_well_status(double T) : temperature{T} {}
};

using Well_status = std::variant<Closed_well_status, Injection_well_status,
                                 Production_well_status>;

struct Well_control {
   Operating_conditions operating_conditions;
   Well_status status;
   template <typename WellStatus>
   constexpr bool check_status() const noexcept {
      return std::holds_alternative<WellStatus>(status);
   }
   template <typename WellOperatingConditions>
   constexpr bool check_operating_conditions() const noexcept {
      return std::holds_alternative<WellOperatingConditions>(
          operating_conditions);
   }
   bool is_closed() const { return check_status<Closed_well_status>(); }
   bool is_injecting() const { return check_status<Injection_well_status>(); }
   bool is_producing() const { return check_status<Production_well_status>(); }
   bool is_flowrate_operating() const {
      return check_operating_conditions<Flowrate_operating_conditions>();
   }
   bool is_pressure_operating() const {
      return check_operating_conditions<Pressure_operating_conditions>();
   }
};

struct Well {
   std::size_t id;
   bool is_multi_segmented = false;
   Well_geometry geometry;
   Well_control control;
   bool is_closed() const { return control.is_closed(); }
   bool is_injecting() const { return control.is_injecting(); }
   bool is_producing() const { return control.is_producing(); }
   bool operates_on_flowrate() const { return control.is_flowrate_operating(); }
   bool operates_on_pressure() const { return control.is_pressure_operating(); }
};

}  // end of namespace Well

}  // end of namespace ComPASS
