//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <array>
#include <vector>
#include <list>
#include <iterator>

// FIXME: use C++17 variant
// #include <variant>
#include <mapbox/variant.hpp>


namespace ComPASS
{

	namespace Well
	{

		struct Well_geometry
		{
			typedef int Node_id_type;
			typedef std::array<Node_id_type, 2> Segment_type;
			double radius;
			std::vector<Segment_type> segments;
		};

		struct Undefined_operating_conditions
		{};

		struct Flowrate_operating_conditions
		{
			double flowrate;
			double pressure_limit;
			Flowrate_operating_conditions(double Q, double Plim) :
				flowrate{ Q },
				pressure_limit{ Plim } {}
		};

		struct Pressure_operating_conditions
		{
			double pressure;
			double flowrate_limit;
			Pressure_operating_conditions(double P, double Qlim) :
				pressure{ P },
				flowrate_limit{ Qlim } {}
		};

		typedef mapbox::util::variant<Undefined_operating_conditions, Flowrate_operating_conditions, Pressure_operating_conditions> Operating_conditions;

		struct Stopped_well_status
		{};

		struct Production_well_status
		{};

		struct Injection_well_status
		{
			// FIXME: Add PhysicalState structures, including temperature
			double temperature;
			Injection_well_status(double T) :
				temperature{ T } {}
		};

		typedef mapbox::util::variant<Stopped_well_status, Injection_well_status, Production_well_status> Well_status;

		struct Well_control
		{
			Operating_conditions operating_conditions;
			Well_status status;
			template <typename WellStatus>
			bool check_status() const { return status.is<WellStatus>(); }
			template <typename WellOperatingConditions>
			bool check_operating_conditions() const { return operating_conditions.is<WellOperatingConditions>(); }
			bool is_stopped() const { return check_status<Stopped_well_status>(); }
			bool is_injecting() const { return check_status<Injection_well_status>(); }
			bool is_producing() const { return check_status<Production_well_status>(); }
			bool is_flowrate_operating() const { return check_operating_conditions<Flowrate_operating_conditions>(); }
			bool is_pressure_operating() const { return check_operating_conditions<Pressure_operating_conditions>(); }
		};

		struct Well {
			Well_geometry geometry;
			Well_control control;
			bool is_stopped() const { return control.is_stopped(); }
			bool is_injecting() const { return control.is_injecting(); }
			bool is_producing() const { return control.is_producing(); }
			bool operates_on_flowrate() const { return control.is_flowrate_operating(); }
			bool operates_on_pressure() const { return control.is_pressure_operating(); }
		};

	} // end of namespace Well

} // end of namespace ComPASS
