#pragma once

#include <array>
#include <vector>
#include <list>
#include <iterator>

#include <boost/variant.hpp>
// FIXME: use C++17 variant
// #include <variant>


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

		typedef boost::variant<Undefined_operating_conditions, Flowrate_operating_conditions, Pressure_operating_conditions> Operating_conditions;

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

		typedef boost::variant<Stopped_well_status, Injection_well_status, Production_well_status> Well_status;

		template <typename HeldType>
		struct Check_variant : public boost::static_visitor<bool> {
			template <typename T>
			result_type operator()(const T&) const { return false; }
			result_type operator()(const HeldType&) const { return true; }
		};

		struct Well_control
		{
			Operating_conditions operating_conditions;
			Well_status status;
			template <typename WellStatus>
			bool check_status() const { return boost::apply_visitor(Check_variant<WellStatus>(), status); }
			template <typename WellOperatingConditions>
			bool check_operating_conditions() const { return boost::apply_visitor(Check_variant<WellOperatingConditions>(), operating_conditions); }
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
