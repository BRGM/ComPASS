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

		struct Flowrate_operating_conditions
		{
			double flowrate;
			double pressure_limit;
		};

		struct Pressure_operating_conditions
		{
			double pressure;
			double flowrate_limit;
		};

		typedef boost::variant<Flowrate_operating_conditions, Pressure_operating_conditions> Operating_conditions;

		struct Stopped_well_control
		{};

		struct Injection_well_control
		{
			Operating_conditions operating_conditions;
			// FIXME: Add PhysicalState structures, including temperature
			double temperature;
		};

		struct Production_well_control
		{
			Operating_conditions operating_conditions;
		};

		typedef boost::variant<Stopped_well_control, Injection_well_control, Production_well_control> Well_control;

		struct Is_stopped_well : public boost::static_visitor<bool> {
			template <typename T>
			result_type operator()(const T&) const { return false; }
			result_type operator()(const Stopped_well_control&) const { return true; }
		};

		struct Is_injection_well : public boost::static_visitor<bool> {
			template <typename T>
			result_type operator()(const T&) const { return false; }
			result_type operator()(const Injection_well_control&) const { return true; }
		};

		struct Is_production_well : public boost::static_visitor<bool> {
			template <typename T>
			result_type operator()(const T&) const { return false; }
			result_type operator()(const Production_well_control&) const { return true; }
		};

		struct Is_pressure_operating : public boost::static_visitor<bool> {
			template <typename T>
			result_type operator()(const T&) const { return false; }
			result_type operator()(const Injection_well_control& control) const {
				return boost::apply_visitor(Is_pressure_operating(), control.operating_conditions);
			}
			result_type operator()(const Production_well_control& control) const {
				return boost::apply_visitor(Is_pressure_operating(), control.operating_conditions);
			}
			result_type operator()(const Pressure_operating_conditions&) const { return true; }
		};

		struct Is_flowrate_operating : public boost::static_visitor<bool> {
			template <typename T>
			result_type operator()(const T&) const { return false; }
			result_type operator()(const Injection_well_control& control) const {
				return boost::apply_visitor(Is_flowrate_operating(), control.operating_conditions);
			}
			result_type operator()(const Production_well_control& control) const {
				return boost::apply_visitor(Is_flowrate_operating(), control.operating_conditions);
			}
			result_type operator()(const Flowrate_operating_conditions&) const { return true; }
		};

		struct Well {
			Well_geometry geometry;
			Well_control control;
			bool is_stopped() const { return boost::apply_visitor(Is_stopped_well(), control); }
			bool is_injecting() const { return boost::apply_visitor(Is_injection_well(), control); }
			bool is_producing() const { return boost::apply_visitor(Is_production_well(), control); }
			bool is_flowrate_operating() const { return boost::apply_visitor(Is_flowrate_operating(), control); }
			bool is_pressure_operating() const { return boost::apply_visitor(Is_pressure_operating(), control); }
		};

	} // end of namespace Well

} // end of namespace ComPASS
