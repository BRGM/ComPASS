//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <mapbox/variant.hpp>

#include "Tuple_utils.h"

// cf. http://en.cppreference.com/w/cpp/utility/variant/monostate
//struct monostate { };
//constexpr bool operator<(monostate, monostate) noexcept { return false; }
//constexpr bool operator>(monostate, monostate) noexcept { return false; }
//constexpr bool operator<=(monostate, monostate) noexcept { return true; }
//constexpr bool operator>=(monostate, monostate) noexcept { return true; }
//constexpr bool operator==(monostate, monostate) noexcept { return true; }
//constexpr bool operator!=(monostate, monostate) noexcept { return false; }

template <typename TupleType>
struct VariantFromTuple;

template <typename ... Ts>
struct VariantFromTuple<std::tuple<Ts...>>
{
	typedef mapbox::util::variant<Ts...> type;
};

template <typename TupleType>
using VariantFromTupleSet = typename VariantFromTuple<typename TupleSet<TupleType>::type>::type;

template <typename TupleType, bool>
struct MakeVariantFromTuple;

template <typename TupleType>
struct MakeVariantFromTuple<TupleType, true>
{
	static_assert(!is_single_type<TupleType>(), "Tuple is single type, don't use variant!");
	typedef VariantFromTupleSet<TupleType> type;
};

template <typename TupleType>
struct MakeVariantFromTuple<TupleType, false>
{
	static_assert(is_single_type<TupleType>(), "Tuple is not single type!");
	typedef typename SplitTuple<TupleType>::First_element_type type;
};

template <typename TupleType>
using VariantFromTupleIfNeeded = typename MakeVariantFromTuple<TupleType, !is_single_type<TupleType>()>::type;

template <typename ... Ts>
using VariantIfNeeded = VariantFromTupleIfNeeded<std::tuple<Ts...>>;
