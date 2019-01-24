//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

struct CTVector
{
#ifdef _THERMIQUE_
    static constexpr std::size_t npv = ComPASS_NUMBER_OF_COMPONENTS + 1;
#else // _THERMIQUE_
    static constexpr std::size_t npv = ComPASS_NUMBER_OF_COMPONENTS;
#endif // _THERMIQUE_
    double values[npv];
};
