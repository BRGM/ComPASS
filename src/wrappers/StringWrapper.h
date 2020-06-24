//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#pragma once

#include <string>

struct StringWrapper {
   const char* pointer;
   const size_t length;
   StringWrapper() = delete;
   StringWrapper(const StringWrapper&) = delete;
   void operator=(const StringWrapper&) = delete;
   StringWrapper(const std::string& s)
       : pointer(s.c_str()), length(s.length()) {}
};
