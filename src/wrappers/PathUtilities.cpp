//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both the terms
// of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
// and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <iostream>
#include <string>

#include <Python.h>

extern "C"
{
  // *** WARNING *** WARNING *** WARNING *** WARNING ***
  // WARNING: The following assumes that a python interpreter is running.
  //          Otherwise encapsulation between Py_Initialize() / Py_Finalize() is needed.
  void c_make_directory(const char * s)
  {
    //std::cerr << "%% creating directory: *|" << s << "|*" << std::endl;
    std::string python_code;
    python_code.append("from os import makedirs\n");
    //python_code.append("print('%% python attempt at: ").append(s).append("')\n");
    python_code.append("try:\n");
    python_code.append("  makedirs('").append(s).append("')\n");
    //python_code.append("  print('%% python created: ").append(s).append("')\n");
    python_code.append("except FileExistsError:\n");
    // CHECKME: Silently ignores if target directory already exists this may be necessary due to concurrent calls to the routine
    python_code.append("  pass\n");
    PyRun_SimpleString(python_code.c_str());
  }
}
