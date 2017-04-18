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
    std::string python_code;
    python_code.append("from os import makedirs\n");
    python_code.append("if not os.path.exists('").append(s).append("'):\n");
	python_code.append("    try:\n");
	python_code.append("        makedirs('").append(s).append("')\n");
	python_code.append("    except FileExistsError:\n");
	// CHECKME: Silently ignores if target directory already exists this may be necessary due to concurrent calls to the routine
	python_code.append("        pass\n");
	PyRun_SimpleString(python_code.c_str());
  }
}
