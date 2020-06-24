#include <Python.h>

int main(int argc, char *argv[]) {
   Py_Initialize();
   PyRun_SimpleString(
       "from os import makedirs\n"
       "makedirs('foodir/foosubdir')\n");
   Py_Finalize();
}
