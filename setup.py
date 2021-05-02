# This is adapted from https://github.com/molpopgen/cmake_example/tree/subdir_and_python

import os
import re
import sys
from pathlib import Path
import platform
import tempfile
import subprocess
from datetime import datetime

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import setuptools_git_versioning as sgv
from distutils.version import LooseVersion
from multiprocessing import cpu_count

config_cmake_path = Path("cmake/config.cmake")
package_name = "ComPASS"
version_file = "version_info"

with Path(f"./{package_name}/{version_file}").open("w") as f:
    print(
        f"{package_name} {sgv.version_from_git()} (branch: {sgv.get_branch()})", file=f,
    )
    print(f"built {datetime.now().isoformat()}", file=f)
    print(f"with python {platform.python_version()} on {platform.node()}", file=f)
    print(f"build system info: {platform.platform()}", file=f)


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        # depends option does not seem to work
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):

        # check that cmake is present on the path and recent enough
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def _retrieve_cmake_build_type(self):
        # We first check if a cmake configuration file is present
        if config_cmake_path.exists():
            with tempfile.NamedTemporaryFile(mode="w+") as cmake_script:
                with config_cmake_path.open() as config_cmake:
                    for line in config_cmake:
                        l = line.strip()
                        if len(l) > 0 and not l.startswith("#"):
                            print(l, file=cmake_script)
                print("message(STATUS ${CMAKE_BUILD_TYPE})", file=cmake_script)
                cmake_script.flush()
                cmake_cmd = ["cmake", "-P", cmake_script.name]
                cfg = subprocess.check_output(cmake_cmd, text=True)
                assert cfg.startswith("--"), "unexpected build type"
                cfg = cfg[2:].strip()
                if len(cfg) > 0:
                    return cfg
        if self.debug:
            return "Debug"
        return "Release"

    def build_extension(self, ext):

        cmake_args = []

        # Check if we have a specific configuration file to override CMake cache
        if config_cmake_path.exists():
            cmake_args += ["-C", f"{config_cmake_path.absolute()}"]

        cmake_args += ["-DPYTHON_EXECUTABLE=" + sys.executable]

        # this is where setuptools will build the module
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath("")))
        setup_build_directory, setup_package_directory = os.path.split(extdir)
        assert setup_package_directory == self.package
        cmake_args += ["-DCMAKE_INSTALL_PREFIX=" + setup_build_directory]
        build_args = ["--target", "install"]  # we will build the install target

        # config defaults to Debug on Windows if not explicitly passed to the build tool
        cfg = self._retrieve_cmake_build_type()
        build_args += ["--config", cfg]

        # '--' will pass remaining options to the native build tool
        build_args += ["--"]

        if platform.system() == "Windows":
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["-j%d" % cpu_count()]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )

        # this is the build directory for compiled modules
        build_dir = os.path.abspath(self.build_temp) + "-cmake"
        if not os.path.exists(build_dir):
            os.makedirs(build_dir)

        # cmake call to set up the native build
        cmake_cmd = ["cmake", ext.sourcedir] + cmake_args
        subprocess.check_call(cmake_cmd, cwd=build_dir, env=env)

        # cmake is then used to call the native build tool
        # just re-run cmake --build . --target install
        build_cmd = ["cmake", "--build", "."] + build_args
        subprocess.check_call(build_cmd, cwd=build_dir)


setup(
    name=package_name,
    version_config=True,
    author="various contributors",
    author_email="anr-charms@brgm.fr",
    description="A parallel multiphase multicomponents simulator.",
    long_description="",
    url="https://charms.gitlabpages.inria.fr/ComPASS/",
    license="GPLv3/CeCILLv2.1",
    packages=[
        package_name,
        f"{package_name}.eos",
        f"{package_name}.ghosts",
        f"{package_name}.io",
        f"{package_name}.linalg",
        f"{package_name}.petrophysics",
        f"{package_name}.petrophysics.models",
        f"{package_name}.physics",
        f"{package_name}.physics.water2ph",
        f"{package_name}.schemes",
        f"{package_name}.simulation",
        f"{package_name}.utils",
        f"{package_name}.wells",
    ],
    ext_package=package_name,
    ext_modules=[
        CMakeExtension("water2ph")
    ],  # name (first argument) must match the name of the exported pybind11 module
    cmdclass=dict(build_ext=CMakeBuild),
    package_data={package_name: [version_file],},
    zip_safe=False,
    setup_requires=["setuptools-git-versioning"],
)
