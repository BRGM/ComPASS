import sys
import re
from pathlib import Path

import PackageInfo


def get_wheel_filename(package_directory, prereltag=None):
    from wheel.bdist_wheel import safer_name, safer_version

    package = PackageInfo.PackageInfo(package_directory)
    wheel_filename = "{name}-{version}{prereltag}-{wheeltag}.whl".format(
        name=safer_name(package.name),
        version=safer_version(package.version),
        prereltag="" if prereltag is None else safer_name(prereltag),
        wheeltag=PackageInfo.wheel_tag(),
    )
    return wheel_filename


def check_wheel_path(directory, wheel_filename):
    wheel_path = Path(directory) / wheel_filename
    assert wheel_path.is_file(), f"Could not find {wheel_path}"
    return wheel_path


def check_pip_executable(pipexe):
    pipexe = Path(pipexe)
    assert pipexe.exists()
    return pipexe


if __name__ == "__main__":
    assert len(sys.argv) > 3
    if len(sys.argv) > 4:
        wheel_file = get_wheel_filename(sys.argv[2], sys.argv[3])
    else:
        wheel_file = get_wheel_filename(sys.argv[2])
    pipexe = check_pip_executable(sys.argv[1])
    wheel_path = check_wheel_path(sys.argv[-1], wheel_file)
    print(pipexe, "install", wheel_path)
