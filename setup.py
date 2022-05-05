# dirty workaround around sysconfig.get_platform bug on MacOSX
import os
import sysconfig

platform_tag = sysconfig.get_platform()
if platform_tag.startswith("macosx"):
    assert all(platform_tag.split("-")), "Cannot use platform information!"
    os.environ["_PYTHON_HOST_PLATFORM"] = platform_tag

# Generate info about package compilation
import setuptools_scm as scm
from pathlib import Path
import platform
from datetime import datetime

package_name = "ComPASS"
build_info_file = "version_info"

with Path(f"./{package_name}/{build_info_file}").open("w") as f:
    print(
        f"{package_name} {scm.get_version()}",
        file=f,
    )
    print(f"built {datetime.now().isoformat()}", file=f)
    print(f"with python {platform.python_version()} on {platform.node()}", file=f)
    print(f"build system info: {platform.platform()}", file=f)

# Clean local eos directory (!_ is to keep __init__.py)
for f in Path(f"{package_name}/eos").glob("[!_]*"):
    f.unlink()

from skbuild import setup

setup(
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
    package_data={
        package_name: [build_info_file],
    },
)
