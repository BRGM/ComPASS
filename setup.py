import sys

try:
    from skbuild import setup
except ImportError:
    print("scikit-build is required to build from source!", file=sys.stderr)
    print("Install it running: python -m pip install scikit-build.", file=sys.stderr)
    sys.exit(1)

import setuptools_scm as scm
from pathlib import Path
import platform
from datetime import datetime

package_name = "ComPASS"
build_info_file = "version_info"

# Generate info about package compilation
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

setup(
    name=package_name,
    use_scm_version={"write_to": f"{package_name}/_version.py"},
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
    package_data={
        package_name: [build_info_file],
    },
    setup_requires=["setuptools_scm"],
)
