# dirty workaround around sysconfig.get_platform bug on MacOSX
import os
import sysconfig

platform_tag = sysconfig.get_platform()
if platform_tag.startswith("macosx"):
    assert all(platform_tag.split("-")), "Cannot use platform information!"
    os.environ["_PYTHON_HOST_PLATFORM"] = platform_tag

import setuptools_scm as scm
from pathlib import Path
import platform
from datetime import datetime

print(f"ComPASS {scm.get_version()}")
print(f"built {datetime.now().isoformat()}")
print(f"with python {platform.python_version()} on {platform.node()}")
print(f"build system info: {platform.platform()}")
