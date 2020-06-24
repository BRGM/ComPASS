import sys
import re
from pathlib import Path


def wheel_tag():
    # FIXME: is there a nicer way to do that ???
    from distutils import dist
    from wheel.bdist_wheel import bdist_wheel

    return "-".join(bdist_wheel(dist.Distribution()).get_tag())


class PackageInfo:
    def __init__(self, directory):
        directory = Path(directory)
        assert directory.exists()
        self.directory = directory
        init_file = directory / "__init__.py"
        assert init_file.is_file()
        for name, value in re.findall(
            r"""__([a-z]+)__ = "([^"]+)""", init_file.read_text(encoding="utf8")
        ):
            setattr(self, name, value)

    @property
    def name(self):
        return self.directory.name
