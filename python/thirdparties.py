

# The onlyt thing that this does is that it adds third parties path into system path
# FIXME: Find a more elegant way to manage third parties packages...

import sys
import os

thirdparty_directory = '../thirdparty'
thirdparties_directories = [
        os.path.abspath(os.path.join(thirdparty_directory, s)) for s in
        ('meshtools',)
        ]

for path in thirdparties_directories:
    assert os.path.exists(path)
    assert os.path.isdir(path)
    if not path in sys.path:
        sys.path.append(path)

