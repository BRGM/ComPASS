import sys
import os
from pathlib import Path
import importlib
import site
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-a", "--add",
                  action="store_true", dest="add", default=False,
                  help="add module.pth file")

parser.add_option("--remove",
                  action="store_true", dest="remove", default=False,
                  help="remove module.pth files")

options, args = parser.parse_args()


def add_pth_file_to_user_site(mname, mpath, message=True):
    mpath = Path(mpath)
    if not mpath.is_dir():
        print(f'{mpath} directory does not exist')
        sys.exit(-1)
    user_site = Path(site.getusersitepackages())
    if not user_site.is_dir():
        assert not user_site.exists()
        user_site.mkdir()
    pthfile = user_site / f'{mname}.pth'
    if message:
        print(f'Setting {mname} module path to {mpath} via user site file: {pthfile}')
    with pthfile.open('w') as f:
        print(mpath.resolve(), file=f)


def add(mname, mpath):
    assert sys.version_info.major == 3
    # ModuleNotFoundError is new version 3.6.
    ModuleNotFound = ImportError
    if sys.version_info.minor > 5:
        ModuleNotFound = ModuleNotFoundError
    try:
        mfound = importlib.import_module(mname)
    except ModuleNotFound:
        mfound = None 
    if mfound is not None:
        mfoundpath, mfile = os.path.split(mfound.__file__)
        if mfile=='__init__.py':
            packagepath, packagename = os.path.split(mfoundpath)
            assert packagename==mname
            mfoundpath = packagepath
        if os.path.normpath(mfoundpath)!=os.path.normpath(mpath):
            print('WARNING: changing', mname, 'directory',
                  'from:', mfoundpath, 'to:', mpath)
            add_pth_file_to_user_site(mname, mpath, message=False)
        else:
            print('ComPASS.pth is already set to:', mpath)
    else:
        add_pth_file_to_user_site(mname, mpath)


def remove_pth_file_from_user_site(mname):
    user_site = site.getusersitepackages()
    if os.path.isdir(user_site):
        pthfile = os.path.join(user_site, mname + '.pth')
        if os.path.isfile(pthfile):
            os.remove(pthfile)
        else:
            print('no pth file for', mname)
    else:
        print('no user site directory found')

if options.add:
    if len(args)!=2:
        sys.exit(-1)
    add(*args)
elif options.remove:
    for mname in args:
        remove_pth_file_from_user_site(mname)
