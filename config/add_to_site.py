import sys
import os
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
    if not os.path.isdir(mpath):
        print('directory does not exist')
        sys.exit(-1)
    user_site = site.getusersitepackages()
    if not os.path.isdir(user_site):
        assert not os.path.exists(user_site)
        os.makedirs(user_site)
    if message:
        print('setting', mname, 'module path to', mpath)
    pthfile = os.path.join(user_site, mname + '.pth')
    with open(pthfile, 'w') as f:
        print(mpath, file=f)


def add(mname, mpath):
    try:
        mfound = importlib.import_module(mname)
    except ModuleNotFoundError:
        mfound = None 
    if mfound is not None:
        mfoundpath, _ = os.path.split(mfound.__file__)
        packagepath, packagename = os.path.split(mfoundpath)
        if packagename==mname:
            mfoundpath = packagepath
        if os.path.normpath(mfoundpath)!=os.path.normpath(mpath):
            print('WARNING: changing', mname, 'directory',
                  'from:', mfoundpath, 'to:', mpath)
            add_pth_file_to_user_site(mname, mpath, message=False)
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
