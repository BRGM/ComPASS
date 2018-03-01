import sys
import os
import importlib
import site

def add_pth_file_to_user_site(mname, mpath):
    if not os.path.isdir(mpath):
        print('directory does not exist')
        sys.exit(-1)
    user_site = site.getusersitepackages()
    if not os.path.isdir(user_site):
        assert not os.path.exists(user_site)
        os.makedirs(user_site)

    with open(os.path.join(user_site, mname + '.pth'), 'w') as f:
        print(mpath, file=f)


args = sys.argv[1:]
if len(args)!=2:
    sys.exit(-1)
# module name and path
mname, mpath = args

try:
    mfound = importlib.import_module(mname)
except ModuleNotFoundError:
    mfound = None 

if mfound is not None:
    mfoundpath, _ = os.path.split(mfound.__file__)
    if os.path.normpath(mfoundpath)!=os.path.normpath(mpath):
        print('WARNING: changing', mname, 'directory',
              'from:', mfoundpath, 'to:', mpath)
        add_pth_file_to_user_site(mname, mpath)
else:
    add_pth_file_to_user_site(mname, mpath)
