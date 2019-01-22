import sys
#import glob, os, re
import optparse
#import numpy as np

parser = optparse.OptionParser(
usage = """usage: [options] [path]

Parameters:
  path           path to a simulation script"""
)
parser.add_option("--customize-session",
                  action="store", type="string", dest="customization_file",
                  default="/etc/customize-session",
                  help=optparse.SUPPRESS_HELP # path to the file where session customization will be written
)
parser.add_option("--simulation-process",
                  action="store", type="string", dest="process_file",
                  default="/etc/simulation-process",
                  help=optparse.SUPPRESS_HELP # the user uid that will be given to the compass user
)
parser.add_option("--compass-uid",
                  action="store", type="int", dest="uid", default=None,
                  help="the user uid that will be given to the compass user")
parser.add_option("-p", "--parallel",
                  action="store_true", dest="parallel_run",
                  help="will run as parallel jobs using mpirun with `nproc` procs")
options, args = parser.parse_args()

#print('Writing customization file:', options.customization_file)
with open(options.customization_file, 'w') as f:
    if options.uid:
        print(r'''
if [ `id -u $COMPASS_USER` -ne %d ]
then
    usermod -u %d compass
fi
''' % ((options.uid,)*2) , file=f)

# This is the command to be substituted as the current process
#print('Writing process file:', options.process_file)
with open(options.process_file, 'w') as f:
    cmd = []
    if options.parallel_run:
        cmd.append('mpirun -n `nproc`')
    cmd.append('python3')
    cmd.extend(args)
    print(' '.join(cmd), file=f)
