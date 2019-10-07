import sys
import optparse
import multiprocessing

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
parser.add_option("--parallel",
                  action="store_true", dest="parallel_run",
                  help="will run as parallel job with the number of available procs using mpirun to process simulation script or pytest-xdist if --pytest is present")
parser.add_option("--pytest",
                  action="store_true", dest="pytest_run",
                  help="will run pytest with `nproc` procs if --parallel option is present, this will override some options")
parser.add_option("--bash",
                  action="store_true", dest="bash_session",
                  help="will enter a bash session executing the given command and overriding all other options")
parser.add_option("--postprocess",
                  action="store_true", dest="postprocess_run",
                  help=("will run the ComPASS postprocess script overriding all other"
                        "options but the bash session that is priority"
                        " BEWARE that if you want to pass specific command"
                        " to the postprocess script, you have to provide them with a slash, i.e. the MS way (e.g. /h)"))
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
    if options.bash_session:
        cmd.append('/bin/bash')
        cmd.extend(args)
    elif options.postprocess_run:
        cmd.append('python3 -m ComPASS.postprocess')
        cmd.extend(['-'+s[1:] if s.startswith('/') else s for s in args])
    elif options.pytest_run:
        cmd.append('python3 -m pytest')
        if options.parallel_run:
            cmd.append('-n %d' % multiprocessing.cpu_count())
    else:
        if options.parallel_run:
            cmd.append('mpirun -n %d' % multiprocessing.cpu_count())
        cmd.append('python3')
        cmd.extend(args)
    #print('Command:', ' '.join(cmd))
    print(' '.join(cmd), file=f)
