import ComPASS.mpi as mpi

def banner_message(s):
    return r'''
%s !
%s ! %%s
%s !
''' % (s, s, s)

def message_with_banner(banner, message, abort):
    print(banner_message(banner) % message)
    if abort:
        mpi.abort()

@mpi.on_master_proc
def error(message, abort=True):
    message_with_banner('ERROR', message, abort)

@mpi.on_master_proc
def warning(message, abort=False):
    message_with_banner('WARNING', message, abort)

@mpi.on_master_proc
def deprecation(message, abort=False):
    message_with_banner('DEPRECATION', message, abort)
   