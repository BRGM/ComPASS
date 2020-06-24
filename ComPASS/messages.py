from . import mpi


def banner_message(s):
    return r"""
%s !
%s ! %%s
%s !
""" % (
        s,
        s,
        s,
    )


def message_with_banner(banner, message, abort, explanation=None):
    print(banner_message(banner) % message)
    if explanation is not None:
        print(explanation)
    if abort:
        mpi.abort()


@mpi.on_master_proc
def error(message, abort=True, explanation=None):
    message_with_banner("ERROR", message, abort, explanation)


@mpi.on_master_proc
def warning(message, abort=False, explanation=None):
    message_with_banner("WARNING", message, abort, explanation)


@mpi.on_master_proc
def deprecation(message, abort=False, explanation=None):
    message_with_banner("DEPRECATION", message, abort, explanation)
