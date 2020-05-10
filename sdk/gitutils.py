import sys
import subprocess
import shlex


def run(cmd, verbose=True, output_stderr_on_error=True):
    if verbose:
        print(cmd)
    res = subprocess.run(shlex.split(cmd), capture_output=True)
    if output_stderr_on_error and res.returncode != 0:
        print(res.stderr.decode("utf-8"))
    res.check_returncode()
    return res.stdout.strip().decode("utf-8")


def check_clean_head():
    check = run("git status -uno --short", verbose=False)
    if len(check.strip()) > 0:
        print()
        print("This script needs to start from a clean HEAD state.")
        print("You may want to use: git reset --hard")
        print()
        sys.exit(-1)
