import sys
import os
from pathlib import Path
import click
import psutil


def write_customization_file(filename, uid=None):
    with open(filename, "w") as f:
        if uid is not None:
            print(
                f"""
    if [ `id -u $COMPASS_USER` -ne {uid} ]
    then
        usermod -u {uid} compass
    fi""",
                file=f,
            )


class SessionScript:
    def __init__(self, filename):
        self.filepath = Path(filename)
        self.unlink()

    def __enter__(self):
        self.file = self.filepath.open("a")
        return self.file

    def __exit__(self, *args):
        self.file.close()

    def unlink(self):
        # UPGRADE: from python 3.8, use the missing_ok option
        if self.filepath.exists():
            self.filepath.unlink()


@click.group()
@click.option(
    "--customize-session",
    "customization_file",
    default="/tmp/customize-session",
    hidden=True,
    help="path to the file where session customization will be written",
)
@click.option(
    "--session-script", "session_script", default="/tmp/script", hidden=True,
)
@click.option(
    "--compass-uid",
    "-u",
    "uid",
    default=None,
    type=int,
    help="the user uid that will be given to the compass user",
)
@click.pass_context
def compass(ctx, customization_file, session_script, uid):
    write_customization_file(customization_file, uid)
    ctx.obj = SessionScript(session_script)


@compass.command()
@click.option(
    "--parallel",
    is_flag=True,
    default=False,
    help=f"run in parallel with the number of available procs ({psutil.cpu_count(logical=False)}) using mpirun to process simulation script",
)
@click.argument("script", nargs=1)
@click.pass_obj
def run(session_script, parallel, script):
    """Run a ComPASS python script."""
    with session_script as f:
        cmd = f"python3 {script}"
        if parallel:
            print(f"mpirun -n {psutil.cpu_count(logical=False)} {cmd}", file=f)
        else:
            print(cmd, file=f)


@compass.command(context_settings={"ignore_unknown_options": True})
@click.argument("cmd", nargs=-1, type=click.UNPROCESSED)
@click.pass_obj
def bash(session_script, cmd):
    """Enter a bash session otherwise
    execute arguments as a command in a bash session.
    If your command
    has double quotes use simple quote to escape them."""
    session_script = session_script
    with session_script as f:
        if len(cmd) > 0:
            print(" ".join(cmd), file=f)


@compass.command()
@click.argument("script", nargs=1)
@click.pass_obj
def exec(session_script, script):
    """Execute a given script in a bash session"""
    script = Path(script)
    if script.exists():
        with script.open() as f:
            with session_script as g:
                for line in f:
                    g.write(line)
    else:
        click.echo(f"Could not find: {script}", err=True)


# We can not use compass.add_command(ComPASS.postprocess.postprocess)
# because the comand is fired from the session script
# This is why we capture the help option
@compass.command(context_settings={"ignore_unknown_options": True})
@click.option("--help", "request_help", is_flag=True, default=False)
@click.argument("args", nargs=-1, type=click.UNPROCESSED)
@click.pass_obj
def postprocess(session_script, request_help, args):
    with session_script as f:
        cmd = "python3 -m ComPASS.postprocess"
        if request_help:
            print(f"{cmd} --help", file=f)
        else:
            cmd = [cmd]
            cmd.extend(args)
            print(f"{' '.join(cmd)}", file=f)


@compass.command()
@click.option(
    "--parallel",
    is_flag=True,
    default=False,
    help=f"Use pytest-xdist with the number of available procs ({psutil.cpu_count(logical=False)})",
)
@click.pass_obj
def test(session_script, parallel):
    """Run pytest."""
    with session_script as f:
        pytest = "python3 -m pytest"
        if parallel:
            print(f"{pytest} -n {psutil.cpu_count(logical=False)}", file=f)
        else:
            print(f"{pytest}", file=f)


if __name__ == "__main__":
    compass(prog_name="compass")
