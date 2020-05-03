import sys
import yaml
import subprocess
import shlex
import tempfile


def run(cmd, verbose=True):
    if verbose:
        print(cmd)
    res = subprocess.run(shlex.split(cmd), capture_output=True)
    if res.returncode != 0:
        print(res.stderr)
    res.check_returncode()
    return res.stdout.strip().decode("utf-8")


assert len(sys.argv) == 2, f"Syntax is: {sys.argv[0]} yaml_history_file"
with open(sys.argv[1]) as f:
    commits = yaml.safe_load(f)

shaone_aliases = [commit["sha1"] for commit in commits]
shaone_aliases.insert(0, "HEAD")
sha1 = []
for alias in shaone_aliases:
    sha1.append(run(f"git rev-parse {alias}", verbose=False))

for ancestor, commit in zip(sha1[:-1], sha1[1:]):
    base = run(f"git merge-base {ancestor} {commit}", verbose=False)
    assert base == ancestor, f"{ancestor} is not an ancestor of {commit}"

messages = [commit["message"] for commit in commits]
for k in range(len(messages)):
    ancestor, commit = sha1[k : k + 2]
    print()
    print(f"Going from {ancestor[:6]} to {commit[:6]}")
    modified_files = run(f"git diff --name-status {ancestor} {commit}", verbose=False)
    modified_files = [l.strip().split() for l in modified_files.split("\n")]
    # print(modified_files)
    assert all([f[0][0] in "AMRD" for f in modified_files]), "Unknown status code"
    for f in modified_files:
        status = f[0]
        if status == "M" or status == "A":
            filename = f[1]
            run(f"git checkout {commit} -- {filename}")
        elif status == "D":
            filename = f[1]
            run(f"git rm {filename}")
        elif status[0] == "R":
            source, dest = f[1:]
            run(f"git mv {source} {dest}")
            run(f"git checkout {commit} -- {dest}")
        else:
            assert (
                False
            ), f"Stuck at {commit} with diff '{' '.join(f)}' I do not kwow what to do!"
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(messages[k].encode("utf-8"))
        tmp.flush()
        run(f"git commit --no-verify -F {tmp.name}")
