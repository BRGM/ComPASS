import sys
from subprocess import CalledProcessError
from pathlib import Path

import click

from gitutils import run, check_clean_head


def check_file_paths(files, skip_missing_files=False):
    kept = []
    for path in [Path(path) for path in files]:
        if path.exists():
            kept.append(path.as_posix())
        else:
            if not skip_missing_files:
                raise FileNotFoundError(path.as_posix())
    return kept


def reformat_files(files, commit_message=None, verbose=False, skip_missing_files=False):
    files = check_file_paths(files, skip_missing_files)
    iterations_left = len(files)
    still_files_to_format = iterations_left > 0
    # WARNING: potentially infinite loop if pre-commit fails to format some files
    while iterations_left > 0 and still_files_to_format:
        try:
            run(
                f"pre-commit run --file {' '.join(files)}",
                verbose=False,
                output_stderr_on_error=False,
            )
        except CalledProcessError:
            if verbose:
                print("Files that were reformatted:")
                print(run("git diff --name-only", verbose=False))
        else:
            still_files_to_format = False
        iterations_left -= 1
    for path in files:
        run(f"git add -u {path}", verbose=False)
    need_commit = len(run("git status -uno --short", verbose=False).strip()) > 0
    if need_commit and commit_message is not None:
        try:
            run(
                f"git commit -m '{commit_message}'",
                verbose=False,
                output_stderr_on_error=True,
            )
        except CalledProcessError:
            print()
            print("Something went wrong...")
            print("Aborting.")
            print()
            sys.exit(-1)
    else:
        print("No commit needed.")


def create_format_commit(selection=None, format_source=None):

    check_clean_head()

    configuration_files = [".pre-commit-config.yaml", ".clang-format"]

    if format_source is not None:
        run(f"git checkout {format_source} -- {' '.join(configuration_files)}")

    try:
        precommit_version = run("pre-commit --version", verbose=False)
    except FileNotFoundError:
        print()
        print("pre-commit script not found !")
        print("Install pre-commit with:")
        print(f"{sys.executable} -m pip install pre-commit")
        print("or check: https://pre-commit.com/")
        print()
        sys.exit(-1)
    print(f"Using: {precommit_version}")

    for conffile in configuration_files:
        if not Path(conffile).is_file():
            print()
            print(f"Could not find configuration file: {conffile}")
            print()

    skipped = []
    if selection is None:
        skipped_directories = ["thirdparty"]
        for dirname in skipped_directories:
            path = Path(dirname)
            if path.is_dir():
                skipped.append(path)
        tracked_files = run("git ls-tree -r --name-only HEAD", verbose=False)
        tracked_files = [l.strip() for l in tracked_files.split("\n")]
    else:
        tracked_files = selection

    for conffile in configuration_files:
        if conffile not in tracked_files:
            tracked_files.append(conffile)
    checked_files = []
    for path in tracked_files:
        path = Path(path)
        if not any([skipped_path in path.parents for skipped_path in skipped]):
            checked_files.append(path)
    commit_message = f"Reformat code\n\nUsing {precommit_version}"
    if format_source is not None:
        commit_message += f"\nwith {' '.join(configuration_files)} from {format_source}"
    reformat_files(checked_files, commit_message, skip_missing_files=True)


def collect_affected_files(base, target="HEAD"):
    merge_base = run(f"git merge-base {base} {target}", verbose=False).strip()
    history = run(f"git log --oneline {merge_base}..{target}", verbose=False)
    sha1s = [line.strip().split()[0] for line in history.split("\n")][::-1]
    files = set()
    for sha1 in sha1s:
        files.update(
            [
                filename.strip()
                for filename in run(
                    f'git show --pretty="" --name-only {sha1}', verbose=False
                ).split("\n")
            ]
        )
    return list(files)


def extract_history(base, target="HEAD"):
    merge_base = run(f"git merge-base {base} {target}", verbose=False).strip()
    history = run(f"git log --oneline {merge_base}..{target}", verbose=False)
    sha1s = [line.strip().split()[0] for line in history.split("\n")][::-1]
    messages = [
        run(f"git log --format=%B -n 1 {sha1}", verbose=False) for sha1 in sha1s
    ]
    statuses = [
        [
            status_line.strip().split()
            for status_line in run(
                f'git show --pretty="" --name-status {sha1}', verbose=False
            ).split("\n")
        ]
        for sha1 in sha1s
    ]
    return merge_base, sha1s, messages, statuses


@click.command()
@click.option(
    "--format-source",
    help=(
        "A git valid reference that can provide the format configuration files: "
        ".pre-commit-config.yaml and .clang-format"
    ),
)
@click.argument("base")
@click.argument("target")
def format_history(base, target, format_source):
    """
    Reformat all commits between BASE and TARGET using the pre-commit tool.
    The first commit is a special formating commit that reformat
    all the files in BASE.
    The new commits will be made with a detached HEAD state.
    Don't forget to create a new branch to save the results.
    """
    affected_files = collect_affected_files(base, target)
    start, sha1s, messages, statuses = extract_history(base, target)
    run(f"git checkout --detach {start}")
    print(f"Reformating the merge-base...")
    create_format_commit(selection=affected_files, format_source=format_source)
    for sha1, message, status in zip(sha1s, messages, statuses):
        print()
        print(f"-- Reformating: {sha1}", "-" * 57)
        print()
        print(message)
        kept = []
        for line in status:
            status = line[0]
            if status[0] not in "AMRD":
                raise ValueError(f"Unknown status code: {' '.join(line)}")
            if status == "M" or status == "A":
                kept.append(line[1])
            elif status == "D":
                run(f"git rm {line[1]}", verbose=False)
            else:
                assert status[0] == "R"
                source, dest = line[1:]
                run(f"git mv {source} {dest}", verbose=False)
                kept.append(dest)
        if len(kept) > 0:
            run(f"git checkout {sha1} -- {' '.join(kept)}", verbose=False)
        reformat_files(kept, message)
    print()
    print("-- Done!", "-" * 72)
    print()
    print("You are now in a detached HEAD state.")
    print("Use: git checkout -b new_branch_name")
    print("to create a new branch and checkt it out.")
    print()


base = "v4.1"
target = "HEAD"
format_source = "v4.2"

if __name__ == "__main__":
    format_history()
