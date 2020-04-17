import sys
from multiprocessing import Pool
from subprocess import run, CalledProcessError
from pathlib import Path
import time


def error_diagnosis(filename, result):
    print(
        f"""
===== ERROR === ERROR === ERROR === ERROR === ERROR === ERROR === ERROR =====

{filename} execution failed

with standard output:
=====================

{result.stdout}

and standard error:
===================

{result.stderr}

"""
    )


def time_script(filename, tag=None, timeout=600):
    start = time.time()
    filename = Path(filename)
    assert filename.exists()
    message = f"{'Launching:':13s} {filename.name}"
    if tag is not None:
        message = message + f" which is job {tag}"
    print(message)
    try:
        result = run(
            [sys.executable, filename.as_posix()],
            capture_output=True,
            encoding="utf-8",
            timeout=timeout,
        )
    except TimeoutExpired:
        error_diagnosis(filename, result)
        print(
            f"{'Timed-out:':13s} {filename.name} timed out after {tiemout}s.\n",
            file=sys.stderr,
        )
        raise
    elapsed = time.time() - start
    try:
        result.check_returncode()
    except CalledProcessError:
        error_diagnosis(filename, result)
        raise
    return filename, elapsed


if __name__ == "__main__":
    filenames = []
    for d in sys.argv[1:]:
        path = Path(d)
        if path.exists():
            if path.is_dir():
                filenames.extend(path.glob("*.py"))
            elif path.is_file():
                filenames.append(path)
    nb_scripts = len(filenames)
    print(f"Collected {nb_scripts} scripts.")
    with Pool(maxtasksperchild=1) as pool:
        jobs = [
            pool.apply_async(time_script, (filename,), {"tag": f"{fk+1}/{nb_scripts}"})
            for fk, filename in enumerate(filenames)
        ]
        for jk, job in enumerate(jobs):
            try:
                filename, elapsed = job.get()
            except CalledProcessError:
                sys.exit(-1)
            print(
                f"{'Happy-ending:':13s} {filename} ran in {elapsed:.3f} s (indicative time for job {jk+1}/{nb_scripts})"
            )
