This directory contains yaml files to set
[conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments).

To use them just run (note that `mamba` can replace `conda`):

```shell
conda env create -f file.yml [-n my_fancy_name]
```

Environments come in several flavors.

Two environments are designed for development purposes (you want to compile and modify ComPASS):
- `compass-stable.yml` is an environment with version constraints set, it is known to be stable and is regularly tested in CI.
- `compass-pristine.yml` is an environment with very few constraints set, so it will install the latest package version from conda-forge. It is regularly tested in CI but is not guaranteed to lead to a successful compilation (it is mainly used to detect breaking API changes in package releases).

If you only want to run ComPASS scripts, you can use any of the two previous environments to have ComPASS compiled once (and for all).

If you want the latest version from the main branch, just add the following line to pip requirements:
```
  - git+https://github.com/BRGM/ComPASS.git
```

You can also add a specific version tag, e.g.:
```
  - git+https://github.com/BRGM/ComPASS.git@v4.5.6
```

Be aware that older versions may need a bit of manual tuning to achieve a successful build.
