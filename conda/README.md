This directory contains yaml files to set
[conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments).
To use them just run:

```shell
conda env create -f file.yml [-n my_fancy_name]
```

They come in several flavors:
- `compass.yml` defines an environment where you can compile ComPASS
  and which can serve for development purposes
  (you will have first to git clone the ComPASS repository).
- `compass-latest.yml` is the same environment as above with a fresh
  installation of the latest version of ComPASS that will be downloaded
  from github: at the end of the process you will be ready to
  run a ComPASS script.
- `compass-v4.4.1.yml` as above but with ComPASS v4.4.1: at the end of
  the process you will be ready to run a ComPASS script with this
  specific version.
