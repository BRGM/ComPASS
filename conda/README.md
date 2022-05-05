This directory contains yaml files to set
[conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#managing-environments).
To use them just run:

```shell
conda env create -f file.yml [-n my_fancy_name]
```

They come in two flavors:
- `compass.yml` defines an environment where you can compile ComPASS
  and which can serve for development purposes
  (you will have first to git clone the ComPASS repository).
- `compass-latest.yml` is the same environment as above but it will also
  install the latest version of ComPASS available on github so that
  you will be able to readily use it an run a ComPASS script.
