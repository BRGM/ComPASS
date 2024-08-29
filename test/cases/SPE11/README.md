# 11th SPE CSP

This repository accompanies the 11th Society of Petroleum Engineers
Comparative Solution Project, see https://spe.org/csp/, or the
[github project](https://github.com/Simulation-Benchmarks/11thSPE-CSP).


A team has participated with the ComPASS code (v4).
Only the case _b_ has been implemented.
The results have been obtained with the tag v4.5.7.

## Create the test b mesh

You need the [geometry files from github](https://github.com/Simulation-Benchmarks/11thSPE-CSP/tree/main/geometries),
it creates the spe11b_structured.msh file.
```shell
cd SPE11b/geometries/
python make_structured_mesh.py --variant B -nx 420 -ny 140
```

Then the ``msh2compass.py`` utility convert the gmsh mesh into a compass mesh
(it is called from the SPE11b.py file).

## Execute the SPE11b test

We advise to run in parallel, from the _SPE11b_ directory, execute:
```shell
mpirun -n 128 python SPE11_b.py
```
The outputs are stored in the _output-SPE11_b_ directory.
It also creates _cells_rkt.vtu_ in the same directory (__output-SPE11b_).

## Convert to the SPE11 format

The outputs of the SPE11 benchmark must be in a specific format.
Several utilities have been developped to convert into the good format.

You can execute `sbatch report.sbatch conda_env job_id` on Leto,
where:
* _conda_env_ is the _ComPASS-conda_env_ environment (which must contain Paraview 5.12),
* _job_id_ is the job identifier of the ComPASS execution.

**OR** from the _SPE11b_ directory, execute:
```shell
python -m ComPASS.postprocess -s -C -t "second" \$PATH/output-SPE11_b
pvpython pv_extract.py "\$PATH/output-SPE11_b/paraview/states.pvd" \$PATH/output-SPE11_b/cells_rkt.vtu -o "\$PATH/output-SPE11_b/extracts"
python make_report.py 'b' --src \$PATH/output-SPE11_b/extracts --out \$PATH/output-SPE11_b/report
```
where ``pv_extract`` needs Paraview 5.12.
