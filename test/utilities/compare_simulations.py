import sys
import yaml
from yaml_plot import SimulationLog
from pathlib import Path

dir1, dir2 = sys.argv[1], sys.argv[2]


dir1, dir2 = Path(dir1), Path(dir2)
if not dir1.is_dir():
    raise IOError(f"Could not find {str(dir1)} directory.")
if not dir2.is_dir():
    raise IOError(f"Could not find {str(dir2)} directory.")
datapath1 = dir1 / "simulation_log.yaml"
datapath2 = dir2 / "simulation_log.yaml"
if not datapath1.exists():
    raise IOError(f"Could not find {str(datapath1)} file.")
if not datapath2.exists():
    raise IOError(f"Could not find {str(datapath2)} file.")

with open(datapath1, "r") as f:
    data1 = yaml.safe_load(f)
with open(datapath2, "r") as f:
    data2 = yaml.safe_load(f)

simlog1 = SimulationLog(data1)
simlog2 = SimulationLog(data2)
i = 0
for ts_dict1, ts_dict2 in zip(simlog1.time_steps, simlog2.time_steps):
    if len(ts_dict1) != len(ts_dict2):
        print(f"Timeline divergence at time step {i}")
        print(ts_dict1, ts_dict2)
    i += 1
