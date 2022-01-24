"""
Performs incremental dynamic analysis using a single pair of records with 10 runs.

Note: for full IDA, a suite of records must be increased.
"""
from pathlib import Path
from rcmrf import RCMRF
import pickle

from utils.utils import get_time, get_start_time


start_time = get_start_time()

# Directories
main_dir = Path.cwd()
input_dir = main_dir / "Inputs"
materials_file = input_dir / "materials.csv"
outputsDir = main_dir / "outputs/ida"
loads_file = input_dir / "action.csv"
modal_analysis_path = main_dir / "outputs/modal/MA.json"

section_file = input_dir / "hinge_models.pickle"
with open(section_file, "rb") as f:
    section_file = pickle.load(f)

# GM directory
gmdir = main_dir / "gm_ida"

gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt"]

# RCMRF inputs
hingeModel = "hysteretic"
system = "space"
analysis_type = ["IDA"]
flag3d = True
# Export outputs as a single file (this time :))
export_at_each_step = False
analysis_time_step = 0.01
max_runs = 10

m = RCMRF(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
          analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
          export_at_each_step=export_at_each_step, modal_analysis_path=modal_analysis_path, max_runs=max_runs)

m.wipe()
m.run_model()

# Wipe the model
m.wipe()

# Time it
get_time(start_time)
