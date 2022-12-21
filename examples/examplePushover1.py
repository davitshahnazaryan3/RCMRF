"""
Example file to run static analysis on a 3D building model
Run after running exampleModal.py to use load pattern type 2 (i.e., 1st mode proportional)
If MA.json exists, add to pushover folder, or give path as input
"""

from pathlib import Path
from rcmrf import RCMRF
import pickle

from utils.utils import get_time, get_start_time


# record time to run model
start_time = get_start_time()

# Directories
main_dir = Path.cwd()
input_dir = main_dir / "RCMRF"
materials_file = input_dir / "materials.csv"
outputsDir = main_dir / "outputs/RCMRF_temp"
loads_file = input_dir / "action_ancona.csv"
modal_analysis_path = main_dir / "outputs/RCMRF_temp/MA.json"

section_file = input_dir / "hinge_models_Ancona.pickle"
with open(section_file, "rb") as f:
    section_file = pickle.load(f)

# RCMRF inputs
analysis_type = ["PO"]
flag3d = True
direction = 0
period_assignment = {"x": 1, "y": 0}

m = RCMRF(section_file, loads_file, materials_file, outputsDir, analysis_type=analysis_type, flag3d=flag3d,
          direction=direction, modal_analysis_path=modal_analysis_path, period_assignment=period_assignment)

m.wipe()
m.run_model()

# Time it
get_time(start_time)

# Wipe the model
m.wipe()
