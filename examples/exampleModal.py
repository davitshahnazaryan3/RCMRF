"""
Example file to run static analysis on a 3D building model
"""

from pathlib import Path
from rcmrf import RCMRF
import pickle

from utils.utils import get_time, get_start_time


# record time to run model
start_time = get_start_time()

# Directories
main_dir = Path.cwd()
input_dir = main_dir / "Inputs"
materials_file = input_dir / "materials.csv"
outputsDir = main_dir / "outputs/modal"
loads_file = input_dir / "action.csv"

section_file = input_dir / "hinge_models.pickle"
with open(section_file, "rb") as f:
    section_file = pickle.load(f)

# RCMRF inputs
analysis_type = ["MA"]
flag3d = True

m = RCMRF(section_file, loads_file, materials_file, outputsDir, analysis_type=analysis_type, flag3d=flag3d)

m.wipe()
m.run_model()

# Time it
get_time(start_time)

# Wipe the model
m.wipe()
