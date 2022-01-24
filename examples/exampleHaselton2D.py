"""
Performs static pushover analysis for a 2D RC MRF using Haselton plastic hinge models.
"""

from pathlib import Path
from rcmrf import RCMRF

from utils.utils import get_time, get_start_time


start_time = get_start_time()

# Directories
main_dir = Path.cwd()
input_dir = main_dir / "Inputs"
materials_file = input_dir / "materials.csv"
loads_file = input_dir / "action.csv"
outputsDir = main_dir / "outputs/haselton"

# If csv is provided, variable direction is not mandatory
# If pickle is provided, section is selected based on direction value (default to 0, i.e. X)
section_file = input_dir / "hinge_models2D.csv"

# RCMRF inputs
hingeModel = "Haselton"
analysis_type = ["PO"]
flag3d = False

m = RCMRF(section_file, loads_file, materials_file, outputsDir, analysis_type=analysis_type, hinge_model=hingeModel,
          flag3d=flag3d)

m.wipe()
m.run_model()

# Wipe the model
m.wipe()

# Time it
get_time(start_time)
