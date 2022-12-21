"""
Example file to run static analysis on a 3D building model
Exports results to outputs directory
Writes model to a .tcl file
"""

from pathlib import Path
from rcmrf import RCMRF
import pickle

from utils.utils import get_time, get_start_time, createFolder
from postprocess.get_rendering import plot_model
from postprocess.utils_plotter import export_figure, export_figure_basic


# record time to run model
start_time = get_start_time()

# Directories
main_dir = Path.cwd()
input_dir = main_dir / "Inputs"
materials_file = input_dir / "materials.csv"
export_model_figs = True
outputsDir = main_dir / "outputs/static"
loads_file = input_dir / "action.csv"

section_file = input_dir / "hinge_models.pickle"
with open(section_file, "rb") as f:
    section_file = pickle.load(f)

# RCMRF inputs
analysis_type = ["ST"]
flag3d = True

m = RCMRF(section_file, loads_file, materials_file, outputsDir, analysis_type=analysis_type, flag3d=flag3d,
          tcl_filename="model")

m.wipe()
m.run_model()

# Time it
get_time(start_time)

# Plotting
if export_model_figs:
    createFolder(outputsDir / "figs")

    fig, _ = plot_model("nodes")
    # With inkscape in Path environmental variables
    # export_figure(fig, filename=outputsDir / "figs/nodes", filetype="svg")

    export_figure_basic(fig, filename=outputsDir / "figs/nodes", filetype="pdf")

    fig, _ = plot_model("elements")
    # With inkscape in Path environmental variables
    # export_figure(fig, filename=outputsDir / "figs/elements", filetype="svg")

    export_figure_basic(fig, filename=outputsDir / "figs/elements", filetype="pdf")

# Wipe the model
m.wipe()
