from pathlib import Path
from main import Main
import timeit
import pickle
from visualize.get_rendering import *
from visualize.utils_plotter import export_figure


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def get_time(start_time):
    elapsed = timeit.default_timer() - start_time
    print('Running time: ', truncate(elapsed, 2), ' seconds')
    print('Running time: ', truncate(elapsed / float(60), 2), ' minutes')


start_time = timeit.default_timer()

# Directories
input_dir = Path.cwd()
materials_file = input_dir / "materials.csv"
export_model_figs = False
export_model_to_tcl = False

outputsDir = input_dir / "RCMRF"

site = "Ancona"
if site == "Milano":
    seismicity = "low"
elif site == "Ancona":
    seismicity = "medium"
else:
    seismicity = "high"

loads_file = input_dir / f"action_{site.lower()}.csv"

section_file = outputsDir / f"hinge_models_{site}.pickle"
# section_file = Path.cwd().parents[0] / "tests/hinge_temp.pickle"
with open(section_file, "rb") as f:
    section_file = pickle.load(f)

# GM directory
gmdir = Path.cwd() / "sample/groundMotion"
gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt", "GMR_durs.txt"]

# RCMRF inputs
hingeModel = "hysteretic"
system = "space"
analysis_type = ["MA"]
direction = 0
flag3d = True
export_at_each_step = True
period_assignment = {"x": 1, "y": 0}
periods = [0.95, 0.90]
m = Main(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
         analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
         direction=direction, export_at_each_step=export_at_each_step, period_assignment=period_assignment,
         periods_ida=periods, max_runs=15, tcl_filename=f"model_{seismicity}",
         export_model_to_tcl=export_model_to_tcl)

m.wipe()
m.run_model()

if export_model_figs:
    fig, _ = plot_model("nodes")
    export_figure(fig, filename=outputsDir / "Models/figs/nodes", filetype="svg")

    fig, _ = plot_model("elements")
    export_figure(fig, filename=outputsDir / "Models/figs/elements", filetype="svg")

# Wipe the model
m.wipe()

# Time it
get_time(start_time)
