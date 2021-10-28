from pathlib import Path
from master import Master
import timeit
import pickle


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def get_time(start_time):
    elapsed = timeit.default_timer() - start_time
    print('Running time: ', truncate(elapsed, 2), ' seconds')
    print('Running time: ', truncate(elapsed / float(60), 2), ' minutes')


start_time = timeit.default_timer()

# Directories
input_dir = Path.cwd().parents[0] / ".applications/Loss Validation Manuscript/space/EC8/case11a"
materials_file = input_dir / "materials.csv"
loads_file = input_dir / "action.csv"

outputsDir = input_dir / "RCMRF"
# section_file_x = input_dir / "Cache/framex/hinge_models.csv"
# section_file_y = input_dir / "Cache/framey/hinge_models.csv"
# section_file_gr = input_dir / "Cache/gravity_hinges.csv"
# Directories
# section_file = {"x": section_file_x, "y": section_file_y, "gravity": section_file_gr}

# section_file = input_dir / "Cache/hinge_models.pickle"
# with open(section_file, "rb") as f:
#     section_file = pickle.load(f)

section_file = input_dir / "hinge_models.pickle"
with open(section_file, "rb") as f:
    section_file = pickle.load(f)

# GM directory
gmdir = Path.cwd() / "sample/groundMotion"
gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt", "GMR_durs.txt"]

# RCMRF inputs
hingeModel = "hysteretic"
system = "space"
analysis_type = ["TH"]
direction = 0
flag3d = True
export_at_each_step = True
period_assignment = {"x": 1, "y": 0}
periods = [0.89, 0.87]
m = Master(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
           analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
           direction=direction, export_at_each_step=export_at_each_step, period_assignment=period_assignment,
           periods_ida=periods, max_runs=15)

m.wipe()
m.run_model()

# Wipe the model
m.wipe()

# Time it
get_time(start_time)
