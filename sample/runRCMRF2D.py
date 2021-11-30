from pathlib import Path
from main import Main
import timeit


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def get_time(start_time):
    elapsed = timeit.default_timer() - start_time
    print('Running time: ', truncate(elapsed, 2), ' seconds')
    print('Running time: ', truncate(elapsed / float(60), 2), ' minutes')


start_time = timeit.default_timer()

# Directories
input_dir = Path.cwd().parents[0] / ".applications/LOSS Validation Manuscript/loss2d"
materials_file = input_dir / "materials.csv"
loads_file = input_dir / "action.csv"
frame = "frame15_mod"
outputsDir = input_dir / f"Cache/{frame}/RCMRF"
section_file = input_dir / f"Cache/{frame}/hinge_models.csv"
# GM directory
gmdir = Path.cwd() / "sample/groundMotion"
gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt", "GMR_durs.txt"]

# RCMRF inputs
hingeModel = "Hysteretic"
analysis_type = ["TH"]
flag3d = False
periods = [0.74, 0.1]
m = Main(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
         analysis_type=analysis_type, system="Perimeter", hinge_model=hingeModel, flag3d=flag3d,
         periods_ida=periods, max_runs=15, export_at_each_step=True)

m.wipe()
m.run_model()

# Wipe the model
m.wipe()

# Time it
get_time(start_time)
