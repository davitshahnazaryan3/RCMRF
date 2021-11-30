from pathlib import Path

import multiprocessing as mp

from main import Main
import timeit
import pickle

from utils.utils import get_time, get_start_time
from analysis.multiStripeAnalysis import MultiStripeAnalysis
from visualize.get_rendering import *
from visualize.utils_plotter import export_figure


if __name__ == "__main__":
    start_time = get_start_time()
    # Directories
    input_dir = Path.cwd()
    materials_file = input_dir / "materials.csv"
    export_model_figs = False
    export_model_to_tcl = False

    outputsDir = input_dir / "RCMRF/Medium"

    site = "Ancona"
    if site == "Milano":
        seismicity = "low"
    elif site == "Ancona":
        seismicity = "medium"
    else:
        seismicity = "high"

    loads_file = input_dir / f"action_{site.lower()}.csv"

    section_file = outputsDir.parents[0] / f"hinge_models_{site}.pickle"
    # section_file = Path.cwd().parents[0] / "tests/hinge_temp.pickle"
    with open(section_file, "rb") as f:
        section_file = pickle.load(f)

    # GM directory
    gmdir = Path.cwd() / "groundMotionMSA/Medium"
    gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt", "GMR_durs.txt"]

    # RCMRF inputs
    hingeModel = "hysteretic"
    system = "space"
    analysis_type = ["MSA"]
    direction = 0
    flag3d = True
    export_at_each_step = True
    period_assignment = {"x": 1, "y": 0}
    periods = [0.95, 0.90]
    analysis_time_step = 0.01

    m = Main(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
             analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
             direction=direction, export_at_each_step=export_at_each_step, period_assignment=period_assignment,
             periods_ida=periods, max_runs=15, tcl_filename=f"model_{seismicity}",
             export_model_to_tcl=export_model_to_tcl, analysis_time_step=analysis_time_step)

    m.wipe()
    m.run_model()

    if "MSA" in analysis_type:
        period, damping, omegas = m.get_modal_parameters()

        msa = MultiStripeAnalysis(section_file, loads_file, materials_file, gmdir, damping, omegas, outputsDir,
                                  system=system, hingeModel=hingeModel, flag3d=flag3d,
                                  export_at_each_step=export_at_each_step, pflag=True,
                                  analysis_time_step=analysis_time_step, drift_capacity=10.)

        msa.start_process(m.records)

    # if export_model_figs:
    #     fig, _ = plot_model("nodes")
    #     export_figure(fig, filename=outputsDir / "Models/figs/nodes", filetype="svg")
    #
    #     fig, _ = plot_model("elements")
    #     export_figure(fig, filename=outputsDir / "Models/figs/elements", filetype="svg")

    # Wipe the model
    m.wipe()

    # Time it
    get_time(start_time)
