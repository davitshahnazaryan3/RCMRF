from pathlib import Path

from rcmrf import RCMRF
import pickle

from utils.utils import get_time, get_start_time
from analysis.multiStripeAnalysis import MultiStripeAnalysis
from postprocess.get_rendering import *
from postprocess.utils_plotter import export_figure

if __name__ == "__main__":
    start_time = get_start_time()
    # Directories
    input_dir = Path.cwd()
    materials_file = input_dir / "RCMRF/materials.csv"
    export_model_figs = False
    export_model_to_tcl = False

    site = "LAquila"
    if site == "Milano":
        seismicity = "low"
        seism = "Low"
    elif site == "Ancona":
        seismicity = "medium"
        seism = "Medium"
    else:
        # LAquila
        seismicity = "high"
        seism = "High"

    outputsDir = input_dir / f"RCMRF/{seism}"

    loads_file = input_dir / f"RCMRF/action_{site.lower()}.csv"

    section_file = outputsDir.parents[0] / f"hinge_models_{site} - Copy.pickle"
    # section_file = Path.cwd().parents[0] / "tests/hinge_temp.pickle"
    with open(section_file, "rb") as f:
        section_file = pickle.load(f)

    # GM directory
    gmdir = Path.cwd().parents[2] / f"ProjectsCurrent/NSPerformanceAssessment/sample/groundMotionMSA/{seism}"
    gmfileNames = ["GMR_H1_names.txt", "GMR_H2_names.txt", "GMR_dts.txt", "GMR_durs.txt"]

    # RCMRF inputs
    hingeModel = "hysteretic"
    system = "space"
    analysis_type = ["PO"]
    direction = 1
    flag3d = True
    export_at_each_step = True
    period_assignment = {"x": 1, "y": 0}
    # periods = [1.25, 1.34]
    periods = [1.21, 1.29]
    analysis_time_step = 0.01

    m = RCMRF(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
              analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
              direction=direction, export_at_each_step=export_at_each_step, period_assignment=period_assignment,
              periods_ida=periods, max_runs=15, tcl_filename=None,
              analysis_time_step=analysis_time_step)

    m.wipe()
    m.run_model()

    if "MSA" in analysis_type:
        period, damping, omegas = m.get_modal_parameters()

        msa = MultiStripeAnalysis(section_file, loads_file, materials_file, gmdir, damping, omegas, outputsDir,
                                  system=system, hingeModel=hingeModel, flag3d=flag3d,
                                  export_at_each_step=export_at_each_step, pflag=True,
                                  analysis_time_step=analysis_time_step, drift_capacity=10., gm_scaling_factor=9.81)

        msa.start_process(m.records)

    # if export_model_figs:
    #     fig, _ = plot_model("nodes")
    #     export_figure(fig, filename=outputsDir / "Models/figs/nodes", filetype="svg")
    #
    #     fig, _ = plot_model("elements")
    #     export_figure(fig, filename=outputsDir / "Models/figs/elements", filetype="svg")

    # Wipe the model
    m.wipe()

    print("[SUCCESS] ANALYSIS 1 COMPLETE!!!")

    # --------------------------------------------------------------------------
    # # Directories
    # input_dir = Path.cwd()
    # materials_file = input_dir / "RCMRF/materials.csv"
    # export_model_figs = False
    # export_model_to_tcl = False
    #
    # site = "LAquila"
    # if site == "Milano":
    #     seismicity = "low"
    #     seism = "Low"
    # elif site == "Ancona":
    #     seismicity = "medium"
    #     seism = "Medium"
    # else:
    #     # LAquila
    #     seismicity = "high"
    #     seism = "High"
    #
    # outputsDir = input_dir / f"RCMRF/{seism}"
    #
    # loads_file = input_dir / f"RCMRF/action_{site.lower()}.csv"
    #
    # section_file = outputsDir.parents[0] / f"hinge_models_{site}.pickle"
    # # section_file = Path.cwd().parents[0] / "tests/hinge_temp.pickle"
    # with open(section_file, "rb") as f:
    #     section_file = pickle.load(f)
    #
    # # GM directory
    # gmdir = Path.cwd().parents[2] / f"ProjectsCurrent/NSPerformanceAssessment/sample/groundMotionMSA/{seism}"
    # gmfileNames = ["GMR_H1_names.txt", "GMR_H2_names.txt", "GMR_dts.txt", "GMR_durs.txt"]
    #
    # # RCMRF inputs
    # hingeModel = "hysteretic"
    # system = "space"
    # analysis_type = ["PO"]
    # direction = 0
    # flag3d = True
    # export_at_each_step = True
    # period_assignment = {"x": 1, "y": 0}
    # periods = [0.96, 1.03]
    # # periods = [1.21, 1.29]
    # analysis_time_step = 0.01
    #
    # m = RCMRF(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
    #           analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
    #           direction=direction, export_at_each_step=export_at_each_step, period_assignment=period_assignment,
    #           periods_ida=periods, max_runs=15, tcl_filename=None,
    #           analysis_time_step=analysis_time_step)
    #
    # m.wipe()
    # m.run_model()
    #
    # if "MSA" in analysis_type:
    #     period, damping, omegas = m.get_modal_parameters()
    #
    #     msa = MultiStripeAnalysis(section_file, loads_file, materials_file, gmdir, damping, omegas, outputsDir,
    #                               system=system, hingeModel=hingeModel, flag3d=flag3d,
    #                               export_at_each_step=export_at_each_step, pflag=True,
    #                               analysis_time_step=analysis_time_step, drift_capacity=10., gm_scaling_factor=9.81)
    #
    #     msa.start_process(m.records)
    #
    # # if export_model_figs:
    # #     fig, _ = plot_model("nodes")
    # #     export_figure(fig, filename=outputsDir / "Models/figs/nodes", filetype="svg")
    # #
    # #     fig, _ = plot_model("elements")
    # #     export_figure(fig, filename=outputsDir / "Models/figs/elements", filetype="svg")
    #
    # # Wipe the model
    # m.wipe()
    #
    # # Time it
    # get_time(start_time)
