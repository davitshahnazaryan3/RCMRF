"""
Performs static pushover analysis for a 2D RC MRF using Hysteretic plastic hinge models.
"""

from pathlib import Path

from analysis.multiStripeAnalysis import MultiStripeAnalysis
from rcmrf import RCMRF

from utils.utils import get_time, get_start_time


if __name__ == "__main__":
    start_time = get_start_time()

    # Directories
    main_dir = Path.cwd()
    input_dir = main_dir / "Inputs"
    materials_file = input_dir / "materials.csv"
    loads_file = input_dir / "action_hyst.csv"
    outputsDir = main_dir / "outputs/hysteretic2d"

    # If csv is provided, variable direction is not mandatory
    # If pickle is provided, section is selected based on direction value (default to 0, i.e. X)
    section_file = input_dir / "hinge_models_hyst.csv"

    # RCMRF inputs
    gmdir = main_dir.parents[2] / "IntegratedSeismicDesign/Work/GMs"
    gmfileNames = ["GMR_names1.txt", "GMR_names2.txt", "GMR_dts.txt"]

    hingeModel = "Hysteretic"
    analysis_type = ["MA"]
    flag3d = False
    analysis_time_step = 0.01
    system = "perimeter"
    max_runs = 15
    export_at_each_step = True
    modal_analysis_path = main_dir / "outputs/hysteretic2d/MA.json"

    m = RCMRF(section_file, loads_file, materials_file, outputsDir, gmdir=gmdir, gmfileNames=gmfileNames,
              analysis_type=analysis_type, system=system, hinge_model=hingeModel, flag3d=flag3d,
              export_at_each_step=export_at_each_step, modal_analysis_path=modal_analysis_path, max_runs=max_runs)

    m.wipe()
    m.run_model()

    # Wipe the model
    m.wipe()

    # Time it
    get_time(start_time)
