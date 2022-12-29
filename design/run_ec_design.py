from design.eurocodeDesign import EurocodeDesign
from utils.utils import export_to
from pathlib import Path

"""
1. read_input
    Inputs set manually within Input object
    Ec based on uncracked section properties

2. get_design_gravity_loads
    Gets factored design loads and masses
    
3. read_hazard
    Reads seismic hazard
    
4. get_preliminary_solutions
    Gets preliminary solutions from input files
    
5. run_modal_analysis
    Runs modal analysis with 50% reduction of stiffness
    Based on preliminary solutions
    
6. apply_ec_based_analysis
    Uses 50% stiffness reduction
    Based on preliminary solutions
    
    a. run_elastic_analysis
        Gravity only
    b. run_elastic_analysis
        Gravity + ELFM
    c. postprocess_analysis_results
        For each direction
        Compute critical demands
    d. ensure_symmetry
        For each direction
    e. get_critical_of_both_directions
        Gets critical demands of both directions
    f. apply_capacity_design
        Capacity design requirements on combined demands
    g. design_elements based on final demands for the entire building
        Detailing
    
"""

if __name__ == "__main__":
    from design.model2 import Input

    spectra = {}

    site = "Ancona"
    hazard_path = Path.cwd().parents[0] / "examples/inputs/hazard_ancona.pickle"
    prel_solution_path = Path.cwd().parents[0] / "examples/PBEE/inputs/model2"

    ec = EurocodeDesign("Model2", site, flag3d=False)
    ec.data = Input(flag3d=False)
    ec.get_design_gravity_loads()
    hazard = ec.read_hazard(hazard_path)
    solution = ec.get_preliminary_solutions(prel_solution_path)
    model_periods, modalShape, gamma, mstar = ec.run_modal_analysis(solution)

    hinge_models, details = ec.apply_ec_based_analysis(solution, model_periods, modalShape, hazard)

    spectra[site] = ec.spectra

    store = False
    if store:
        export_to(ec.output_path / f"spectra_{site}", spectra[site], "pickle")

        export_to(ec.output_path / f"details_6", details, "pickle")

        export_to(ec.output_path / f"hinge_models_6", hinge_models["x_seismic"], "csv")

        # for i in hinge_models:
        #     export_to(ec.output_path / f"hinge_models_{site}_{i}", hinge_models[i], "csv")
        # export(ec.data, ec.fstiff, ec.output_path, site)
