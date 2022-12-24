import numpy as np
import pandas as pd

from design.elasticAnalysis import ElasticAnalysis


def run_simple_analysis(direction, solution, yield_sa, sls, data):
    if direction == 0:
        nbays = data.n_bays
    else:
        nbays = len(data.spans_y)

    print("[INITIATE] Starting simplified approximate demand estimation...")
    response = pd.DataFrame({'Mbi': np.zeros(data.nst),
                             'Mci': np.zeros(data.nst),
                             'Mce': np.zeros(data.nst)})

    # gravity acceleration, m/s2
    g = 9.81
    base_shear = yield_sa * solution["Mstar"] * solution["Part Factor"] * g
    masses = data.masses / data.n_seismic
    modes = [sls[str(st + 1)]['phi'] for st in range(data.nst)]
    # lateral forces
    forces = np.zeros(data.nst)
    # shear at each storey level
    shear = np.zeros(data.nst)
    for st in range(data.nst):
        forces[st] = masses[st] * modes[st] * base_shear / sum(map(lambda x, y: x * y, masses, modes))
    for st in range(data.nst):
        shear[st] = sum(fi for fi in forces[st:data.nst])

    # Demands on beams and columns in kNm
    # Assuming contraflexure point at 0.6h for the columns
    for st in range(data.nst):
        if st != data.nst - 1:
            response['Mbi'][st] = 1 / 2 / nbays * data.h[st] / 2 * (shear[st] + shear[st + 1])
        else:
            response['Mbi'][st] = 1 / 2 / nbays * data.h[st] / 2 * shear[st]
        # The following is based on assumption that beam stiffness effects are neglected
        ei_external = solution[f"he{st + 1}"] ** 4
        ei_internal = solution[f"hi{st + 1}"] ** 4
        ei_ratio = ei_internal / ei_external
        ei_total = 2 + ei_ratio * (nbays - 1)
        shear_external = shear[st] / ei_total
        shear_internal = shear_external * ei_ratio
        response['Mci'][st] = 0.6 * data.h[st] * shear_internal
        response['Mce'][st] = 0.6 * data.h[st] * shear_external


def run_opensees_analysis(spans_x, spans_y, heights, cross_sections, concrete_modulus, fstiff, loads, flag3d=False):
    """
    Runs OpenSees analysis
    :param spans_x: list
    :param spans_y: list
    :param cross_sections: dataframe
    :param concrete_modulus: float
    :param fstiff: float
    :param loads: list
    :param flag3d: bool
    :return:
    """
    # call OpenSees object
    op = ElasticAnalysis(spans_x, spans_y, heights, cross_sections, concrete_modulus, fstiff, flag3d=flag3d)

    # create the model
    op.create_model()
    op.define_nodes()
    op.define_geometric_transformations()
    op.create_elements()

    # define masses
    masses = op.define_masses(loads)

    num_modes = 3

    return op.run_modal_analysis(num_modes, masses)
