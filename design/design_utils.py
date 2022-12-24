import numpy as np
from pathlib import Path

from utils.utils import create_folder


def _create_path(case):
    directory = Path.cwd().parents[0]
    output_path = directory / f"sample/{case}"
    create_folder(output_path)
    return output_path


def _hazard(coef, TR, beta_al):
    x = np.linspace(0.005, 3.0, 201)
    k0 = coef['k0']
    k1 = coef['k1']
    k2 = coef['k2']

    # Ground shaking MAFE
    H = 1 / TR
    p = 1 / (1 + 2 * k2 * (beta_al ** 2))
    Hs = float(k0) * np.exp(-float(k2) * np.log(x) ** 2 - float(k1) * np.log(x))
    MAF = np.sqrt(p) * k0 ** (1 - p) * Hs ** p * np.exp(0.5 * p * k1 ** 2 * (beta_al ** 2))
    p = 1 / (1 + 2 * k2 * (np.power(beta_al, 2)))
    lambdaLS = np.sqrt(p) * k0 ** (1 - p) * H ** p * np.exp(0.5 * p * np.power(k1, 2) * (np.power(beta_al, 2)))
    PGA = np.exp((-k1 + np.sqrt(k1 ** 2 - 4 * k2 * np.log(lambdaLS / k0))) / 2 / k2)
    return lambdaLS, PGA, MAF, x


def get_critical_designs(hinge_models_x, hinge_models_y):
    """
    Modify hinge elements of analysis seismic columns to the strongest (larger My) from designs of both directions
    :param hinge_models_x:
    :param hinge_models_y:
    :return:
    """
    external_hinges_x = hinge_models_x[(hinge_models_x["Position"] == "external") &
                                       (hinge_models_x["Element"] == "Column")].reset_index()
    external_hinges_y = hinge_models_y[(hinge_models_y["Position"] == "external") &
                                       (hinge_models_y["Element"] == "Column")].reset_index()

    for index, row in external_hinges_x.iterrows():
        my_x = external_hinges_x["m1"].iloc[index]
        my_y = external_hinges_y["m1"].iloc[index]
        idx_x = external_hinges_x["index"].iloc[index]
        idx_y = external_hinges_y["index"].iloc[index]
        bay_n_x = external_hinges_x["Bay"].iloc[index]
        bay_n_y = external_hinges_y["Bay"].iloc[index]

        if my_x >= my_y:
            hinge_models_y.iloc[idx_y] = external_hinges_x.drop(columns=["index"]).iloc[index]
            # Modify corresponding Bay number
            hinge_models_y.at[idx_y, "Bay"] = bay_n_y

        else:
            hinge_models_x.iloc[idx_x] = external_hinges_y.drop(columns=["index"]).iloc[index]
            hinge_models_x.at[idx_x, "Bay"] = bay_n_x

    return hinge_models_x, hinge_models_y

