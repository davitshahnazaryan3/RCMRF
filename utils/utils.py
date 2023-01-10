import os
import subprocess
import timeit
import json
import pickle
import pandas as pd
import numpy as np


def createFolder(directory):
    """
    Checks whether provided directory exists, if no creates one
    :param directory: str
    :return: None
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)


def check_integrity(system, flag3d, hingeModel):
    """
    Checks whether the input arguments have been supplied properly
    :return: None
    """
    if system not in ('perimeter', 'space'):
        raise ValueError('[EXCEPTION] Wrong system type provided, should be Perimeter or Space')
    if flag3d and hingeModel == "haselton":
        raise ValueError('[EXCEPTION] Currently 3D modelling with Haselton springs not supported!')
    print('[SUCCESS] Integrity of input arguments checked')


def read_text_file(name):
    """
    Reads a text file
    :param name: str                            Name of text file in "*.txt" format
    :return: ndarray                            Read column in np.array format
    """
    return np.loadtxt(name)


def export_figure(figure, **kwargs):
    """
    Saves figure as .emf
    :param figure: fig handle
    :param kwargs: filepath: str                File name, e.g. '*\filename'
    :return: None
    """
    inkscape_path = kwargs.get('inkscape', "C://Program Files//Inkscape//bin//inkscape.exe")
    filepath = kwargs.get('filename', None)
    filetype = kwargs.get('filetype', None)
    if filepath is not None:
        path, filename = os.path.split(filepath)
        filename, extension = os.path.splitext(filename)
        svg_filepath = os.path.join(path, filename + '.svg')
        target_filepath = os.path.join(path, filename + f'.{filetype}')
        figure.savefig(svg_filepath, bbox_inches='tight', format='svg')
        subprocess.call([inkscape_path, svg_filepath, f'--export-{filetype}', target_filepath])
        # os.remove(svg_filepath)


def read_text(name, usecols=None):
    """
    Reads a text file
    :param name: str                            Name of text file in "*.txt" format
    :return: ndarray                            Read column in np.array format
    """
    return np.genfromtxt(name, invalid_raise=False, usecols=usecols)


def get_start_time():
    return timeit.default_timer()


def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier


def get_time(start_time):
    elapsed = timeit.default_timer() - start_time
    print('Running time: ', truncate(elapsed, 2), ' seconds')
    print('Running time: ', truncate(elapsed / float(60), 2), ' minutes')


def export_to(filepath, data, filetype):
    """
    Store results in the database
    :param filepath: str                            Filepath, e.g. "directory/name"
    :param data:                                    Data to be stored
    :param filetype: str                            Filetype, e.g. npy, json, pkl, csv
    :return: None
    """
    if filetype == "npy":
        np.save(f"{filepath}.npy", data)
    elif filetype == "pkl" or filetype == "pickle":
        with open(f"{filepath}.pickle", 'wb') as handle:
            pickle.dump(data, handle)
    elif filetype == "json":
        with open(f"{filepath}.json", "w") as json_file:
            json.dump(data, json_file)
    elif filetype == "csv":
        data.to_csv(f"{filepath}.csv", index=False)


def tuple_to_dict(data):
    out = dict((i, d) for i, d in enumerate(data))
    return out


def export(data, fstiff, path, site):
    """
    export cache to path
    :param data: dict                               IPBSD AnconaSaAvg
    :param fstiff: float                            Stiffness reduction factor
    :param path: str                                Path to directory to export the files
    :param site: str
    :return: None
    """
    # Number of storeys
    nst = data.nst

    # Masses (actual frame mass is exported) (for 2D only)
    masses = np.array(data.masses)

    # Creating a DataFrame for loads
    loads = pd.DataFrame(columns=["Storey", "Pattern", "Load"])

    for st in range(1, nst + 1):
        # Masses (for 2D only)
        loads = loads.append({"Storey": st,
                              "Pattern": "mass",
                              "Load": masses[st - 1] / data.n_seismic}, ignore_index=True)

        # Area loads (for both 2D and 3D)
        q = data.inputs["loads"][st - 1]
        loads = loads.append({"Storey": st,
                              "Pattern": "q",
                              "Load": q}, ignore_index=True)

        q = data.inputs["seismic"][st - 1]
        loads = loads.append({"Storey": st,
                              "Pattern": "seismic",
                              "Load": q}, ignore_index=True)

    # Exporting action for use by a Modeler module
    """
    For a two-way slab assumption, load distribution will not be uniform.
    For now and for simplicity, total load over each directions is computed and then divided by global length to 
    assume a uniform distribution. 
    """
    export_to(path / f"action_{site}", loads, "csv")

    # Materials
    fc = data.fc
    fy = data.fy
    Es = data.elastic_modulus_steel

    # Elastic modulus of uncracked concrete
    Ec = (3320 * np.sqrt(fc) + 6900) * fstiff

    materials = pd.DataFrame({"fc": fc,
                              "fy": fy,
                              "Es": Es,
                              "Ec": Ec}, index=[0])

    # Exporting the materials file for use by a Modeler module
    export_to(path / f"materials_{site}", materials, "csv")


def getIndex(target, data, tol=0.):
    if np.where(data >= target)[0].size == 0:
        return np.nan
    else:
        return np.where(data >= target - tol)[0][0]


def create_folder(directory):
    """
    creates directory
    :param directory: str                   Directory to be created
    :return: None
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print("Error: Creating directory. " + directory)

