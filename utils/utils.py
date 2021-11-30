import os
import timeit
import json
import pickle
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
