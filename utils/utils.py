import os


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
