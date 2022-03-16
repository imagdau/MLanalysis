import numpy as np
import os
import re
import ase
from ase.io import read

# --- Files

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def get_filesEW(folder, extension, verbose=True):
    file_list = list()
    for entry in os.listdir(folder):
        if os.path.isfile(os.path.join(folder, entry)):
            if entry.endswith(extension):
                file_list.append(entry)
    file_list.sort(key=natural_keys)
    if verbose:
        print(f"Files:\n{file_list}")
    return file_list


# --- ASE functions

# simple ase reader for the traj file
def ase_xyz_reader(traj_file,
                   index_slice=slice(None, None, None)):
    """
    Simple translator from a XYZ traj file
    to a ASE atoms traj file.
    :param traj_file:
    :param index_slice:
    :return: ASEatoms traj objectv
    """
    print(f"Done with ase version {ase.__version__}")
    print(f"--- Reading {traj_file}")
    return read(traj_file, index=index_slice, format='xyz')


# adding pbc to the xyz traj read by ase
def add_pbc_to_ase(traj_object,
                   box_file,
                   index_slice=slice(None, None, None),
                   workdir=None):
    """
    Add the CELL and PBC information of the original traj
    to the newly prepped ASE traj. The trajectory passed
    to this function is intended to be already sliced.
    :param traj_object: ASEatoms traj object
    :param box_file: file name of box information
    :param index_slice: file slicer
    :param workdir: additional path for specifying files
    :return: ASEatoms traj object
    """
    if workdir:
        box_file = workdir + box_file
    if isinstance(box_file, str):
        print(f"--- Applying PBC {box_file}\n")
        box_file_sliced = np.loadtxt(box_file)[index_slice]
        for (ts, box) in zip(traj_object, box_file_sliced):
            ts.set_cell([box[0], box[1], box[2]])
            ts.set_pbc([1, 1, 1])
    elif isinstance(box_file, float) or isinstance(box_file, int):
        print(f"--- Applying PBC manually set\n")
        for ts in traj_object:
            ts.set_cell([box_file, box_file, box_file])
            ts.set_pbc([1, 1, 1])
    else:
        print(f"Boxfile variable is {type(box_file)}\n")
        print("to add a PBC correction either specify 'boxfile' as 'int' or 'list'")
        print("No PBC correction added")
    return traj_object

