#!

# ------------------------------------------------------------
#
# Frame by frame analysis given some descriptor based on quippy
#
# ------------------------------------------------------------

import numpy as np
import argparse
import toml
from ase.io import read, write
from quippy import descriptors
from tqdm import tqdm

# ------------------------------------------------------------
# Functions
# ------------------------------------------------------------

# --- .toml file reader
from types import SimpleNamespace

def get_config(filename):
    with open(filename, 'r') as f:
        if filename.endswith('.toml'):
            _config = toml.load(f)
        else:
            raise WrongConfigFormat(
                f"Format of the '{filename}' is not supported. "
                "Available formats: .toml"
            )
        return SimpleNamespace(**_config)

    
# --- frame by frame decsriptor computer
def ds_frame_by_frame(sys_name, ds_obj, 
                      f_range=1, every=1,
                      debug=False):
    ps_list = list()

    # here the range is set to goes from range(1,f+1)
    # and effectively it will output the right amount of
    # frame - necessary due to f-1 in the index=':' selection
    range_ = np.arange(1,f_range+1,every)
    if debug:
        print(len(range_), range_)
    for f in tqdm(range_, desc='Computing descriptor'):
        snap_tmp = read(sys_name, index=str(f-1)+':'+str(f))
        ps_tmp = ds_obj.calc_descriptor(snap_tmp)
        # print(np.shape(ps_tmp))
        ps_list.append(ps_tmp)
    return np.concatenate(ps_list)
    
# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main(config):
    
    # --- TOML VARIABLES
    # - syst variable
    systname = config.system['dirpath'] + config.system['name']
    
    # Here the numbers starts from 1 and not from 0
    # to get just the first frame n_frame=1
    n_frame = config.system['n_frame']
    skip = config.system['every']
    
    # - descriptor variables
    dscr_dict = config.descriptor
    dscr_str = ''
    for string,value in dscr_dict.items():
        if not isinstance(value, str):
            value = str(value)
        if string == 'type':
            dscr_str += value+' '
        else:
            dscr_str += string+'='+value+' '
    # ---
    
    # --- frame by frame computer
    ds = descriptors.Descriptor(dscr_str)
    
    ds_vec = ds_frame_by_frame(sys_name=systname, 
                               ds_obj=ds, 
                               f_range=n_frame, 
                               every=skip)
    
    output_name = config.descriptor['type']+\
                  '_rcut'+str(int(config.descriptor['cutoff']))+\
                  '_n_Z'+\
                  str(config.descriptor['n_Z'])+\
                  '_Z'+config.descriptor['Z'][1:-1]+\
                  '_Nframe'+str(n_frame)+'every'+str(skip)
    
    if config.system['saveFile']:
        np.save(config.system['outheader']+'_'+output_name, ds_vec)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="config", type=str,
                        help="config file")
    args = parser.parse_args()
    main(get_config(args.config))
