#!

# ------------------------------------------------------------
#
# Descriptor utilities
#
# ------------------------------------------------------------
import numpy as np
import copy
from tqdm import tqdm
from quippy import descriptors

# ------------------------------
#
# Functions
#
# ------------------------------


def get_soap_parameters(param_dict, debug=False):
    dscr_str = ''
    for string,value in param_dict.items():
        if not isinstance(value, str):
            value = str(value)
        if string == 'type':
            dscr_str += value+' '
        else:
            dscr_str += string+'='+value+' '
    if debug:
        print(dscr_str)
    return dscr_str


def params2qs(orig_params, multi, double=False):
    """convert param dicts to quippy string"""
    params = copy.deepcopy(orig_params)
    params["n_species"] = len(params["species_Z"])
    Zstr = "{"
    for Z in params["species_Z"]:
        Zstr += str(Z) + " "
    Zstr = Zstr[:-1]+"}"
    params["species_Z"] = Zstr
    for key, value in multi.items():
        s = "{"
        for i in range(0, params["n_species"]):
            s += str(value) + " "
        params[key] = s[:-1] + "}"
    s = params['type']+" "
    for key, value in params.items():
        if key != "Zs":
            s += key + "=" + str(value) + " "
    return s


def get_turbo_parameters(param_dict):
    """convert to string for quippy turbo_soap"""
    dscr_str = ''
    for Z in param_dict["Zs"]:
        param_dict["central_index"] = param_dict["species_Z"].index(Z)+1
        qs = params2qs(param_dict, param_dict["multi"])
        dscr_str += qs + ' '
    dscr_str = dscr_str[:]
    return dscr_str


def descrSaveName(param_dict,frame_tuple):
    interval_str = '-'.join(map(str,frame_tuple))
    if param_dict['type'] == 'soap':
        output_name = param_dict['type']+\
        '_rcut'+str(param_dict['cutoff'])+\
        '_nmax'+str(param_dict['n_max'])+\
        '_lmax'+str(param_dict['l_max'])+\
        '_n_Z'+str(param_dict['n_Z'])+\
        '_Z'+param_dict['Z'][1:-1]+\
        '_'+interval_str+\
        ''
    elif param_dict['type'] == 'soap_turbo':
        output_name = 'turbo_'+\
        'rcutH'+str(param_dict["rcut_hard"])+\
        '_rcutS'+str(param_dict["rcut_soft"])+\
        '_nmax'+str(param_dict['multi']['alpha_max'])+\
        '_lmax'+str(param_dict['l_max'])+\
        '_n_Z'+str(len(param_dict['Zs']))+\
        '_Z'+str(param_dict['Zs'][0])+\
        '_'+interval_str+\
        ''
    return output_name


# ------------------------------
#
# Classes
#
# ------------------------------


class descriptorInit:
    def __init__(self,param_dict):
        self.mode = param_dict['type']
        self.paramDict = param_dict
        
    def __info__(self):
        return print(self.paramDict)
    
    
class SOAPdescriptor(descriptorInit):
    def __init__(self, param_dict):
        super().__init__(param_dict)
        
        self.param_dict = param_dict
        self.dscr_str = get_soap_parameters(param_dict)
        self.dscr_obj = descriptors.Descriptor(self.dscr_str)
    
    def ChunkEval(self, ase_traj, frame_tuple, chunk):
        b,e,s = frame_tuple
        times = int(e/chunk)
        range_chunks = [np.arange(0+chunk*u,chunk+chunk*u+1,s) for u in range(times)]
        for c,range_ in tqdm(enumerate(range_chunks), desc='Computing descriptor (chunks)'):
            print(f"Chunk ({c+1}): {range_[0]} - {range_[-1]}")
            chunk_tmp = ase_traj[range_[0]:range_[-1]]
            ps_tmp = self.dscr_obj.calc_descriptor(chunk_tmp)
            # saving to file
            frame_tpl_tmp = (range_[0],range_[-1],1)
            file_name_tmp = descrSaveName(param_dict,frame_tpl_tmp)
            np.save(file_name_tmp,ps_tmp)
            
    # def ChunkEvalEmbedded(ase_traj, embedObj, frame_tuple, chunk):
    #     pass
            
    def DirectEval(self, ase_traj):
        return np.array(self.dscr_obj.calc_descriptor(ase_traj))
    
    
class TURBOdescriptor(descriptorInit):
    def __init__(self, param_dict):
        super().__init__(param_dict)
        
        self.param_dict = param_dict
        self.dscr_str = get_turbo_parameters(param_dict)
        self.dscr_obj = descriptors.Descriptor(self.dscr_str)
    
    def ChunkEval(self, ase_traj, frame_tuple, chunk):
        b,e,s = frame_tuple
        times = int(e/chunk)
        range_chunks = [np.arange(0+chunk*u,chunk+chunk*u+1,s) for u in range(times)]
        for c,range_ in tqdm(enumerate(range_chunks), desc='Computing descriptor (chunks)'):
            print(f"Chunk ({c+1}): {range_[0]} - {range_[-1]}")
            chunk_tmp = ase_traj[range_[0]:range_[-1]]
            ps_tmp = self.dscr_obj.calc_descriptor(chunk_tmp)
            # saving to file
            frame_tpl_tmp = (range_[0],range_[-1],1)
            file_name_tmp = descrSaveName(self.param_dict,frame_tpl_tmp)
            np.save(file_name_tmp,ps_tmp)
              
    def DirectEval(self, ase_traj):
        return np.array(self.dscr_obj.calc_descriptor(ase_traj))