import numpy as np
import string
from multiprocessing import Pool
import argparse
import numpy as np
import copy

def params2qs(orig_params, multi, double=False):
   """convert param dicts to quippy string"""
   #update species from Zs
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
   s = "soap_turbo "
   for key, value in params.items():
      if key != "Zs":
         s += key + "=" + str(value) + " "
   if double:
      s = s.replace("{", "{{")
      s = s.replace("}", "}}")
   return s

def params2fullqs(params, multi):
   """convert to string for gap_fit"""
   gap_string = ""
   for Z in params["Zs"]:
      params["central_index"] = params['species_Z'].index(Z)+1
      qs = params2qs(params, multi, double=True)
      gap_string += qs + " : \n"
   gap_string = gap_string[:-2]
   return gap_string

params = {}
params["Zs"] = [1,3]
params["species_Z"] = [1,3,6,8,9,15]
params["l_max"] = 4
params["n_species"] = 6
params["rcut_hard"] = 6
params["rcut_soft"] = 5.0
params["basis"] = "poly3gauss"
params["scaling_mode"] = "polynomial"
params["radial_enhancement"] = 1
multi = {
  "alpha_max":8,
  "atom_sigma_r":0.2,
  "atom_sigma_t":0.2,
  "atom_sigma_r_scaling":0.1,
  "atom_sigma_t_scaling":0.1,
  "amplitude_scaling":1.0,
  "central_weight":1.0
}

print(params2fullqs(params, multi))
