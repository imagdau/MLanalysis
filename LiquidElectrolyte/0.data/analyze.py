from ase.io import read, write
from quippy import descriptors
import numpy as np

#data is a link to /data/ibm26/MLdata on womble

traj = read('test_traj_10frames.xyz', index=':1')
print(np.unique(traj[0].numbers))

ds_str = 'soap cutoff=6.0 cutoff_transition_width=1.0  n_max=8 l_max=4 atom_sigma=0.5  n_Z=6 Z={1 3 6 8 9 15} n_species=6 species_Z={1 3 6 8 9 15}'
ds = descriptors.Descriptor(ds_str)
ps = ds.calc(traj)

print(ps[0]['data'].shape)
