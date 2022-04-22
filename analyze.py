from ase.io import read, write
from quippy import descriptors
import numpy as np

#data is a link to /data/ibm26/MLdata on womble

traj = read('data/LiquidElectrolyte/traj_2.1.extxyz', ':2')
print(np.unique(traj[0].numbers))

ds_str = 'soap cutoff=6.0 cutoff_transition_width=1.0  n_max=8 l_max=4 atom_sigma=0.5  n_Z=6 Z={1 3 6 8 9 15} n_species=6 species_Z={1 3 6 8 9 15}'
ds = descriptors.Descriptor(ds_str)
ps = ds.calc(traj)

print(ps[0]['data'])

import sys
sys.path.append('../../software/python-atoms/')
import anaAtoms as aa

moltraj = aa.extract_molecs(traj, fct={'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1})
menvs = aa.mol_envs(moltraj, ['EC', 'EMC', 'Li', 'PF6'], Rcut=6.0, returnEnvs=True)
moltraj[0].arrays['molEnv']

idxs = traj[0].arrays['molID']
mtype = moltraj[0].arrays['molSym']

for idx, m in enumerate(mtype):
    if m=='EC':
        traj[0].numbers[idxs==idx]+=15
