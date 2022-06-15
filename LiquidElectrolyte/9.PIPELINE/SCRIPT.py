#
import numpy as np
import time
import pipeTrajectory as pT
import pipeDescriptor as pD

# -----------------------------------------------------------------
#
#   INPUT DICTs
#
# -----------------------------------------------------------------

# TRAJ
traj_dict = dict(
    dirname='/mnt/c/Users/andre/Documents/Work/1.Cambridge/0.systems/1.ioanSystems/data/LiquidElectrolyte/',
    #'/data/ibm26/MLdata/LiquidElectrolyte/',
    sysname='traj_2.1.xyz',
    #'traj_2.1.xyz',
    read_frame_tuple = (0,500,1),
    
    traj_species_dict = dict(
        rcut_correction = {'H':1,'C':1,'O':1,'Li':0.1,'P':1,'F':1},
        molecular_species = ['EC', 'EMC', 'Li', 'PF6']),
    
    zshift_tuple = ('EC', [6, 8]),
    
    unwrap_dict = dict(
        species = ['Li','PF6'],
        method = 'hybrid'),

)

# DESCR (from JD routine)
turbo_param_multi_dict = dict(
    alpha_max = 8,
    atom_sigma_r = 0.2,
    atom_sigma_t = 0.2,
    atom_sigma_r_scaling = 0.1,
    atom_sigma_t_scaling = 0.1,
    amplitude_scaling = 1.0,
    central_weight = 1.0,
)

turbo_param_dict = dict(
    type = "soap_turbo",
    Zs = [3],
    species_Z = [1,3,6,8,9,15,21,23],
    l_max = 0,
    n_species = 8,
    rcut_hard = 4.5,
    rcut_soft = 3.5,
    basis = "poly3gauss",
    scaling_mode = "polynomial",
    radial_enhancement = 1,
    multi = turbo_param_multi_dict,
)

descr_dict = dict(
    chunkEval = 500,
    save_file = True,
)

# -----------------------------------------------------------------
#
#   MAIN PIPE
#
# -----------------------------------------------------------------

# --- 1. Traj

# - Init the analysis
trajObj = pT.TrajLoader(**traj_dict)

# - Reading the amounts of frames needed
t0 = time.time()
traj_read = trajObj.readTraj()
t1 = time.time()
print(f"{np.round(t1-t0, 1)}s")


# --- 2. Descriptor

print(f"\n# --- Computing the descriptor ---\n")

# - Computing the descriptor
t0 = time.time()
descrObj = pD.TURBOdescriptor(turbo_param_dict)
print(f"Paramenters:\n{descrObj.dscr_str}")
soap_descr = descrObj.ChunkEval(traj_read, frame_tuple=trajObj.readFrames, chunk=descr_dict['chunkEval'])
print(f"{np.round(t1-t0, 1)}s")
print(f"\nData matrix: {soap_descr.shape}")


# --- 3. Low Dim Embedding
# !!! TODO
    
print("\n### END ###\n")