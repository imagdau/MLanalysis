#!

# ------------------------------------------------------------
#
# Traj utilities
#
# ------------------------------------------------------------
import numpy as np
import math
from ase.io import read, write
from tqdm import tqdm
import anaAtoms as aA

# ------------------------------
#
# Functions
#
# ------------------------------

def heuristic_unwrapping(w,box):
    # init the coordinates
    u = np.empty((len(w),3))
    # set the same starting point
    u[0] = w[0]
    # computing the diff of the w traj
    difw = np.diff(w, axis=0)
    # ---
    # Eq. 1 Heuristic method
    for i,dw in enumerate(difw):
        u[i+1] = w[i+1] - np.floor((w[i+1] - u[i])/box[i+1] + 0.5)*box[i+1]
    return u
    

def displacement_unwrapping(w,box):
    # init the coordinates
    u = np.empty((len(w),3))
    # set the same starting point
    u[0] = w[0]
    # computing the diff of the w traj
    difw = np.diff(w, axis=0)
    # ---
    # Eq 2. Displacemnet method
    for i,dw in enumerate(difw):
        u[i+1] = u[i] + (w[i+1] - w[i]) - np.floor((w[i+1] - w[i])/box[i+1] + 0.5)*box[i+1]
    return u
        
        
def hybrid_unwrapping(w,box):
    # init the coordinates
    u = np.empty((len(w),3))
    # set the same starting point
    u[0] = w[0]
    # computing the diff of the w traj
    difw = np.diff(w, axis=0)
    # ---
    # Eq. 12 Hybrid mehtod
    for i,dw in enumerate(difw):
        u[i+1] = u[i]+(w[i+1]-w[i])-np.floor((w[i+1]-w[i])/box[i+1]+0.5)*box[i+1]-np.floor((w[i]-u[i])/box[i]+0.5)*(box[i+1]-box[i])
    return u


def get_MSD(xyz):
    # assuming shape (Nat,Nframe,XYZ)
    msd = np.empty(np.shape(xyz))
    xyz0 = [t[0] for t in xyz]
    for i,t0 in enumerate(xyz0):
        msd[i] = np.array([np.linalg.norm(t-t0,axis=1)**2 for t in xyz[i]])
    return msd


# ------------------------------
#
# Classes
#
# ------------------------------


class trajectoryHandler:
    def __init__(self, dirname, sysname, read_frame_tuple=None):
        self.dirname = dirname
        self.sysname = sysname
        self.readFrames = read_frame_tuple

class TrajLoader(trajectoryHandler):
    def __init__(self, dirname, sysname, read_frame_tuple, 
                 traj_species_dict, zshift_tuple, unwrap_dict):
        super().__init__(dirname, sysname, read_frame_tuple)

        print('\n# --- Init the trajectory')
        # just read the first frame to get the quantities
        frameZero = read(self.dirname+self.sysname, index='0:1')
        self.RcutCorrectionDict = traj_species_dict['rcut_correction']
        self.UnwrapDict = unwrap_dict
        
        # extraction of the CM position for each "molecule"
        print('\n# --- Extracting molecular information')
        frameZeroMol = aA.extract_molecs(frameZero[:1], fct=self.RcutCorrectionDict)
        self.atMolID = frameZero[0].arrays['molID']
        self.molSym = frameZeroMol[0].arrays['molSym']
        self.Znumbers = frameZero[0].numbers
        # VVB thing
        self.__OldZnumbers = frameZero[0].numbers
        
        # shifting mumbojumbo
        if zshift_tuple:
            self._mol_shifted = zshift_tuple[0]
            self._Z_shifted = zshift_tuple[1]
            shift_ = np.max(self.Znumbers)

            # loop inside the moltypes and select the right type
            for idx, m in enumerate(self.molSym):
                if m == self._mol_shifted:
                    # create a mask
                    mask = self.atMolID == idx
                    for i,ture in enumerate(mask):
                        if ture:
                            if self.Znumbers[i] in self._Z_shifted:
                                self.Znumbers[i] += shift_
                            else:
                                pass
                            
            
    def readFrame(self,n_frame,Zdiff=True):
        ase_frame = read(self.dirname+self.sysname, index=f'{n_frame}:{n_frame+1}')[0]
        if Zdiff and hasattr(self, 'Znumbers'):
            ase_frame.numbers = self.Znumbers
        return ase_frame
    
    
    def readTraj(self,Zdiff=True):
        print(f"\n--- Loading traj {self.readFrames} ---\n")
        b,e,s = self.readFrames
        ase_traj = read(self.dirname+self.sysname, index=f'{b}:{e}:{s}')
        
        if self.UnwrapDict:
            _ = self.trajUnwrapper(frame_tuple=self.readFrames, 
                                   traj=ase_traj, 
                                   **self.UnwrapDict)
        
        if Zdiff and hasattr(self, 'Znumbers'):
            for snap in ase_traj:
                snap.numbers = self.Znumbers
        return ase_traj
    
    
    def trajUnwrapper(self,frame_tuple,species,method='hybrid',traj=None):
        print(f"\n--- Unwrapping the trajs {species} (COM) ---\n")
        UNWRAP_FUNC = dict(
            heuristic = heuristic_unwrapping,
            displacement = displacement_unwrapping,
            hybrid = hybrid_unwrapping,
        )
        
        mol_traj = aA.extract_molecs_molID(traj, molID=self.atMolID, fct=self.RcutCorrectionDict)

        # assuming cubic
        box_val = [atom.cell[0][0] for atom in mol_traj]
        # ref species for the unwrap
        ref_spec_idx = [[idx for idx,spec in enumerate(self.molSym) if spec == MOL] for MOL in species]
        
        unwrapped_coords = dict()
        for s,idxmask in enumerate(ref_spec_idx):
            
            unwrapped_coords[species[s]] = list()
            for idx in tqdm(idxmask, desc=f'Unwrapping: {species[s]}'):
                
                w = np.array([ts[idx].position for ts in mol_traj])
                unwrapped_coords[species[s]].append(UNWRAP_FUNC[method](w,box_val))
        
        # saving the trajs
        interval_str = '-'.join(map(str,frame_tuple))
        for key,item in unwrapped_coords.items():
            save_name_tmp = method+'_unwrap_traj_'+key+'_'+interval_str
            np.save(save_name_tmp, item)
        
        return unwrapped_coords
        

    # def save_traj(self,saveName):
    #     write()
            
