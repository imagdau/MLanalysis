#!

# ------------------------------------------------------------
#
# Apply dimnesionality reduction method to a high dim dataset
#
# ------------------------------------------------------------

import numpy as np
import numpy as np
import argparse
import toml
from sklearn.decomposition import PCA, KernelPCA
from umap import UMAP

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

    
# --- sampling of dataù
# !!! TODO ?
# this might be useful if we are dealing with
# high amount of data, to ease the fit we could
# sample in some way some data and then .fit()
# only on that and later predict on the complete
# !!!
def shuffle(X, Y=None, n=None):
    l = np.arange(X.shape[0])
    random.shuffle(l)
    if Y is None:
        return X[l[:n], :]
    elif Y is None and n is None:
        return X[l, :]
    elif n is None:
        return X[l, :], Y[l]
    else:
        return X[l[:n], :], Y[l[:n]]

    
# --- fit a dimreduction model
def get_embedding(data, 
                  method, 
                  method_dict,
                  fit=True):
    
    # avaliable general methods
    EMBED_METHODS = dict(
    lpca = PCA,
    kpca = KernelPCA,
    umap = UMAP
    )
    
    if method in EMBED_METHODS.keys():
        embed_model = EMBED_METHODS[method](**method_dict)
    else:
        embed_model = method(**method_dict)
    
    print(f"Data shape: {np.shape(X)}")
    print(f"Method selected: {embed_model}")
    
    if fit:
        return embed_model.fit(data)
    else:
        return embed_model

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main(config):
    
    # --- TOML VARIABLES
    # - syst variable
    systname = config.system['dirpath'] + config.system['name']
    
    # --- assumed to work with the .npy pyhton compressed format
    X = np.load(systname)
    if len(X.shape) > 2:
        X = np.concatenate(X)
    
    # --- set up the model
    embedModel = get_embedding(**config.dimred_method)
    
    # --- fit and transform
    embedModel.fit(X)
    # - transfomr
    Xemb = embedModel.transform(X)
    
    # --- save data
    if isinstance(config.dimred_method['method'], str):
        save_name = config.system['name'].replace('soap_', 
                                                  config.dimred_method['method']+'_')
    else:
        save_name = config.system['name'].replace('soap_', 
                                                  'customEmb_')
        
    np.save(save_name, Xemb)
    if config.dimred_method['method'] == 'lpca':
        var = np.cumsum(dimredModel.explained_variance_ratio_)
        np.savetxt('lpca_variace.dat', var)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="config", type=str,
                        help="config file")
    args = parser.parse_args()
    main(get_config(args.config))
