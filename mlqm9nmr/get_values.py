import numpy as np 
import os 
import bz2 

# Index 
def get_index():

    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'data')

    path = os.path.join(data_dir, 'index.txt')
    index = np.loadtxt(path, dtype=int)
    
    return list(index)


# Training descriptors
def get_training_descriptors(descriptor):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'data')

    if descriptor == 'acm':
        path = os.path.join(data_dir, 'aCM_4.dat.bz2')
        with bz2.open(path, "rt") as f:
            di = np.loadtxt(f)

    elif descriptor == 'acm_rbf':
        path = os.path.join(data_dir, 'aCM_RBF_4.dat.bz2')
        with bz2.open(path, "rt") as f:
            di = np.loadtxt(f)

    elif descriptor == 'abob':
        path = os.path.join(data_dir, 'aBoB_4.dat.bz2')
        with bz2.open(path, "rt") as f:
            di = np.loadtxt(f)

    elif descriptor == 'abob_rbf':
        path = os.path.join(data_dir, 'aBoB_RBF_4.dat.bz2')
        with bz2.open(path, "rt") as f:
            di = np.loadtxt(f)

    else:
        raise ValueError(f"Unknown descriptor: {descriptor}")

    return di


# coefficients
def get_coefficient(descriptor):
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, 'data')

    if descriptor == 'acm':
        path = os.path.join(data_dir, 'ci_aCM_4.dat')
        ci = np.loadtxt(path)

    if descriptor == 'acm_rbf':
        path = os.path.join(data_dir, 'ci_aCM_RBF_4.dat')
        ci = np.loadtxt(path)

    if descriptor == 'abob':
        path = os.path.join(data_dir, 'ci_aBoB_4.dat')
        ci = np.loadtxt(path)
    
    if descriptor == 'abob_rbf':
        path = os.path.join(data_dir, 'ci_aBoB_RBF_4.dat')
        ci = np.loadtxt(path)

    return ci
