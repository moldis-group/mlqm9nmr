import numpy as np 
import matplotlib.pyplot as plt 

from mlqm9nmr import read_xyz
from mlqm9nmr import acm,acm_rbf,abob,abob_rbf
from mlqm9nmr import get_training_descriptors
from mlqm9nmr import get_coefficient

# mPW1PW91/6-311+G(2d,p) // B3LYP/6–31G(2df,p) 
tms = 186.9704  

def calc_nmr(xyz_file,descriptor,di_path):

    if descriptor == 'acm':      sig = 2643.9900
    if descriptor == 'acm_rbf':  sig = 3097.5820
    if descriptor == 'abob':     sig =   65.5336
    if descriptor == 'abob_rbf': sig = 1682.5215

    if descriptor == 'acm':      p05,p25,p50,p75,p95 = 1297.20, 1612.01, 1833.97, 2058.76, 2367.35
    if descriptor == 'abob':     p05,p25,p50,p75,p95 =   13.78,   31.68,   45.45,   60.84,   81.97
    if descriptor == 'acm_rbf':  p05,p25,p50,p75,p95 = 1010.39, 1683.87, 2334.47, 3039.01, 3923.13
    if descriptor == 'abob_rbf': p05,p25,p50,p75,p95 =  493.72,  851.42, 1166.87, 1528.98, 1932.51

    if di_path == 'bz2': di = get_training_descriptors(descriptor)
    else: di = np.loadtxt(di_path)
    
    ci = get_coefficient(descriptor)

    zs,rs = read_xyz(xyz_file)

    if descriptor == 'acm' or descriptor == 'abob':
        if len(zs) > 29: raise ValueError('The number of atoms is more than 29.')
    if 6 not in set(zs): raise ValueError('No carbon atom present in the dataset.')

    if descriptor == 'acm':      dq = acm(4,6,zs,rs)
    if descriptor == 'acm_rbf':  dq = acm_rbf(4,6,zs,rs)
    if descriptor == 'abob':     dq = abob(4,6,zs,rs)
    if descriptor == 'abob_rbf': dq = abob_rbf(4,6,zs,rs)

    cs = [] # chemical shifts

    for i in range(np.shape(dq)[0]):
        
        dqt = np.sum(np.abs(dq[i] - di), axis=1)  # L1 distance
        kqt = np.exp(-dqt / sig)
    
        nmr_prd = np.sum(ci * kqt)
        cs.append(tms - nmr_prd)

        mini = np.min(dqt)

        if mini < 1e-8: print('Molecule is present in the training dataset')

        if   mini < p05: print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (<p5)')
        elif mini < p25: print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (<p25)')
        elif mini < p50: print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (<p50)')
        elif mini < p75: print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (<p75)')
        else:            print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (>p75)')

    return cs


def plot_nmr(cs):

    fig,ax=plt.subplots(figsize=(6,2))

    cs_min = 0
    cs_max = 200

    if min(cs) < cs_min: cs_min = min(cs) - 10.0
    if max(cs) > cs_max: cs_max = max(cs) + 10.0

    x = np.linspace(cs_min,cs_max,2**13)
    f = np.zeros(len(x))

    for i in range(len(cs)):
        ax.plot([cs[i], cs[i]], [0, 5], color='black')
    
    ax.set_xlim([cs_max, cs_min])
    ax.set_ylim([0, 6])

    ax.set_xticks([0,25,50,75,100,125,150,175,200])
    ax.set_xticklabels([0,25,50,75,100,125,150,175,200],fontsize=12)
    ax.set_yticks([])

    ax.set_xlabel('ppm',fontsize=12)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.savefig('chemical_shift.pdf',format='pdf',bbox_inches='tight')

    return plt.show()