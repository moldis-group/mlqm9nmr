import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

from mlqm9nmr import read_xyz
from mlqm9nmr import acm,acm_rbf,abob,abob_rbf
from mlqm9nmr import gau_cs
from mlqm9nmr import get_coefficient

# mPW1PW91/6-311+G(2d,p) @ B3LYP/6–31G(2df,p) 
tms = 186.9704  

def calc_nmr(xyz_file,di_path,descriptor):

    if descriptor == 'acm':      sig = 2643.9900
    if descriptor == 'acm_rbf':  sig = 3381.9080
    if descriptor == 'abob':     sig = 65.5336
    if descriptor == 'abob_rbf': sig = 1695.3555

    if descriptor == 'acm':      p25,p50,p75,p90 = 1617.04, 1828.13, 2066.72, 2248.00
    if descriptor == 'abob':     p25,p50,p75,p90 =   31.35,   45.40,   61.37,   75.66
    if descriptor == 'acm_rbf':  p25,p50,p75,p90 = 1682.29, 2360.82, 3075.64, 3549.91
    if descriptor == 'abob_rbf': p25,p50,p75,p90 =  855.39, 1167.80, 1547.30, 1813.36

    di = np.loadtxt(di_path)
    ci = get_coefficient(descriptor)

    zs,rs = read_xyz(xyz_file)

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

        df = pd.DataFrame({'dij': dqt})
        p25q = df['dij'].quantile(0.25)

        if p25q < p25: print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (Good)')
        else:          print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm (Moderate)')

    return cs


def plot_nmr(cs):

    fig,ax=plt.subplots(figsize=(8,6))

    cs_min = 0
    cs_max = 175

    if min(cs) < cs_min: cs_min = min(cs)
    if max(cs) > cs_max: cs_max = max(cs)

    x = np.linspace(cs_min,cs_max,2**13)
    f = np.zeros(len(x))

    for cs in cs:
        ax.scatter(cs,-0.1,marker='x',color='red')
        f += gau_cs(x,cs)

    for i in range(int(max(f))+1):
        ax.plot([cs_min,cs_max],[i+1,i+1],color='grey',alpha=0.1,linestyle='--')

    ax.plot(x,f,color='black')

    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(top=False, right=False, left=False)
    ax.tick_params(axis='x', labelsize=16)
    ax.set_yticklabels([])

    ax.set_xlim([cs_min,cs_max])
    ax.invert_xaxis()

    plt.savefig('chemical_shift.pdf',format='pdf',bbox_inches='tight')

    return plt.show()