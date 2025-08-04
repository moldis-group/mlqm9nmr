import numpy as np 
import matplotlib.pyplot as plt 

from qm9nmr-ml import acm,acm_rbf,abob,abob_rbf
from qm9nmr-ml import gau_cs

def plot_nmr(cs):

    fig,ax=plt.subplots(figsize=(8,6))

    cs_min = 0
    cs_max = 175

    if min(cs_prd) < cs_min: cs_min = min(cs_prd)
    if max(cs_prd) > cs_max: cs_max = max(cs_prd)

    x = np.linspace(cs_min,cs_max,2**13)
    f = np.zeros(len(r))

    for cs in cs_prd:
        ax.scatter(cs,-0.1,marker='x',color='red')
        f += gau_cs(r,cs)

    for i in range(int(max(gr))+1):
        ax.plot([cs_min,cs_max],[i+1,i+1],color='grey',alpha=0.1,linestyle='--')

    ax.plot(x,f,color='black')

    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(top=False, right=False, left=False)
    ax.tick_params(axis='x', labelsize=16)
    ax.set_yticklabels([])

    ax.set_xlim([cs_min,cs_max])
    ax.invert_xaxis()

    return plt.show()