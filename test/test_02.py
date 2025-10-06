from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'
descriptor = 'acm_rbf_4'
path       = 'aCM_RBF_4.dat'

cs = calc_nmr(filename,descriptor,di_path=path)
plot_nmr(cs)