from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'
descriptor = 'acm_rbf_4'

cs = calc_nmr(filename,descriptor,di_path='bz2')
plot_nmr(cs)