from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'     # test file name
descriptor = 'abob'         # OPTIONS: 'acm', 'acm_rbf', 'abob', or 'abob_rbf'

cs = calc_nmr(filename,descriptor,di_path='bz2')

plot_nmr(cs)