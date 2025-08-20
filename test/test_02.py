from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'     # test file name
descriptor = 'abob'         # OPTIONS: 'acm', 'acm_rbf', 'abob', or 'abob_rbf'
path       = 'aBoB_4.dat'   # absolute or relative path for the di file

cs = calc_nmr(filename,descriptor,di_path=path)

plot_nmr(cs)