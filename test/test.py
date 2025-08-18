from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'     # test file name
descriptor = 'abob'         # OPTIONS: acm, acm_rbf, abob, abob_rbf

cs = calc_nmr(filename,descriptor)

plot_nmr(cs)