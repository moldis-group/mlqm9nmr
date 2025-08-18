from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'     # test file name
descriptor = 'abob'         # OPTIONS: acm, acm_rbf, abob, abob_rbf
di_path    = 'di_abob4.txt' # absolute or relative path 

cs = calc_nmr(filename,di_path,descriptor)
plot_nmr(cs)