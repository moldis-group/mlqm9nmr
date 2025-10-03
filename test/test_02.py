from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'drug12_07.xyz'
descriptor = 'abob_rbf'
path       = 'aBoB_RBF_4.dat'

cs = calc_nmr(filename,descriptor,di_path=path)
plot_nmr(cs)