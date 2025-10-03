from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'drug12_07.xyz'
descriptor = 'abob_rbf'

cs = calc_nmr(filename,descriptor,di_path='bz2')
plot_nmr(cs)