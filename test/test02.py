from mlqm9nmr import calc_nmr
from mlqm9nmr import plot_nmr

filename   = 'test.xyz'

# OPTIONS: acm, acm_rbf, abob, abob_rbf
descriptor = 'acm' 

# absolute or relative path.
di_path = '/home/surajit/QM9_134k_descriptor/di_train_100k_qm9_acm4_Cnmr_3_cos2.0.txt' 

cs = calc_nmr(filename,descriptor)
plot_nmr(cs)