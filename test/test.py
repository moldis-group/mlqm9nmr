import numpy as np 

filename = 'test.xyz' # File name

descriptor = 'acm(4)' # acm(4), acm_rbf(4), abob(4), abob_rbf(4)

di  = np.loadtxt('/home/surajit/QM9_134k_descriptor/di_train_100k_qm9_conabob4_Cnmr_0.05_0.5_6_0.05_3_cos2.0.txt')
ci  = np.loadtxt('../QM9_134k_NMR/data/ci_files/ci_100000_qm9_conabob4_Cnmr_0.05_0.5_6_0.05_3_cos2.0.txt')
sig = 1695.3555 # con-abob4
tms = 186.9704  # ppm at mPW1PW91/6-311+G(2d,p) @ B3LYP/6–31G(2df,p) level


xyzfile = 'test.xyz'
Zs,Rs   = read_xyz(xyzfile)

dq = con_abob(bob_N,6,Zs,Rs,r,sig_w,func,para,beta) # Target atom: 13C-NMR

cs_prd = [] # chemical shifts

for i in range(np.shape(dq)[0]):
    
    dqt = np.sum(np.abs(dq[i] - di), axis=1)  # L1 distance
    kqt = np.exp(-dqt / sig)
   
    nmr_prd = np.sum(ci * kqt)
    cs_prd.append(tms - nmr_prd)
    print(f'C{i+1:d} {tms - nmr_prd:10.2f} ppm')