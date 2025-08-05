import numpy as np 

r_min = 0.5
r_max = 6
dr    = 0.05
r     = np.arange(r_min,r_max+dr,dr)
sig   = 0.05 
beta  = 3
func  = 'cos'
para  = 2.0