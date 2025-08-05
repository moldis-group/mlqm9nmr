import numpy as np 

# Reading xyz file
def read_xyz(filename):

    with open(filename,'r') as file:
        num = file.readline()
        nam = file.readline()
        xyz = file.readlines()
    file.close()

    Zs = []
    Rs = []

    for i in range(int(num)):
        a,x,y,z=xyz[i].split()
        if   a.lower() == 'h': atno = 1
        elif a.lower() == 'c': atno = 6
        elif a.lower() == 'n': atno = 7
        elif a.lower() == 'o': atno = 8
        elif a.lower() == 'f': atno = 9
        else: raise ValueError(f'Molecule has {a} atom, which is not included in the training.')
        
        Zs.append(atno)
        Rs.append(np.array([float(x),float(y),float(z)]))

    return np.array(Zs), np.array(Rs)


# Normalized Gaussian function
def ngau_func(r,zi,zj,ri,rj,sig,beta):
    gr = 1/r**beta
    N = 1/(sig*np.sqrt(2*np.pi))
    rij = np.linalg.norm(ri-rj)
    return N * (zi*zj/rij) * np.exp(- ((r-rij)**2) / (2*sig**2) ) * gr


# Scaling function
def scal_func(r,beta):
    return 1/r**beta


# Damping function
def damp_func(func,para,d):
    if func == 'exp' : return np.exp(-para*d)
    if func == 'pol' : return 1./(1+d)**para
    if func == 'cos' : return 0.5 * ( 1 + np.cos(np.pi*d/para) )


# Gaussian for assigning chemical shifts
def gau_cs(r,cs):
    sw = 0.5
    return np.exp(- ((r-cs)**2) / (2*sw**2) ) 
