import numpy as np 

# Discrete aCM descriptor
def acm(cm_N,target_atomic_number,mol,func,para,beta,max_atom):

    n       = len(mol.nuclear_charges)
    bags    = []
    indices = []
    nei_dis = []
    hybrid  = []

    for i in range(n):
        ri = mol.coordinates[i]

        distances = []
        hybd_count = 0 # 1:s 2:sp 3:sp2 4:sp3

        for j in range(n):
            rj = mol.coordinates[j]

            if i != j:
                distance = np.linalg.norm(ri-rj)
                if distance < 1.6: # Bond length C-C 
                    hybd_count += 1
                distances.append((j, distance))
        distances.sort(key=lambda x: x[1], reverse=False)

        index = []
        distx = []
        index.append(i)
        distx.append(0.0)
        for j, distance in distances:
            index.append(j)
            distx.append(distance)

        indices.append(index)
        nei_dis.append(distx)
        hybrid.append(hybd_count)

    
    for k in range(n):
        acm = np.zeros( [max_atom,max_atom] )

        for i1,i in enumerate(indices[k]):
            zi = mol.nuclear_charges[i]
            ri = mol.coordinates[i]

            for j1,j in enumerate(indices[k]):
                zj = mol.nuclear_charges[j]
                rj = mol.coordinates[j]

                rij = np.linalg.norm(ri-rj)

                if i == j:
                    acm[i1][i1] = 0.5*zi**2.4 
                else:
                    acm[i1][j1] = zi*zj/rij * gr(rij,beta) # scaling function: 1/r^{beta}
                    
        vec_acm = []

        for i in range(max_atom):
            for j in range(i,max_atom):
                vec_acm.append(acm[i][j])
        
        bags.append(vec_acm)

    
    descriptor = []

    # Fill zero bag
    zero_bag = []
    for x in range(int(max_atom*(max_atom+1)/2)): zero_bag.append(0.0) 


    for i in range(n):
        if mol.nuclear_charges[i] == target_atomic_number:

            neighbour = np.size(indices[i])
            hy_neighbours = hybrid[i]

            if cm_N == 0: desc_bag = bags[i] # CM0

            if cm_N == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]]  # CM0+CM1
            
            if cm_N == 2: 
                if   neighbour >= 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]]  # CM0+CM1+CM2 
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag
            
            if cm_N == 3: 
                if   neighbour >= 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] # CM0+CM1+CM2+CM3
                elif neighbour == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag                                                             + zero_bag
            
            if cm_N == 4: 
                if   neighbour >= 4: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]]# CM0+CM1+CM2+CM3+CM4
                elif neighbour == 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + zero_bag
                elif neighbour == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag                                                              + zero_bag
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag                                                             + zero_bag                                                              + zero_bag
            
            if cm_N == 5: 
                if   neighbour >= 5: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]] + [x*damp_func(func,para,nei_dis[i][4]) for x in bags[indices[i][4]]] # CM0+CM1+CM2+CM3+CM4+CM5
                elif neighbour == 4: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]] + zero_bag
                elif neighbour == 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + zero_bag            + zero_bag
                elif neighbour == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag            + zero_bag            + zero_bag
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag                                                             + zero_bag            + zero_bag            + zero_bag
            
            if cm_N == 'hy':
                if   hy_neighbours >= 4: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]]# CM0+CM1+CM2+CM3+CM4
                elif hy_neighbours == 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + zero_bag
                elif hy_neighbours == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag            + zero_bag
                elif hy_neighbours == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag            + zero_bag            + zero_bag
            
            descriptor.append(desc_bag)

    return descriptor

# Continuous CM descriptor
def acm_rbf(cm_N,target_atomic_number,mol,r,sig,func,para,beta):

    n_atom = len(mol.nuclear_charges)
    n_r    = len(r)

    descs   = []
    indices = []
    nei_dis = []
    hybrid  = []

    for i in range(n_atom):
        zi = mol.nuclear_charges[i]
        ri = mol.coordinates[i]

        di = gau(r,0,0,1,2,sig,beta)

        distances  = []
        hybd_count = 0 # 1:s 2:sp 3:sp2 4:sp3

        for j in range(n_atom):
            zj = mol.nuclear_charges[j]
            rj = mol.coordinates[j]

            if i != j:
                distance = np.linalg.norm(ri-rj)
                if distance < 1.6: # Bond length C-C 
                    hybd_count += 1
                distances.append((j, distance))
        distances.sort(key = lambda x: x[1])

        index = []
        distx = []
        for j,distance in distances:
            index.append(j)
            distx.append(distance)
            
        for j in index:
            zj = mol.nuclear_charges[j]
            rj = mol.coordinates[j]

            di += gau(r,zi,zj,ri,rj,sig,beta)       

        descs.append(list(di)) 
        indices.append(index)
        nei_dis.append(distx)
        hybrid.append(hybd_count)

    descriptor = []
    
    # Fill zero bag
    zero_bag = []
    for x in range(n_r): zero_bag.append(0.0) # only one bag

    for i in range(n_atom):
        if mol.nuclear_charges[i] == target_atomic_number:

            neighbour = len(indices[i])
            hy_neighbours = hybrid[i]

            if cm_N == 0: desc_bag = descs[i] # BoB0

            if cm_N == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]]  # BoB0+BoB1
            
            if cm_N == 2: 
                if   neighbour >= 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]]  # BoB0+BoB1+BoB2 
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag
            
            if cm_N == 3: 
                if   neighbour >= 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] # BoB0+BoB1+BoB2+BoB3
                elif neighbour == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag                                                             + zero_bag
            
            if cm_N == 4: 
                if   neighbour >= 4: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]]# BoB0+BoB1+BoB2+BoB3+BoB4
                elif neighbour == 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + zero_bag
                elif neighbour == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag                                                              + zero_bag
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag                                                             + zero_bag                                                              + zero_bag
            
            if cm_N == 5: 
                if   neighbour >= 5: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]] + [x*damp_func(func,para,nei_dis[i][4]) for x in descs[indices[i][4]]] # BoB0+BoB1+BoB2+BoB3+BoB4+BoB5
                elif neighbour == 4: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]] + zero_bag
                elif neighbour == 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + zero_bag            + zero_bag
                elif neighbour == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag            + zero_bag            + zero_bag
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag                                                             + zero_bag            + zero_bag            + zero_bag
            
            if cm_N == 'hy':
                if   hy_neighbours >= 4: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]]# BoB0+BoB1+BoB2+BoB3+BoB4
                elif hy_neighbours == 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + zero_bag
                elif hy_neighbours == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag            + zero_bag
                elif hy_neighbours == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag            + zero_bag            + zero_bag
            
            descriptor.append(desc_bag)

    return np.array(descriptor)

# Discrete aBoB descriptor
def abob(bob_N,target_atomic_number,mol,func,para,beta):

    n       = np.size(mol.nuclear_charges)
    bags    = []
    indices = []
    nei_dis = []
    hybrid  = []

    bag_len    = 139 # Total bag length for QM9(5k) dataset

    bag_11_len = 19 
    bag_66_len = 8  
    bag_77_len = 5  
    bag_88_len = 3  
    bag_99_len = 2  
    bag_16_len = 20 
    bag_17_len = 17 
    bag_18_len = 18 
    bag_19_len = 9  
    bag_67_len = 8  
    bag_68_len = 8  
    bag_69_len = 8  
    bag_78_len = 6  
    bag_79_len = 5  
    bag_89_len = 3  

    for i in range(n): 
        zi = mol.nuclear_charges[i]
        ri = mol.coordinates[i]
    
        bag_11 = []   # HH   
        bag_66 = []   # CC     
        bag_77 = []   # NN    
        bag_88 = []   # OO    
        bag_99 = []   # FF    
        bag_16 = []   # HC  
        bag_17 = []   # HO   
        bag_18 = []   # HN  
        bag_19 = []   # HF  
        bag_67 = []   # CO  
        bag_68 = []   # CN  
        bag_69 = []   # CF    
        bag_78 = []   # NO  
        bag_79 = []   # NF    
        bag_89 = []   # OF

        distances = []
        hybd_count = 0 # 1:s 2:sp 3:sp2 4:sp3

        for j in range(n):
            zj = mol.nuclear_charges[j]
            rj = mol.coordinates[j]

            if i != j:
                distance = np.linalg.norm(ri-rj)
                if distance < 1.6: # Bond length C-C 
                    hybd_count += 1
                distances.append((j, distance))
        distances.sort(key = lambda x: x[1])

        index = []
        distx = []
        for j,distance in distances:
            index.append(j)
            distx.append(distance)
        
        for j in index:
            zj = mol.nuclear_charges[j]
            rj = mol.coordinates[j]

            rij = np.linalg.norm(ri-rj)
            val = zi*zj/rij * gr(rij,beta) # scaling function: 1/r^{beta}

            if zi == 1:
                if zj == 1: bag_11.append(val)
                elif zj == 6: bag_16.append(val)
                elif zj == 7: bag_17.append(val)
                elif zj == 8: bag_18.append(val)
                elif zj == 9: bag_19.append(val)

            elif zi == 6: 
                if zj == 1: bag_16.append(val)
                elif zj == 6: bag_66.append(val)
                elif zj == 7: bag_67.append(val)
                elif zj == 8: bag_68.append(val)
                elif zj == 9: bag_69.append(val)

            elif zi == 7:
                if zj == 1: bag_17.append(val)
                elif zj == 6: bag_67.append(val)
                elif zj == 7: bag_77.append(val)
                elif zj == 8: bag_78.append(val)
                elif zj == 9: bag_79.append(val)

            elif zi == 8: 
                if zj == 1: bag_18.append(val)
                elif zj == 6: bag_68.append(val)
                elif zj == 7: bag_78.append(val)
                elif zj == 8: bag_88.append(val)
                elif zj == 9: bag_89.append(val)

            elif zi == 9:
                if zj == 1: bag_19.append(val)
                elif zj == 6: bag_69.append(val)
                elif zj == 7: bag_79.append(val)
                elif zj == 8: bag_89.append(val)
                elif zj == 9: bag_99.append(val)

        bag_11 = sorted(bag_11, reverse=True)
        bag_66 = sorted(bag_66, reverse=True)
        bag_77 = sorted(bag_77, reverse=True)
        bag_88 = sorted(bag_88, reverse=True)
        bag_99 = sorted(bag_99, reverse=True)
        bag_16 = sorted(bag_16, reverse=True)
        bag_17 = sorted(bag_17, reverse=True)
        bag_18 = sorted(bag_18, reverse=True)
        bag_19 = sorted(bag_19, reverse=True)
        bag_67 = sorted(bag_67, reverse=True)
        bag_68 = sorted(bag_68, reverse=True)
        bag_69 = sorted(bag_69, reverse=True)
        bag_78 = sorted(bag_78, reverse=True)
        bag_79 = sorted(bag_79, reverse=True)
        bag_89 = sorted(bag_89, reverse=True)

        for _ in range(bag_11_len - np.size(bag_11)): bag_11.append(0)
        for _ in range(bag_66_len - np.size(bag_66)): bag_66.append(0)
        for _ in range(bag_77_len - np.size(bag_77)): bag_77.append(0)
        for _ in range(bag_88_len - np.size(bag_88)): bag_88.append(0)
        for _ in range(bag_99_len - np.size(bag_99)): bag_99.append(0)
        for _ in range(bag_16_len - np.size(bag_16)): bag_16.append(0)
        for _ in range(bag_17_len - np.size(bag_17)): bag_17.append(0)
        for _ in range(bag_18_len - np.size(bag_18)): bag_18.append(0)
        for _ in range(bag_19_len - np.size(bag_19)): bag_19.append(0)
        for _ in range(bag_67_len - np.size(bag_67)): bag_67.append(0)
        for _ in range(bag_68_len - np.size(bag_68)): bag_68.append(0)
        for _ in range(bag_69_len - np.size(bag_69)): bag_69.append(0)
        for _ in range(bag_78_len - np.size(bag_78)): bag_78.append(0)
        for _ in range(bag_79_len - np.size(bag_79)): bag_79.append(0)
        for _ in range(bag_89_len - np.size(bag_89)): bag_89.append(0)

        temp = (
            bag_11 + bag_66 + bag_77 + bag_88 + bag_99 + 
            bag_16 + bag_17 + bag_18 + bag_19 + bag_67 + 
            bag_68 + bag_69 + bag_78 + bag_79 + bag_89     
        )

        indices.append(index)
        bags.append(temp)
        nei_dis.append(distx)
        hybrid.append(hybd_count)

    descriptor = []

    zero_bag = []
    for x in range(bag_len): zero_bag.append(0.0)


    for i in range(n):

        if mol.nuclear_charges[i] == target_atomic_number:

            neighbour = np.size(indices[i])
            hy_neighbours = hybrid[i]

            if bob_N == 0: desc_bag = bags[i] # BoB0

            if bob_N == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]]  # BoB0+BoB1
            
            if bob_N == 2: 
                if   neighbour >= 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]]  # BoB0+BoB1+BoB2 
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag
            
            if bob_N == 3: 
                if   neighbour >= 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] # BoB0+BoB1+BoB2+BoB3
                elif neighbour == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag                                                             + zero_bag
            
            if bob_N == 4: 
                if   neighbour >= 4: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]]# BoB0+BoB1+BoB2+BoB3+BoB4
                elif neighbour == 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + zero_bag
                elif neighbour == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag                                                              + zero_bag
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag                                                             + zero_bag                                                              + zero_bag
            
            if bob_N == 5: 
                if   neighbour >= 5: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]] + [x*damp_func(func,para,nei_dis[i][4]) for x in bags[indices[i][4]]] # BoB0+BoB1+BoB2+BoB3+BoB4+BoB5
                elif neighbour == 4: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]] + zero_bag
                elif neighbour == 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + zero_bag            + zero_bag
                elif neighbour == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag            + zero_bag            + zero_bag
                elif neighbour == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag                                                             + zero_bag            + zero_bag            + zero_bag
            
            if bob_N == 'hy':
                if   hy_neighbours >= 4: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in bags[indices[i][3]]]# BoB0+BoB1+BoB2+BoB3+BoB4
                elif hy_neighbours == 3: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in bags[indices[i][2]]] + zero_bag
                elif hy_neighbours == 2: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in bags[indices[i][1]]] + zero_bag            + zero_bag
                elif hy_neighbours == 1: desc_bag = bags[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in bags[indices[i][0]]] + zero_bag            + zero_bag            + zero_bag
            
            descriptor.append(desc_bag)

    return descriptor

# Continuous aBoB descriptor
def abob_rbf(bob_N,target_atomic_number,mol,r,sig,func,para,beta):

    n_atom = len(mol.nuclear_charges)
    n_r    = len(r)

    descs   = []
    indices = []
    nei_dis = []
    hybrid  = []

    for i in range(n_atom):
        zi = mol.nuclear_charges[i]
        ri = mol.coordinates[i]

        di_bag_11 = gau(r,0,0,1,2,sig,beta)
        di_bag_66 = gau(r,0,0,1,2,sig,beta)
        di_bag_77 = gau(r,0,0,1,2,sig,beta)
        di_bag_88 = gau(r,0,0,1,2,sig,beta)
        di_bag_99 = gau(r,0,0,1,2,sig,beta)
        di_bag_16 = gau(r,0,0,1,2,sig,beta)
        di_bag_17 = gau(r,0,0,1,2,sig,beta)
        di_bag_18 = gau(r,0,0,1,2,sig,beta)
        di_bag_19 = gau(r,0,0,1,2,sig,beta)
        di_bag_67 = gau(r,0,0,1,2,sig,beta)
        di_bag_68 = gau(r,0,0,1,2,sig,beta)
        di_bag_69 = gau(r,0,0,1,2,sig,beta)
        di_bag_78 = gau(r,0,0,1,2,sig,beta)
        di_bag_79 = gau(r,0,0,1,2,sig,beta)
        di_bag_89 = gau(r,0,0,1,2,sig,beta)

        distances  = []
        hybd_count = 0 # 1:s 2:sp 3:sp2 4:sp3

        for j in range(n_atom):
            zj = mol.nuclear_charges[j]
            rj = mol.coordinates[j]

            if i != j:
                distance = np.linalg.norm(ri-rj)
                if distance < 1.6: # Bond length C-C 
                    hybd_count += 1
                distances.append((j, distance))
        distances.sort(key = lambda x: x[1])

        index = []
        distx = []
        for j,distance in distances:
            index.append(j)
            distx.append(distance)
            
        for j in index:
            zj = mol.nuclear_charges[j]
            rj = mol.coordinates[j]

            bag = f'{zi}_{zj}'

            if   bag == '1_1'                :  di_bag_11 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '6_6'                :  di_bag_66 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '7_7'                :  di_bag_77 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '8_8'                :  di_bag_88 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '9_9'                :  di_bag_99 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '1_6' or bag == '6_1':  di_bag_16 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '1_7' or bag == '7_1':  di_bag_17 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '1_8' or bag == '8_1':  di_bag_18 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '1_9' or bag == '9_1':  di_bag_19 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '6_7' or bag == '7_6':  di_bag_67 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '6_8' or bag == '8_6':  di_bag_68 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '6_9' or bag == '9_6':  di_bag_69 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '7_8' or bag == '8_7':  di_bag_78 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '7_9' or bag == '9_7':  di_bag_79 += gau(r,zi,zj,ri,rj,sig,beta)
            elif bag == '8_9' or bag == '9_8':  di_bag_89 += gau(r,zi,zj,ri,rj,sig,beta)  
            else: assert False, f"{bag} bag is missing"

        di = (
            list(di_bag_11)+list(di_bag_66)+list(di_bag_77)+list(di_bag_88)+list(di_bag_99)+
            list(di_bag_16)+list(di_bag_17)+list(di_bag_18)+list(di_bag_19)+list(di_bag_67)+
            list(di_bag_68)+list(di_bag_69)+list(di_bag_78)+list(di_bag_79)+list(di_bag_89)
        )

        descs.append(list(di)) 
        indices.append(index)
        nei_dis.append(distx)
        hybrid.append(hybd_count)

    descriptor = []

    # Fill zero bag
    zero_bag = []
    for x in range(n_r*15): zero_bag.append(0.0) # 15 bags in total
        
    for i in range(n_atom):
        if mol.nuclear_charges[i] == target_atomic_number:

            neighbour = len(indices[i])
            hy_neighbours = hybrid[i]

            if bob_N == 0: desc_bag = descs[i] # BoB0

            if bob_N == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]]  # BoB0+BoB1
            
            if bob_N == 2: 
                if   neighbour >= 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]]  # BoB0+BoB1+BoB2 
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag
            
            if bob_N == 3: 
                if   neighbour >= 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] # BoB0+BoB1+BoB2+BoB3
                elif neighbour == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag                                                             + zero_bag
            
            if bob_N == 4: 
                if   neighbour >= 4: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]]# BoB0+BoB1+BoB2+BoB3+BoB4
                elif neighbour == 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + zero_bag
                elif neighbour == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag                                                              + zero_bag
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag                                                             + zero_bag                                                              + zero_bag
            
            if bob_N == 5: 
                if   neighbour >= 5: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]] + [x*damp_func(func,para,nei_dis[i][4]) for x in descs[indices[i][4]]] # BoB0+BoB1+BoB2+BoB3+BoB4+BoB5
                elif neighbour == 4: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]] + zero_bag
                elif neighbour == 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + zero_bag            + zero_bag
                elif neighbour == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag            + zero_bag            + zero_bag
                elif neighbour == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag                                                             + zero_bag            + zero_bag            + zero_bag
            
            if bob_N == 'hy':
                if   hy_neighbours >= 4: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + [x*damp_func(func,para,nei_dis[i][3]) for x in descs[indices[i][3]]]# BoB0+BoB1+BoB2+BoB3+BoB4
                elif hy_neighbours == 3: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + [x*damp_func(func,para,nei_dis[i][2]) for x in descs[indices[i][2]]] + zero_bag
                elif hy_neighbours == 2: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + [x*damp_func(func,para,nei_dis[i][1]) for x in descs[indices[i][1]]] + zero_bag            + zero_bag
                elif hy_neighbours == 1: desc_bag = descs[i] + [x*damp_func(func,para,nei_dis[i][0]) for x in descs[indices[i][0]]] + zero_bag            + zero_bag            + zero_bag
            
            descriptor.append(desc_bag)

    return np.array(descriptor)