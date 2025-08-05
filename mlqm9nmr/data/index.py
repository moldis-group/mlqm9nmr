import numpy as np 
from sklearn.model_selection import train_test_split

index = np.arange(831925)

n_train = 100000
n_test  =  50000
di_train,di_test,pi_train,pi_test = train_test_split(index, index, train_size=n_train, test_size=n_test, random_state=42)

np.savetxt('index.txt',di_train,fmt='%d')