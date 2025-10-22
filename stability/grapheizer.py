import numpy as np 
from matplotlib import pyplot as plt 
data=np.load("result.npz")['res']
dt=np.load("result.npz")['dt_range']
accuratevalue=data[0]

MSE=[np.square(np.subtract(accuratevalue,sample)).mean() for sample in data]
print(MSE)
# now its plotting time

