import numpy as np


def compute(data, id=-1):
    arr = np.array(data)

    res = [arr.mean(), arr.var()]

    if id == -1:
        print("Mean : " + str(res[0]))
        print("Variance : " + str(res[1]))
    else:
        print("Mean of data " + str(id) + " : " + str(res[0]))
        print("Variance of data " + str(id) + " : " + str(res[1]))


data1 = [174.761, 157.15, 171.373]
data2 = [461.834, 474.535, 476.77]

compute(data1, 1)
compute(data2, 2)
