import subprocess # to launch executable
from concurrent.futures import ProcessPoolExecutor  # to map the thread pool 
from multiprocessing import cpu_count

import numpy as np
import time
def fdtd(x):
    mse = subprocess.run(["../../hpc_project","1", f"{x}"], capture_output=True, text=True)
    array = np.array([float(value) for value in mse.stdout.split("\n") if value.strip() != ''])
    print('.',end='')
    return array
# warmup
for _ in range(3):
    fdtd(0.3)
start_time = time.time()
fdtd(0.3)
timepertask=time.time() - start_time
#numberofcore=cpu_count()
# ntask*timepertask/cpu_count=totaltime
# number of task for an hour entire
import math

ntask=math.trunc(  (1*60*60)//(timepertask*cpu_count()) )

print(f'max task for one hour : {ntask}')

if __name__ == '__main__':
    dt_range=np.linspace(1e-3,2,ntask)
    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        res = executor.map(fdtd,dt_range)
    res=np.array(list(res))
    print(res[:,-1])
    np.savez_compressed("result.npz",res=res,dt_range=dt_range)



