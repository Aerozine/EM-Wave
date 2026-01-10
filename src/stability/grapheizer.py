import numpy as np 
from matplotlib import pyplot as plt 
data=np.load("result.npz")['res']
dt=np.load("result.npz")['dt_range']
idx = np.argmax(dt > 0.2)   # first index where dt > 0.5
accuratevalue = data[idx]

MSE=np.array([np.square(np.subtract(accuratevalue,sample)).mean() for sample in data])
print(MSE)
# now its plotting time
mask = np.isfinite(MSE) & (MSE > 0)
MSE_plot = np.clip(MSE[mask], 1e-100, 1e20)
plt.plot(dt[mask], MSE_plot)
plt.yscale('log')

plt.xlabel(r"Courant safety factor $S_c$")  # Label x-axis with units
plt.ylabel("Mean Squared Error over time of a point \n at coord (42,42) clipped up to $10^{20}$")  # Label y-axis
plt.title("Convergence of FDTD Solution according to Courant safety factor")  # Add title
plt.grid(True, which='both', ls='--', alpha=0.5)  # Grid for easier reading
plt.tight_layout()  # Adjust layout
plt.savefig("MSE_plot.svg", format='svg', dpi=300, transparent=True)

