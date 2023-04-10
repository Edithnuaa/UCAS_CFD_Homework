import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('wave.dat')

print(data)
dim = data.shape[0]
x = np.zeros(dim)
u = np.zeros(dim)
p = np.zeros(dim)
rho = np.zeros(dim)

for i, data_i in enumerate(data):
    x[i] = data[i][0]
    u[i] = data[i][1]
    p[i] = data[i][2]
    rho[i] = data[i][3]

fig, ax = plt.subplots(figsize=(8,5),layout='constrained')
ax.plot(x, u, label='u')
ax.plot(x, p, label='p')    
ax.plot(x, rho, label='$\\rho$')
ax.set_xlabel('$x$')
ax.set_ylabel('values')
ax.set_title('Sod t=0.14s')
ax.legend()
plt.savefig('region.png', dpi=400)