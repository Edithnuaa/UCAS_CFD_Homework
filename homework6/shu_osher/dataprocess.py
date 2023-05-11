import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('shu.dat')

dim = data.shape[0]
x = np.zeros(dim)
u = np.zeros(dim)
p = np.zeros(dim)
rho = np.zeros(dim)

for i, data_i in enumerate(data):
    x[i] = data[i][0]
    u[i] = data[i][1]
    rho[i] = data[i][2]
    p[i] = data[i][3]


exact = np.loadtxt('exact.dat')

dim = exact.shape[0]
exactx = np.zeros(dim)
exactu = np.zeros(dim)
exactp = np.zeros(dim)
exactrho = np.zeros(dim)

for i, data_i in enumerate(exact):
    exactx[i] = exact[i][0]
    exactu[i] = exact[i][1]
    exactp[i] = exact[i][2]
    exactrho[i] = exact[i][3]

fig, ax = plt.subplots(figsize=(8,5),layout='constrained')
ax.plot(exactx, exactu, label='u.Exact', linestyle='dashed')
ax.plot(exactx, exactp, label='p.Exact', linestyle='dashed')    
ax.plot(exactx, exactrho, label='$\\rho$.Exact', linestyle='dashed')

ax.plot(x, u, label='u')
ax.plot(x, p, label='p')    
ax.plot(x, rho, label='$\\rho$')

ax.set_xlabel('$x$')
ax.set_ylabel('values')
ax.set_title('Shu Osher t=0.1s')
ax.legend()
plt.savefig('WENO_Shu_Osher.png', dpi=400)