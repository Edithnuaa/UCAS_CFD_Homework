import numpy as np 
import matplotlib.pyplot as plt


filename = 'func_wave.dat'

data = np.loadtxt(filename)
time = np.zeros(data.shape[0])
func = np.zeros(data.shape[0])

for i, data_i in enumerate(data):
    time[i] = data_i[0]
    func[i] = data_i[1]

fig, ax = plt.subplots(figsize=(8,5), layout='constrained')
ax.plot(time, func)
ax.set_xlabel('$p^*$')
ax.set_ylabel('$F(p^*)$')
ax.set_title('$Function ~ F(p^*)$')
plt.savefig('Fpstar.png', dpi=400)