import numpy as np
import matplotlib.pyplot as plt
from numpy import pi

def LU(delta_x, u_array):
    # 计算空间分布各网格点的LU，等于du的负数。
    beta = np.array([-0.016601569950581],dtype='float64')           # beta参数
    u_copy = np.concatenate(([0], u_array, u_array[1:], u_array[1:]))     # 列出三个波长上的网格点
    dim = u_array.shape[0]                                   # 一个波长上网格点的个数
    u_out = np.zeros(dim, dtype='float64')                   # 输出值的存储向量
    for i in np.arange(0, dim):                              # 循环求出个网格点上的du值
        u_out[i] = ((-1/12-beta)*u_copy[i-3+dim]+(1/2+5*beta)*u_copy[i-2+dim]+(-3/2-10*beta)*u_copy[i-1+dim]+
            (np.round(5/6, decimals=16)+10*beta)*u_copy[i+dim]+(1/4-5*beta)*u_copy[i+1+dim]+beta*u_copy[i+2+dim])/delta_x 
    return -u_out

def RKmethod(delta_x, u_n, delta_t):
    # 3步Runge-Kutta方法
    u1 = u_n + delta_t * LU(delta_x, u_n)
    u2 = 0.75 * u_n + (u1 + delta_t * LU(delta_x, u1)) * 0.25 
    u_n1 = u_n / 3.0 + (2.0 * (u2 + delta_t * LU(delta_x, u2))) / 3.0
    return u_n1

def validation():
    x = np.linspace(0,2*pi,20) # 精确点
    
    nodes = np.linspace(0, 2*pi, 20, dtype='float64') # 取网格点
    u = np.sin(nodes,dtype='float64')   # 初始值
    delta_x = nodes[1] - nodes[0]
    L2erro = np.zeros(5001)
    time = 50 
    t = 0.0
    delta_t = np.array([0.01], dtype='float64')             # 时间步长
    filename = open('log.txt', mode='w') # 打开文件夹，用于存储每一次时间步长下的u
    while t < time:
        print('Time = ', t, 'u = ', u, file=filename) # 存储时间t和u
        u = RKmethod(delta_x, u, delta_t)          # 进行Runge-Kutta迭代
        t += delta_t
        y = np.sin(x-t) # 精确值
        L2erro[int(t/0.01)] = np.sqrt(np.sum(delta_x*np.square(y-u)))
    print('Time = ', t, 'u = ', u, file=filename)     # 存储最后一次计算的u
    filename.close()

    # 画图
    # 注释部分为数值验证图
#    fig, ax = plt.subplots(layout='constrained')
#    ax.plot(x, y, label='Exact.solution')
#    ax.scatter(nodes, u, color='b')
#    ax.plot(nodes, u, label='Num.solution')
#   ax.set_xlabel('$x$')
 #   ax.set_ylabel('$u$')
 #   ax.set_title('$t = 50$')
#    ax.legend()
 #   plt.savefig('validation50.png')


    # 以下为L2误差图
    fig2, ax2 = plt.subplots(layout='constrained')
    t4e = np.linspace(0,50,5001,dtype='float')
    ax2.scatter(t4e, L2erro, label='L2error')
    ax2.set_xlabel('$t$')
    ax2.set_ylabel("$L2error$")
    ax2.legend()
    ax2.set_title('L2 error')
    plt.savefig('L2.png')
validation()
