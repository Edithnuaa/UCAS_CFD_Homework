import numpy as np
import matplotlib.pyplot as plt 
from numpy import pi

def du_approxi(beta, delta_x, u):
    # 计算差分方程的值
    u_array = np.zeros(6)
    for i, u_i in enumerate(u_array): 
        u_array[i] = np.sin(u + (i - 3) * delta_x)
    return ((-1/12-beta)*u_array[0]+(1/2+5*beta)*u_array[1]+(-3/2-10*beta)*u_array[2]+
            (5/6+10*beta)*u_array[3]+(1/4-5*beta)*u_array[4]+(beta)*u_array[5])/delta_x 

def erro_deltax(delta_x):
    # 计算delta_x对应的误差
    beta = -0.010968
    u = np.arange(0, 2*pi, delta_x)
    dim = u.shape
    du = np.zeros(dim)
    for i, du_i in enumerate(du):
        du[i] = du_approxi(beta, delta_x, u[i])
    return np.sum(abs(np.cos(u) - du))

def main():

    delta_x = np.linspace(0.1, 1, 100)
    erro = np.zeros(100)
    for i, erro_i in enumerate(erro):
        erro[i] = erro_deltax(delta_x[i])
    print(erro)
    fit = np.polyfit(np.log(delta_x), np.log(erro), 1)
    fit_erro = fit[1] + fit[0]*np.log(delta_x)

    fig, ax = plt.subplots(layout='constrained')
    ax.plot(np.log(delta_x), np.log(erro), label='diff')
    ax.plot(np.log(delta_x), fit_erro, label='Fitness')
    ax.annotate('$\ln{(erro)}=3.835\ln{(\Delta{x})}-4.597$', xy=(np.log(0.5), fit[1] + fit[0]*np.log(0.5)), 
                xytext=(-2, -4), arrowprops=dict(facecolor='black', shrink=0.05))
    ax.set_xlabel('$\ln{(\Delta{x})}$')
    ax.set_ylabel('$\ln{(error)}$')
    ax.legend()
    plt.savefig('Precision.png',dpi=500)

main()