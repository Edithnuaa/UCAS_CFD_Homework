import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

def kr(beta, alpha):
    return (-1/12-beta)*np.cos(3*alpha)+(1/2+5*beta)*np.cos(2*alpha)+(-3/2-10*beta)*np.cos(alpha)+ \
        (5/6+10*beta)+(1/4-5*beta)*np.cos(alpha)+beta*np.cos(2*alpha)

def ki(beta, alpha):
    return (1/12+beta)*np.sin(3*alpha)-(1/2+5*beta)*np.sin(2*alpha)+(3/2+10*beta)*np.sin(alpha)+ \
        (1/4-5*beta)*np.sin(alpha)+beta*np.sin(2*alpha)

def objkr(beta, alpha):
    # 仅以kr函数作为目标函数进行优化
    return np.power(kr(beta, alpha), 2)

def objki(beta, alpha):
    # 仅以ki函数作为目标函数进行优化
    return np.power((ki(beta, alpha)-alpha),2)

def objfunc(beta, alpha):
    # 整体优化的目标函数，是kr和ki的函数
    # 目标函数 = kr函数的平方 +（ki函数-alpha）的平方
    return objki(beta, alpha)+objkr(beta, alpha)

def sum_err(beta, obj):
    # 对于一个确定的beta值，返回目标函数各alpha取值的和。
    alpha = np.linspace(0, 2, 100)
    sum_err = 0.0
    for item in alpha:
        sum_err += obj(beta, item)
    return np.sqrt(sum_err)

def opt_algo(obj):
    # 二分法求以obj为目标函数时，使得目标函数最小的beta值
    beta_up = 1
    beta_down = -1
    while abs((beta_up - beta_down)) > 1e-9:
        err_up = sum_err(beta_up, obj)
        err_down = sum_err(beta_down, obj)
        if (err_up < err_down):
            beta_down = beta_up - (beta_up - beta_down) / 2
        else:
            beta_up = beta_down + (beta_up - beta_down) / 2
    return beta_up


def optimization():
    
    optki = opt_algo(objki)
    optkr = opt_algo(objkr)
    print("The value of beta when Ki is optimized = ", optki)
    print("The value of beta when Kr is optimized = ", optkr)

    fig, axs = plt.subplots(2, 2, figsize=(10,8),layout='constrained')
    alpha = np.linspace(0, 3, 100)
    axs[0][0].plot(alpha, ki(optki, alpha), label='Optimized')
    axs[0][0].plot(alpha, alpha, label='Standard')
    axs[0][0].set_xlabel("$ \\alpha $")
    axs[0][0].set_ylabel("$K_i$")
    axs[0][0].set_title("Characteristics of $K_i$ with optimized $K_i$", fontsize=10)
    axs[0][0].text(0.5,2,'$\\beta=-0.01097$')
    axs[0][0].legend()

    axs[0][1].plot(alpha, kr(optki, alpha), label='Optimized')
    axs[0][1].set_xlabel("$\\alpha$")
    axs[0][1].set_ylabel("$K_r$")
    axs[0][1].set_title("Characteristics of $K_r$ with optimized $K_i$", fontsize=10)
    axs[0][1].text(0.5,1,'$\\beta=-0.01097$')
    axs[0][1].legend()

    axs[1][0].plot(alpha, ki(optkr, alpha), label='Optimized')
    axs[1][0].plot(alpha, alpha, label='Standard')
    axs[1][0].set_xlabel("$\\alpha$")
    axs[1][0].set_ylabel("$K_i$")
    axs[1][0].set_title("Characteristics of $K_i$ with optimized $K_r$", fontsize=10)
    axs[1][0].text(0.5,2,'$\\beta=-0.08333$')
    axs[1][0].legend()

    axs[1][1].plot(alpha, kr(optkr, alpha), label='Optimized')
    axs[1][1].set_xlabel("$\\alpha$")
    axs[1][1].set_ylabel("$K_r$")
    axs[1][1].set_title("Characteristics of $K_r$ with optimized $K_r$", fontsize=10)
    axs[1][1].text(0.5,0.5e-8,'$\\beta=-0.08333$')
    axs[1][1].legend()

    plt.savefig("hw3.png",dpi=500) 

optimization()

