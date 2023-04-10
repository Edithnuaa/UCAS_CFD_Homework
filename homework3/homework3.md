### 第一题

构造差分格式并进行Fourier分析：

![1680120569219](image/homework3/1680120569219.png)

![1680120584148](image/homework3/1680120584148.png)

根据Taylor公式展开后求出的差分格式，推断出该格式具有四阶精度。

根据Fourier分析结果，将 $a_6 = \frac{\beta}{\Delta x}$ 作为自由参数，其中 $\beta$ 是需要进行优化的参数。

分别对 $K_i$ 和 $K_r$ 进行优化，目标函数分别选为 $\sum (K_i - \alpha)$ 和 $\sum(K_r)$ 。使用二分法求出使得目标函数最小化的最优解。代码文件：optimization.py。结果如下图：

![1680121036807](image/homework3/1680121036807.png)

以下分析均选取 $\beta=-0.01097$ 为参数值。



下面使用 $ u=sin(x) $ 作为测试函数，测量该差分格式的精度。代码文件：precisoin.py。结果如下图：

![1680121176781](image/homework3/1680121176781.png)

拟合后的直线斜率为3.835，接近4，与Taylor展开推导得到的四阶精度近似。



下面进行数值验证，选取20个网格点，使用Runge-Kutta法进行显式时间推进，使用上述差分格式进行空间推进，代码文件：validatoin.py，当 $t=20s$ 时，L2误差为0.014070756311295866，结果如下图：

![1680121464210](image/homework3/1680121464210.png)

当 $t=50s$ 时，L2误差为：0.05480585000830214，结果如下图：

![1680121491122](image/homework3/1680121491122.png)

L2误差关于时间的关系如下图：

![1680123549675](image/homework3/1680123549675.png)

L2误差始终保持在很小的范围内。

### 第二题

![1680123869058](image/homework3/1680123869058.png)

### 代码附录

#### optimizatoin.py

```python
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


```

#### precision.py

```python
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
```


#### validation.py

```python
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

```
