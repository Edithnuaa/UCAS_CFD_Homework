module Sod_FVS_NND
    implicit none
    save
    integer, parameter::dim = 200
    real, dimension(3,dim)::u                                                  ! 守恒变量
    real::u1 = 0, rho1 = 1, p1 = 1, u2 = 0, rho2 = 0.125, p2 = 0.1             ! sod激波管问题初值
    real, dimension(dim)::vel, rho, p, energy, c                               ! 一个时间点上各网格点上的速度、密度、压力、能量、声速
    real, dimension(dim)::x                                                    ! 空间离散网格点的物理位置
    real, parameter::gama=1.4                                                  ! 空气比热比
    real, parameter::delta_x = 2.0/199.0                                       ! 网格点之间的间距
    real, parameter::delta_t = 0.001                                           ! 三阶Runge-Kutta法使用的时间步长
    real, parameter::time = 0.14                                               ! 仿真时间，达到该时间点时仿真结束
    real, dimension(3, dim)::fplus, fminus                                     ! 使用FVS方法进行分裂时，f(U)被分裂成f+和f-，fplus代表f+，fminus代表f-
    real, dimension(3, dim)::partial_fx                                        ! 使用差分格式计算的f(U)对x的偏导数
    integer, parameter::flag = 2                                               ! 方法选择参数，1:使用steger-warming，2:使用Lax-Friedrichs
    contains

    subroutine init()
        ! 初始化
        ! 完成t=0时刻所有网格点上的所有物理量的初始化
        implicit none
        integer::i, j 

        do j = 1, dim
            x(j) = -1.0 + (j-1)*delta_x
        end do


        do i = 1, dim
            if (x(i)<0) then
                vel(i) = u1
                rho(i) = rho1
                p(i) = p1
                energy(i) = (p1/(gama-1)+rho1*(u1**2)/2.0)/rho(i)
                c(i) = sqrt(gama*p1/rho1)
                u(1,i) = rho1
                u(2,i) = rho1*u1
                u(3,i) = (p1/(gama-1)+rho1*(u1**2)/2.0)
            else
                vel(i) = u2
                rho(i) = rho2
                p(i) = p2
                energy(i) = (p2/(gama-1)+rho2*(u2**2)/2.0)/rho(i)
                c(i) = sqrt(gama*p2/rho2)
                u(1,i) = rho2
                u(2,i) = rho2*u2
                u(3,i) = (p2/(gama-1)+rho2*(u2**2)/2.0)
            end if
        end do

    end subroutine


    function fu(u)
        ! 计算f(u)
        ! u不一定指代守恒变量U，可以是RK迭代过程中传进来的中间矢量
        implicit none
        real, dimension(3, dim), intent(in)::u
        real, dimension(3, dim)::fu
        integer::j
        real, dimension(dim)::rho_temp, vel_temp, energy_temp, p_temp

        do j = 1, dim
            ! 反解对应的基础物理量，密度、速度、能量、压强
            rho_temp(j) = u(1, j)
            vel_temp(j) = u(2, j)/rho_temp(j)
            energy_temp(j) = u(3, j)/rho_temp(j)
            p_temp(j) = (u(3,j)-(rho_temp(j)*vel_temp(j)**2)/2)*(gama-1)
            ! 根据反解的物理量，求解f(U)
            fu(1, j) = rho_temp(j)*vel_temp(j)
            fu(2, j) = rho_temp(j)*(vel_temp(j)**2)+p_temp(j)
            fu(3, j) = rho_temp(j)*energy_temp(j)*vel_temp(j)+p_temp(j)*vel_temp(j)
        end do

    end function


    function dfu(u, flag)
        ! 计算f对x的偏导数
        implicit none
        real, dimension(3, dim)::dfu
        real, dimension(3, dim), intent(in)::u
        integer, intent(in)::flag

        if (flag == 1) then
            call steger_warming(u)  ! 使用Steger-warming方法进行FVS分裂
            call nnd                ! 使用NND格式对fplus和fminus的偏导数进行计算，并求和得出f(U)对x的偏导数
            dfu = partial_fx / delta_x   ! 除以网格间距
        else if (flag == 2) then
            call lax_friedrichs(u)  ! 使用Lax-Friedrichs方法进行FVS分裂
            call nnd                ! 使用NND格式对fplus和fminus的偏导数进行计算，并求和得出f(U)对x的偏导数
            dfu = partial_fx / delta_x   ! 除以网格间距
        end if
        
    end function

    subroutine steger_warming(u)
        ! 使用steger_warming方法进行流通量矢量分裂
        implicit none
        real, dimension(3, dim), intent(in)::u
        real, dimension(dim)::rho_temp, vel_temp, c_temp, energy_temp, p_temp
        integer::i, j
        real, dimension(3, dim)::lambda, lambda_plus, lambda_minus
        real, parameter::epsilon=0.1
        real::w_minus, w_plus

        do j = 1, dim
            ! 反解各网格点上的物理量
            rho_temp(j) = u(1, j)
            vel_temp(j) = u(2, j)/rho_temp(j)
            energy_temp(j) = u(3, j)/rho_temp(j)
            p_temp(j) = (u(3,j)-(rho_temp(j)*vel_temp(j)**2)/2)*(gama-1)
            c_temp(j) = sqrt(gama*p_temp(j)/rho_temp(j))
        end do

        do j = 1, 200
            lambda(1, j) = vel_temp(j)
            lambda(2, j) = vel_temp(j) - c_temp(j)
            lambda(3, j) = vel_temp(j) + c_temp(j)
        end do

        do j = 1, 200
            do i = 1, 3
            ! 求特征值+和特征值-
            lambda_plus(i, j) = (lambda(i, j) + sqrt(lambda(i, j)**2+epsilon**2))/2
            lambda_minus(i, j) = (lambda(i, j) - sqrt(lambda(i, j)**2+epsilon**2))/2
            end do
        end do

        do j = 1, 200
            ! 求分裂后的fplus和fminus（过程不同）
            w_plus = ((3-gama)*(lambda_plus(2,j)+lambda_plus(3,j))*(c_temp(j)**2))/(2*(gama-1))
            w_minus = ((3-gama)*(lambda_minus(2,j)+lambda_minus(3,j))*(c_temp(j)**2))/(2*(gama-1))
            fplus(1, j) = (2*(gama-1)*lambda_plus(1, j)+lambda_plus(2, j)+lambda_plus(3, j))*rho_temp(j)/(2*gama)
            fplus(2, j) = ((2*(gama-1)*lambda_plus(1,j)*vel_temp(j)+lambda_plus(2,j)*(vel_temp(j)- &
                c_temp(j))+lambda_plus(3,j)*(vel_temp(j)+c_temp(j)))*rho_temp(j))/(2*gama)
            fplus(3, j) = ((gama-1)*lambda_plus(1,j)*vel_temp(j)**2+(lambda_plus(2,j)*(vel_temp(j)- &
                c_temp(j))**2)/2+(lambda_plus(3,j)*(vel_temp(j)+c_temp(j))**2)/2+w_plus)*rho_temp(j)/(2*gama)
            fminus(1, j) = (2*(gama-1)*lambda_minus(1,j)+lambda_minus(2,j)+lambda_minus(3,j))*rho_temp(j)/(2*gama)
            fminus(2, j) = (2*(gama-1)*lambda_minus(1,j)*vel_temp(j)+lambda_minus(2,j)*(vel_temp(j)- &
                c_temp(j)+lambda_minus(3,j)*(vel_temp(j)+c_temp(j))))*rho_temp(j)/(2*gama)
            fminus(3, j) = (((gama-1)*lambda_minus(1,j)*vel_temp(j)**2+(lambda_minus(2,j)*(vel_temp(j)- &
                c_temp(j))**2)/2+(lambda_minus(3,j)*(vel_temp(j)+c_temp(j))**2)/2+w_minus)*rho_temp(j))/(2*gama)
        end do

    end subroutine

    subroutine lax_friedrichs(u)
        ! 使用L-F方法进行流通量矢量分裂
        implicit none
        real, dimension(3, dim), intent(in)::u
        real, dimension(dim)::rho_temp, vel_temp, c_temp, energy_temp, p_temp
        real::lambda_star

        integer::j

        do j = 1, dim
            ! 反解各网格点上的物理量
            rho_temp(j) = u(1, j)
            vel_temp(j) = u(2, j)/rho_temp(j)
            energy_temp(j) = u(3, j)/rho_temp(j)
            p_temp(j) = (u(3,j)-(rho_temp(j)*vel_temp(j)**2)/2)*(gama-1)
            c_temp(j) = sqrt(gama*p_temp(j)/rho_temp(j))
        end do

        ! 全局Lax-Friedrichs分裂，求出该时间点上最大的lambda*
        lambda_star = maxval(abs(vel_temp)+c_temp)

        fplus = (fu(u)+lambda_star*u)/2
        fminus = (fu(u)-lambda_star*u)/2

    end subroutine

    subroutine nnd()
        !  对fplus和fminus使用NND差分格式
        implicit none

        integer::i, j

        do i = 1, 3
            do j = 3, 198
                partial_fx(i, j) = (fplus(i, j)+(min_mod(fplus(i, j)-fplus(i, j-1), fplus(i, j+1)-fplus(i, j)))/2 + & 
                fminus(i, j+1)-(min_mod(fminus(i, j+1)-fminus(i, j), fminus(i, j+2)-fminus(i, j+1)))/2)- &
                (fplus(i, j-1)+(min_mod(fplus(i, j-1)-fplus(i, j-2), fplus(i, j)-fplus(i, j-1)))/2 + &
                fminus(i, j)-(min_mod(fminus(i,j)-fminus(i, j-1), fminus(i, j+1)-fminus(i, j)))/2)
            end do
        end do

    end subroutine

    function min_mod(val1, val2)
        ! 计算min mod函数
        ! val1和val2符号相同，取绝对值小的
        ! val1和val2符号相反，取0
        implicit none
        real, intent(in)::val1, val2
        real::min_mod

        if ((val1/val2)<0) then
            min_mod = 0
        else
            if (abs(val1) < abs(val2)) then
                min_mod = val1
            else
                min_mod = val2
            end if 
        end if
    end function

    subroutine rungekutta()
        ! Runge-Kutta迭代与数据保存
        implicit none
        real::t=0.0
        real, dimension(3, dim)::rk_u1, rk_u2
        
        integer::j

        call init()

        ! 保存初始化数据，确保初始化正确 
        ! 输出到当前文件夹下的'init.dat'文件中
        open(unit=100, file='init.dat', status='old', &
            action='write')
        200 format (f10.3, 3(5x, f10.6))
        do j = 1, 200
            write(100,200) x(j), vel(j), rho(j), p(j)
        end do

        ! 在迭代开始前打开数据存储文件，文件名为‘sod.dat'
        ! 最终数据存储在该文件中
        open(unit=101, file='sod.dat', status='old', &
            action='write')
        
        ! 开始迭代
        do
            rk_u1 = u + delta_t*(-dfu(u, flag))
            rk_u2 = u*3.0/4.0 + (rk_u1 + delta_t*(-dfu(rk_u1, flag)))/4.0
            u = u/3.0 + (2.0*(rk_u2 + delta_t*(-dfu(rk_u2, flag))))/3.0
            ! 每一次时间步长迭代后更新各网格点上的物理量 
            do j = 1, 200
                rho(j) = u(1, j)
                vel(j) = u(2, j)/rho(j)
                energy(j) = u(3, j)/rho(j)
                p(j) = (u(3,j)-(rho(j)*vel(j)**2)/2)*(gama-1)
                c(j) = sqrt(gama*p(j)/rho(j))
            end do
            t = t + delta_t
            if (t>=time) exit
        end do

        ! 输出对应time时刻各网格点上的物理量
        201 format (f10.3, 3(5x, f15.6))
        do j = 1, 200
            write(101, 201) x(j), vel(j), rho(j), p(j)
        end do

    end subroutine

end module Sod_FVS_NND