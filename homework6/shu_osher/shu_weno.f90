module Shu_WENO
    implicit none
    save
    integer, parameter::cells = 2000
    real, dimension(3,cells+1)::u                                                  ! 守恒变量
    real, dimension(3, cells+1)::h_flux, h_flux_plus, h_flux_minus
    real, parameter::u1 = 2.629, rho1 = 3.857, p1 = 10.333, a_coffe=0.2, omega_coffe=5, u2 = 0, p2 = 1 ! Shu-Osher激波管问题初值
    real, dimension(cells+1)::vel, rho, p, energy, c                               ! 一个时间点上各网格点上的速度、密度、压力、能量、声速
    real, dimension(cells+1)::x                                                    ! 空间离散网格点的物理位置
    real, dimension(3, cells+1)::fplus, fminus                                     ! 使用FVS方法进行分裂时，f(U)被分裂成f+和f-，fplus代表f+，fminus代表f-
    real, dimension(3, cells+1)::dfplus, dfminus
    real, parameter::gama=1.4                                                  ! 空气比热比
    real, parameter::delta_x = 10.0/cells                                       ! 网格点之间的间距
    real, parameter::delta_t = 0.0001                                           ! 三阶Runge-Kutta法使用的时间步长
    real, parameter::time = 1.8                                               ! 仿真时间，达到该时间点时仿真结束
    contains

    subroutine init()
        ! 初始化
        ! 完成t=0时刻所有网格点上的所有物理量的初始化
        implicit none
        integer::i, j 

        do j = 1, cells+1
            x(j) = -5.0 + (j-1)*delta_x
        end do

        do i = 1, cells+1
            if (x(i)<-4) then
                vel(i) = u1
                rho(i) = rho1
                p(i) = p1
                energy(i) = (p(i)/(gama-1)+rho(i)*(vel(i)**2)/2.0)/rho(i)
                c(i) = sqrt(gama*p(i)/rho(i))
                u(1,i) = rho(i)
                u(2,i) = rho(i)*vel(i)
                u(3,i) = (p(i)/(gama-1)+rho(i)*(vel(i)**2)/2.0)
            else
                vel(i) = u2
                rho(i) = 1+a_coffe*sin(omega_coffe*x(i))
                p(i) = p2
                energy(i) = (p(i)/(gama-1)+rho(i)*(vel(i)**2)/2.0)/rho(i)
                c(i) = sqrt(gama*p(i)/rho(i))
                u(1,i) = rho(i)
                u(2,i) = rho(i)*vel(i)
                u(3,i) = (p(i)/(gama-1)+rho(i)*(vel(i)**2)/2.0)
            end if
        end do

    end subroutine


    function fu(u)
        ! 计算f(u)
        ! u不一定指代守恒变量U，可以是RK迭代过程中传进来的中间矢量
        implicit none
        real, dimension(3, cells+1), intent(in)::u
        real, dimension(3, cells+1)::fu
        integer::j
        real, dimension(cells+1)::rho_temp, vel_temp, energy_temp, p_temp

        do j = 1, cells+1
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


    function dfu(u)
        ! 计算f对x的偏导数
        implicit none
        real, dimension(3, cells+1), intent(in)::u
        real, dimension(3, cells+1)::dfu
        integer::i, j

        call lax_friedrichs(u)  ! 使用Lax-Friedrichs方法进行FVS分裂

        call h_flux_calc_plus(fplus)
        call h_flux_calc_minus(fminus)
    
        do i = 1, 3
            do j = 1, cells+1
                if (j>1) then
                    dfplus(i, j) = (h_flux_plus(i, j)-h_flux_plus(i,j-1))/delta_x
                    dfminus(i, j) = (h_flux_minus(i, j)-h_flux_minus(i,j-1))/delta_x
                else if (j==1) then
                    dfplus(i, j) = 0
                    dfminus(i, j) = 0
                end if
            end do
        end do

        dfu = (dfplus + dfminus)
        
    end function

    subroutine h_flux_calc_plus(f_tar)
        implicit none
        integer::i, j
        real, dimension(3, cells+1), intent(in)::f_tar
        real, parameter::epsilon=1e-6, para_p=5      ! WENO参数
        real, dimension(3)::alpha, omega, q, is
        ! 通量系数
        real, parameter, dimension(3)::para_c = [0.1, 0.6, 0.3]
        integer::m

        do i=1,3
            do j=1,cells+1
            if (j>=3 .and. j<=cells+1-2) then

                is(1) = 13.0*((f_tar(i,j-2)-2*f_tar(i,j-1)+f_tar(i,j))**2)/12.0+&
                        ((f_tar(i,j-2)-4*f_tar(i,j-1)+3*f_tar(i,j))**2)/4.0
                is(2) = 13.0*((f_tar(i,j-1)-2*f_tar(i,j)+f_tar(i,j+1))**2)/12.0+&
                        ((f_tar(i,j-1)-f_tar(i,j+1))**2)/4.0
                is(3) = 13.0*((f_tar(i,j)-2*f_tar(i,j+1)+f_tar(i,j+2))**2)/12.0+&
                        ((3*f_tar(i,j)-4*f_tar(i,j+1)+f_tar(i,j+2))**2)/4.0

                do m = 1, 3
                    alpha(m) = para_c(m)/((epsilon+is(m))**para_p)
                end do 

                do m = 1,3
                    omega(m) = alpha(m)/(alpha(1)+alpha(2)+alpha(3))
                end do 

                q(1) = f_tar(i,j-2)/3.0-7.0*f_tar(i,j-1)/6.0+11.0*f_tar(i,j)/6.0
                q(2) = -f_tar(i,j-1)/6.0+5.0*f_tar(i,j)/6.0+f_tar(i,j+1)/3.0
                q(3) = f_tar(i,j)/3.0+5.0*f_tar(i,j+1)/6.0-f_tar(i,j+2)/6.0

                h_flux_plus(i,j) = omega(1)*q(1)+omega(2)*q(2)+omega(3)*q(3)
        
            else if (j == cells+1-1) then

                is(1) = 13.0*((f_tar(i,j-2)-2*f_tar(i,j-1)+f_tar(i,j))**2)/12.0+&
                        ((f_tar(i,j-2)-4*f_tar(i,j-1)+3*f_tar(i,j))**2)/4.0
                is(2) = 13.0*((f_tar(i,j-1)-2*f_tar(i,j)+f_tar(i,j+1))**2)/12.0+&
                        ((f_tar(i,j-1)-f_tar(i,j+1))**2)/4.0
                alpha(1) = para_c(1)/((epsilon+is(1))**para_p)
                alpha(2) = para_c(2)/((epsilon+is(2))**para_p)
                omega(1) = alpha(1)/(alpha(1)+alpha(2))
                omega(2) = alpha(2)/(alpha(1)+alpha(2))
                q(1) = f_tar(i,j-2)/3.0-7.0*f_tar(i,j-1)/6.0+11.0*f_tar(i,j)/6.0
                q(2) = -f_tar(i,j-1)/6.0+5.0*f_tar(i,j)/6.0+f_tar(i,j+1)/3.0
                h_flux_plus(i,j) = omega(1)*q(1)+omega(2)*q(2)
        
            else if (j == cells+1) then
                q(1) = f_tar(i,j-2)/3.0-7.0*f_tar(i,j-1)/6.0+11.0*f_tar(i,j)/6.0
                h_flux_plus(i,j) = q(1)
        
            else if (j==1) then
                q(3) = f_tar(i,j)/3.0+5.0*f_tar(i,j+1)/6.0-f_tar(i,j+2)/6.0 
                h_flux_plus(i,j) = q(3)
        
            else if (j==2) then 
                is(2) = 13.0*((f_tar(i,j-1)-2*f_tar(i,j)+f_tar(i,j+1))**2)/12.0+&
                        ((f_tar(i,j-1)-f_tar(i,j+1))**2)/4.0
                is(3) = 13.0*((f_tar(i,j)-2*f_tar(i,j+1)+f_tar(i,j+2))**2)/12.0+&
                        ((3*f_tar(i,j)-4*f_tar(i,j+1)+f_tar(i,j+2))**2)/4.0
                alpha(2) = para_c(2)/((epsilon+is(2))**para_p)
                alpha(3) = para_c(3)/((epsilon+is(3))**para_p)
                omega(2) = alpha(2)/(alpha(2)+alpha(3))
                omega(3) = alpha(3)/(alpha(2)+alpha(3))
                q(2) = -f_tar(i,j-1)/6.0+5.0*f_tar(i,j)/6.0+f_tar(i,j+1)/3.0
                q(3) = f_tar(i,j)/3.0+5.0*f_tar(i,j+1)/6.0-f_tar(i,j+2)/6.0
                h_flux_plus(i,j) = omega(2)*q(2)+omega(3)*q(3)
            end if
        end do
    end do
        
    end subroutine h_flux_calc_plus

    subroutine h_flux_calc_minus(f_tar)
        implicit none
        integer::i, j
        real, dimension(3, cells+1), intent(in)::f_tar
        real, parameter::epsilon=1e-6, para_p=2      ! WENO参数
        real, dimension(3)::alpha, omega, q, is
        ! 通量系数
        real, parameter, dimension(3)::para_c = [0.3, 0.6, 0.1]
        integer::m

        do i=1,3
            do j=1,cells+1
            if (j>=2 .and. j<=cells+1-3) then

                is(1) = 13.0*((f_tar(i,j+1)-2*f_tar(i,j)+f_tar(i,j-1))**2)/12.0+ &
                        ((3*f_tar(i,j+1)-4*f_tar(i,j)+3*f_tar(i,j-1))**2)/4.0
                is(2) = 13.0*((f_tar(i,j+2)-2*f_tar(i,j+1)+f_tar(i,j))**2)/12.0+&
                        ((f_tar(i,j+2)-f_tar(i,j))**2)/4.0
                is(3) = 13.0*((f_tar(i,j+3)-2*f_tar(i,j+2)+f_tar(i,j+1))**2)/12.0+&
                        ((f_tar(i,j+3)-4*f_tar(i,j+2)+3*f_tar(i,j+1))**2)/4.0

                do m = 1, 3
                    alpha(m) = para_c(m)/((epsilon+is(m))**para_p)
                end do 

                do m = 1,3
                    omega(m) = alpha(m)/(alpha(1)+alpha(2)+alpha(3))
                end do 

                q(1) = -f_tar(i,j-1)/6.0+5.0*f_tar(i,j)/6.0+f_tar(i,j+1)/3.0
                q(2) = f_tar(i,j)/3.0+5.0*f_tar(i,j+1)/6.0-f_tar(i,j+2)/6.0
                q(3) = 11.0*f_tar(i,j+1)/6.0-7.0*f_tar(i,j+2)/6.0+f_tar(i,j+3)/3.0

                h_flux_minus(i,j) = omega(1)*q(1)+omega(2)*q(2)+omega(3)*q(3)
        
            else if (j == cells+1) then
                h_flux_minus(i,j) = f_tar(i,j)
        
            else if (j == cells+1-1) then
                q(1) = -f_tar(i,j-1)/6.0+5.0*f_tar(i,j)/6.0+f_tar(i,j+1)/3.0
                h_flux_minus(i,j) = q(1)

            else if (j == cells+1-2) then
                is(1) = 13.0*((f_tar(i,j+1)-2*f_tar(i,j)+f_tar(i,j-1))**2)/12.0+&
                        ((3*f_tar(i,j+1)-4*f_tar(i,j)+3*f_tar(i,j-1))**2)/4.0
                is(2) = 13.0*((f_tar(i,j+2)-2*f_tar(i,j+1)+f_tar(i,j))**2)/12.0+&
                        ((f_tar(i,j+2)-f_tar(i,j))**2)/4.0
                alpha(1) = para_c(1)/((epsilon+is(1))**para_p)
                alpha(2) = para_c(2)/((epsilon+is(2))**para_p)
                omega(1) = alpha(1)/(alpha(1)+alpha(2))
                omega(2) = alpha(2)/(alpha(1)+alpha(2))
                q(1) = -f_tar(i,j-1)/6.0+5.0*f_tar(i,j)/6.0+f_tar(i,j+1)/3.0
                q(2) = f_tar(i,j)/3.0+5.0*f_tar(i,j+1)/6.0-f_tar(i,j+2)/6.0
                h_flux_minus(i,j) = omega(1)*q(1)+omega(2)*q(2)
        
            else if (j==1) then
                is(2) = 13.0*((f_tar(i,j+2)-2*f_tar(i,j+1)+f_tar(i,j))**2)/12.0+&
                        ((f_tar(i,j+2)-f_tar(i,j))**2)/4.0
                is(3) = 13.0*((f_tar(i,j+3)-2*f_tar(i,j+2)+f_tar(i,j+1))**2)/12.0+&
                        ((f_tar(i,j+3)-4*f_tar(i,j+2)+3*f_tar(i,j+1))**2)/4.0
                alpha(2) = para_c(2)/((epsilon+is(2))**para_p)
                alpha(3) = para_c(3)/((epsilon+is(3))**para_p)
                omega(2) = alpha(2)/(alpha(2)+alpha(3))
                omega(3) = alpha(3)/(alpha(2)+alpha(3))
                q(2) = f_tar(i,j)/3.0+5.0*f_tar(i,j+1)/6.0-f_tar(i,j+2)/6.0
                q(3) = 11.0*f_tar(i,j+1)/6.0-7.0*f_tar(i,j+2)/6.0+f_tar(i,j+3)/3.0
                h_flux_minus(i,j) = omega(2)*q(2)+omega(3)*q(3)
        
            end if
        end do
    end do
        
    end subroutine h_flux_calc_minus

    subroutine lax_friedrichs(u)
        ! 使用L-F方法进行流通量矢量分裂
        implicit none
        real, dimension(3, cells+1), intent(in)::u
        real, dimension(cells+1)::rho_temp, vel_temp, c_temp, energy_temp, p_temp
        real::lambda_star

        integer::j

        do j = 1, cells+1
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

    subroutine rungekutta()
        ! Runge-Kutta迭代与数据保存
        implicit none
        real::t=0.0
        real, dimension(3, cells+1)::rk_u1, rk_u2
        
        integer::j

        call init()

        ! 保存初始化数据，确保初始化正确 
        ! 输出到当前文件夹下的'init.dat'文件中
        open(unit=100, file='init.dat', status='old', &
            action='write')
        200 format (f10.3, 3(5x, f10.6))
        do j = 1, cells+1
            write(100,200) x(j), vel(j), rho(j), p(j)
        end do

        ! 在迭代开始前打开数据存储文件，文件名为‘shu.dat'
        ! 最终数据存储在该文件中
        open(unit=101, file='shu.dat', status='old', &
            action='write')
        
        ! 开始迭代
        do
            rk_u1 = u + delta_t*(-dfu(u))
            rk_u2 = u*3.0/4.0 + (rk_u1 + delta_t*(-dfu(rk_u1)))/4.0
            u = u/3.0 + (2.0*(rk_u2 + delta_t*(-dfu(rk_u2))))/3.0
            ! 每一次时间步长迭代后更新各网格点上的物理量 
            do j = 1, cells+1
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
        do j = 1, cells+1
            write(101, 201) x(j), vel(j), rho(j), p(j)
        end do

    end subroutine

end module Shu_WENO