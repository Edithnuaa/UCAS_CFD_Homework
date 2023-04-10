    module RH_relation
        implicit none
        save
        real::p1=1, rho1=1, u1=0, p2=0.1, rho2=0.125, u2=0     ! p1，rho1，u1，左侧的压力，密度，速度；p2，rho2，u2，右侧的压力，密度，速度 
        real, parameter::gama = 1.4                            ! 气体比热比
        real::z_left = 0, z_left_head = 0, z_left_tail = 0                ! 左侧激波速度、膨胀波波前速度、膨胀波波后速度
        real::z_right = 0, z_right_head = 0, z_right_tail = 0              ! 右侧激波速度、膨胀波波前速度、膨胀波波后速度
        real::p_star = 0, u_star = 0, rho_star_left = 0, rho_star_right = 0    ! 中心区的压力、速度和中心区间断接触两侧的密度
        integer::Lflag, Rflag, midflag = 0                     ! flag = 1, 代表激波；flag = 0, 代表膨胀波
        real::c_star_left = 0, c_star_right = 0                        ! 中心区接触间断两侧的声速
        integer::i
        real, dimension(201)::x = [(-1.0 + i*0.01, i = 0, 200)]! x坐标数组
        real, dimension(201)::u = 0, p = 0, rho = 0                    ! x数组内的速度、压力和密度
        real::t=0.14                                           ! 时间
        real::rarefaction_left_down = 0, rarefaction_left_up = 0       ! 左侧稀疏波的区域上限和下限
        real::rarefaction_right_down = 0, rarefaction_right_up = 0     ! 右侧稀疏波的区域上线和下限
        real::shock_left = 0, shock_right = 0                          ! 左侧和右侧激波的位置
        real::mid_wave = 0                                         ! 接触间断的位置
        real::c1,c2


        contains 

        function f_wave(p, rho, p_star)
            ! 波速方程
            implicit none 
            real,intent(in)::p, rho, p_star
            real::f_wave
            if (p_star >= p) then ! p_star >= p
                f_wave = (p_star - p)/(rho*sqrt(gama*p/rho)*sqrt((((gama+1)*p_star)/(2*gama*p))+((gama-1)/(2*gama))))
            else if (p_star < p) then ! p_star < p  
                f_wave = ((2*sqrt(gama*p/rho)/(gama-1)))*((p_star/p)**((gama-1)/(2*gama))-1)
            end if
        return
        end function f_wave

        function df_wave(p, rho, p_star)
            ! 波速方程的导数
            implicit none
            real, intent(in)::p, rho, p_star
            real::df_wave
            real::df_mid
        
            df_mid = rho*sqrt(gama*p/rho)*sqrt((((gama+1)*p_star)/(2*gama*p))+((gama-1)/(2*gama)))
        
            if (p_star >= p) then
                df_wave = 1/(df_mid) + ((p_star-p)*rho*sqrt(gama*p/rho)*(gama+1))/ & 
                    (4*gama*p*sqrt((((gama+1)*p_star)/(2*gama*p))+((gama-1)/(2*gama)))*(df_mid**2)) 
            else if (p_star < p) then
                df_wave = (2*sqrt(gama*p/rho)*(gama-1)*(p_star/p)**((-1-gama)/(2*gama)))/((gama-1)*2*gama*p)
            end if
        end function df_wave

        function func_wave(p_star)
            ! 计算F(p_star)函数
            ! p1- 压力1，p2- 压力2，rho1- 密度1，rho2- 密度2，p_star- p*
            implicit none
            real::func_wave
            real, intent(in)::p_star
            func_wave = f_wave(p1, rho1, p_star) + f_wave(p2, rho2, p_star)
        end function func_wave

        subroutine func_wave_draw()
            ! 计算并保存func_wave函数的数据
            implicit none
            real::p_star = 0.0
            real::f_temp
            integer::func_wave_stat
            character(len=80)::func_wave_msg
            open(unit=10, file='func_wave.dat', status='old', &
                action='write', iostat=func_wave_stat, iomsg=func_wave_msg)
            101 format (f10.2, 5x, f10.6)
            do
                f_temp = func_wave(p_star)
                write(10, 101) p_star, f_temp
                p_star = p_star + 0.01
                if (abs(p_star-2.0)<1e-4) exit 
            end do 
        end subroutine func_wave_draw

        subroutine Newton2p_star(p_star)
            ! 牛顿法求F(p_star)的零点
            implicit none
            real, intent(out)::p_star
            real::p_init = 1.0                                     ! 迭代初始值
            do
                p_star = p_init - (f_wave(p1, rho1, p_init)+f_wave(p2, rho2, p_init)+u1+u2)/ &
                        (df_wave(p1,rho1,p_init)+df_wave(p2,rho2,p_init))
                if (abs(p_star-p_init)<0.00001) exit        
                p_init = p_star
            end do
        
        end subroutine Newton2p_star

        subroutine wave_type()
            ! 判断两侧波的类型：激波 or 膨胀波
            implicit none
            real::F_0, F_p1, F_p2

            F_p1 = func_wave(p1)
            F_p2 = func_wave(p2)    
            F_0 = func_wave(0.0)

            if ((u1-u2) >= F_p1) then
                write(*,*) "两侧都为激波"
                Lflag = 1
                Rflag = 1
            else if (((u1-u2) < F_p1) .and. ((u1-u2) >= F_p2)) then
                write(*,*) "左侧波为膨胀波，右侧为激波"
                Lflag = 0
                Rflag = 1
            else if (((u1-u2) < F_p2) .and. ((u1-u2) >= F_0)) then
                write(*,*) "两侧均为膨胀波" 
                Lflag = 0
                Rflag = 0
            else if ((u1-u2) < F_0) then
                write(*,*) "两侧为膨胀波，中间为真空区"
                Lflag = 0
                Rflag = 0
                midflag = 1
            else
                write(*,*) "Invalid"
            end if

        end subroutine wave_type

        subroutine wave()
            ! 计算两侧波的速度、密度
            implicit none
            real::a_left, a_right
            c1 = sqrt(gama*p1/rho1)
            c2 = sqrt(gama*p2/rho2)
           
            call wave_type()

            if (Lflag == 1) then
                a_left = rho1*sqrt(gama*p1/rho1)*sqrt(((gama+1)*p_star/(2*gama*p1))+((gama-1)/(2*gama)))
                rho_star_left = (rho1*a_left)/(a_left-rho1*(u1-u_star))
                z_left = u1 - a_left/rho1
            else if (Lflag == 0 .and. midflag ==0) then
                c_star_left = c1 + ((gama-1)*(u1-u_star))/2
                rho_star_left = gama*p_star/(c_star_left**2)
                z_left_head = u1 - c1
                z_left_tail = u_star - c_star_left 
                write(*,*) "rho_star_left = ", rho_star_left
                write(*,*) "左侧膨胀波, z_head = ", z_left_head, "z_tail = ", z_left_tail
            else if (Lflag == 0 .and. midflag == 1) then
                z_left_head = u1 - sqrt(gama*p1/rho1)
                z_left_tail = u1 - 2*sqrt(gama*p1/rho1)/(gama-1) 
            end if

            if (Rflag == 1) then
                a_right = rho2*c2*(sqrt(((gama+1)*p_star)/(2*gama*p2)+((gama-1)/(2*gama))))
                rho_star_right = (rho2*a_right)/(a_right+rho2*(u2-u_star)) 
                z_right = u2 + a_right/rho2
                write(*,*) "rho_star_right = ", rho_star_right
                write(*,*) "右侧激波, z = ", z_right
            else if (Rflag == 0 .and. midflag == 0) then
                c_star_right = sqrt(gama*p2/rho2) + (gama-1)*(u2-u_star)/2
                rho_star_right = gama*p_star/(c_star_right**2)
                z_right_head = u2 + sqrt(gama*p2/rho2)
                z_right_tail = u_star + c_star_right
            else if (Rflag == 0 .and. midflag == 1) then
                z_right_head = u2 + sqrt(gama*p2/rho2)
                z_right_tail = u2 - 2*sqrt(gama*p2/rho2)/(gama-1)
            end if

        end subroutine wave

        subroutine region_calc()
            implicit none
            integer::j, w
            integer::datastat
            character(len=80)::datamsg
            real::c
            c1 = sqrt(gama*p1/rho1)
            c2 = sqrt(gama*p2/rho2)

            if (Lflag == 0 .and. Rflag == 0) then
                ! 左右两侧均为膨胀波
                rarefaction_left_down = z_left_head * t
                rarefaction_left_up = z_left_tail * t
                rarefaction_right_down = z_right_tail * t
                rarefaction_right_up = z_right_head * t
                mid_wave = u_star * t
                do j = 1, 201
                    if (x(j)<=rarefaction_left_down) then
                        u(j) = u1
                        p(j) = p1
                        rho(j) = rho1
                    else if (x(j)<=rarefaction_left_up) then
                        u(j) = x(j)/t + c_star_left
                        p(j) = p1*(c_star_left/sqrt(gama*p1/rho1))**(2*gama/(gama-1))
                        rho(j) = gama*p(j)/c_star_left**2
                    else if (x(j)<=mid_wave) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_left
                    else if (x(j)<=rarefaction_right_down) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_right
                    else if (x(j)<=rarefaction_right_up) then
                        u(j) = x(j) / t - c_star_right
                        p(j) = p2*(c_star_right/sqrt(gama*p2/rho2))**(2*gama/(gama-1))
                        rho(j) = gama*p(j)/c_star_right**2
                    else
                        u(j) = u2
                        p(j) = p2
                        rho(j) = rho2
                    end if
                end do
            else if (Lflag == 0 .and. Rflag == 1) then
                ! 左侧为膨胀波，右侧为激波
                rarefaction_left_down = z_left_head * t
                rarefaction_left_up = z_left_tail * t
                mid_wave = u_star * t
                shock_right = z_right * t
                write(*,*) "左侧为膨胀波: (", rarefaction_left_down, " ,", rarefaction_left_up, ")"
                write(*,*) "右侧为激波: location = ", shock_right
                write(*,*) "接触间断: loacation = ", mid_wave
                do j = 1, 201
                    if (x(j)<=rarefaction_left_down) then
                        u(j) = u1
                        p(j) = p1
                        rho(j) = rho1
                    else if (x(j)<=rarefaction_left_up) then
                        c = (gama-1)*(u1-x(j)/t)/(gama+1)+2*c1/(gama+1)
                        u(j) = x(j)/t + c 
                        p(j) = p1*((c/c1)**(2*gama/(gama-1)))
                        rho(j) = gama*p(j)/(c**2)
                    else if (x(j)<=mid_wave) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_left
                    else if (x(j)<=shock_right) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_right
                    else
                        u(j) = u2
                        p(j) = p2
                        rho(j) = rho2
                    end if
                end do
            else if (Lflag == 1 .and. Rflag == 0) then
                ! 左侧为激波，右侧为膨胀波
                rarefaction_right_down = z_right_tail * t
                rarefaction_right_up = z_right_head * t
                mid_wave = u_star * t
                shock_left = z_left * t
                do j = 1, 201
                    if (x(j)<=shock_left) then
                        u(j) = u1
                        p(j) = p1
                        rho(j) = rho1
                    else if (x(j)<=mid_wave) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_left
                    else if (x(j)<=rarefaction_right_down) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_left
                    else if (x(j)<=rarefaction_right_up) then
                        u(j) = x(j) / t - c_star_right
                        p(j) = p2*(c_star_right/sqrt(gama*p2/rho2))**(2*gama/(gama-1))
                        rho(j) = gama*p(j)/c_star_right**2
                    else
                        u(j) = u2
                        p(j) = p2
                        rho(j) = rho2
                    end if
                end do
            else if (Lflag == 1 .and. Rflag == 1) then
                ! 左右两侧均为激波
                mid_wave = u_star * t
                shock_left = z_left * t
                shock_right = z_right * t
                do j = 1, 201
                    if (x(j)<=shock_left) then
                        u(j) = u1
                        p(j) = p1
                        rho(j) = rho1
                    else if (x(j)<=mid_wave) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_left
                    else if (x(j)<=shock_right) then
                        u(j) = u_star
                        p(j) = p_star
                        rho(j) = rho_star_left
                    else
                        u(j) = u2
                        p(j) = p2
                        rho(j) = rho2
                    end if
                end do
            end if
            open(unit=30, file='wave.dat', status='old', &
                action='write', iostat=datastat, iomsg=datamsg)
            301 format (f5.2, 5x, f10.6, 5x, f15.6, 5x, f15.6)
            do w = 1, 201
                write(30, 301) x(w), u(w), p(w), rho(w)
           end do
        end subroutine region_calc

        subroutine RH_equation() 
            implicit none
            call Newton2p_star(p_star)    
            write(*,*) "p_star = ", p_star
            u_star = (u1+u2+f_wave(p2, rho2, p_star)-f_wave(p1, rho1, p_star))/2
            write(*,*) 'u_star = ', u_star
            call wave()
            call region_calc()
        end subroutine RH_equation
    end module RH_relation