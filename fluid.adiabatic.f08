! .....This program calculates the time evolution of a self-gravitating
! .....spherical gas cloud, using the Lagrangian formulation of hydrodynamics.

program fluid
    implicit none
    real(kind=8)               :: M_core
    real(kind=8)               :: E_int, E_kin, E_grav, E_tot
    integer,      parameter    :: zones       = 500  
    real(kind=8), parameter    :: Au          = 1.496d+13
    real(kind=8), parameter    :: a           = 3.0
    real(kind=8), parameter    :: gamma       = 7./5
    real(kind=8), parameter    :: V_0         = 2.0944d+15         ! constant initial density
    real(kind=8), parameter    :: T_0         = 50                 ! constant initial temperature
    real(kind=8), parameter    :: R_tot       = 1d6 * Au           ! total size
    real(kind=8), parameter    :: eta         = .2 * .5                 ! CLF number
    real(kind=8), parameter    :: R_core      = 3.0d+11
    real(kind=8), parameter    :: M_sun       = 2.0d+33
    real(kind=8), parameter    :: pi          = 3.14159
    real(kind=8), parameter    :: G           = 6.67d-8
    real(kind=8), parameter    :: m_proton    = 1.6598d-24
    real(kind=8), parameter    :: wm          = 2.                    ! ?
    real(kind=8), parameter    :: k_boltzmann = 1.38046d-16
    real(kind=8), parameter    :: comp        = k_boltzmann / (m_proton * wm)
    real(kind=8), parameter    :: A_rad       = 7.56d-15
    real(kind=8), parameter    :: Ch          = (5./2.) * 8.317d+7    ! ?
    integer,      parameter    :: nprint      = 1000
    real(kind=8), parameter    :: E_SN        = 1.d+50                 ! set to 1e51 ergs later

!    bounds important to allow 0 indexing
    real(kind=8), dimension(0:zones)        :: U_old, U_new
    real(kind=8), dimension(0:zones+1)      :: R_old, R_new
    real(kind=8), dimension(zones+1)        :: V_old, V_new
    real(kind=8), dimension(zones+1)        :: P_old, P_new
    real(kind=8), dimension(zones+1)        :: E_old, E_new
    real(kind=8), dimension(zones+1)        :: Q_old, Q_new
    real(kind=8), dimension(0:zones+1)      :: T_old, T_new
    real(kind=8), dimension(0:zones+1)      :: C_ad                ! adiabatic sound speed
    real(kind=8), dimension(0:zones+1)      :: Mass
    real(kind=8), dimension(0:zones+1)      :: M_int

    integer      :: i, k

    integer      :: info_unit, velos_unit, modfile_unit
    character(len=*), parameter    ::    info_fmt  = '(i6,7(1pE15.5E3))'
    character(len=*), parameter    ::    velos_fmt = '(2(i6), 8(1pE14.5E3))'

!   OLD formatting -- doesn't properly handle double precision
!     character(len=*), parameter    ::    info_fmt  = '(i6,7(1p1e14.5))'
!     character(len=*), parameter    ::    velos_fmt = '(2(i6),8(1p1d13.4))'

    real(kind=8) :: time_total

!    declare sundry temporary parameters
    real(kind=8) :: c, del_P, delta, delta_time, term_1, term_2
    real(kind=8) :: T_a, T_b, E_b, V_b

    integer      :: timesteps
    logical      :: interactive   = .False.
    logical      :: from_previous
    logical      :: with_gravity





    open (newunit=info_unit,    file='info',    status='unknown')
    open (newunit=velos_unit,   file='velos',   status='unknown')
    open (newunit=modfile_unit, file='modfile', status='unknown')

    if (interactive) then
        !Currently not working with iPython binding
        write(*,*) 'Enter the number of timesteps:'
            read(*,*) timesteps 

        write(*,*) 'Restart from model file (True) or new simulation (False)?'
            read(*,*) from_previous

        write(*,*) 'Exclude gravity? Yes (True) or no (False)?'
            read(*,*) with_gravity
    else
        timesteps         = 100000
        from_previous     = .False.
        with_gravity      = .False.
    endif


    if (.NOT.from_previous) then

!         U_old(0:zones:2) = 0.0
!         R_old(0:zones:2) = [((real(i) / real(zones))**2.0 * R_tot, i=0,zones,2)]
!         V_old(1:zones-1:2) = V_0
!         T_old(1:zones-1:2) = T_0
!         P_old(1:zones-1:2) = [(comp*T_old(i)*(1/V_old(i)) &
!                     +(.333333*A_rad*T_old(i)**4), i=1,zones-1,2)]
!         c_ad(1:zones-1:2) = [(sqrt(gamma*(P_old(i)/(1/V_old(i)))), i=1, zones-1, 2)]
!         E_old(1:zones-1:2) = [((cH/wm)*T_old(i) &
!                     +A_rad*V_old(i)*T_old(i)**4, i=1, zones-1, 2)]
!         Q_old(i:zones-1:2) = 0.9

!         R_old(0) = 3.0E+11     ! define this inner radius elsewhere

        do i = 0, zones, 2

            if (i.eq.0) then

                U_old(i) = 0.0
                R_old(i) = 3.0d+11

            else if (i.eq.2) then



                U_old(i)      = 0.0
                R_old(i)      = (real(i) / real(zones))**2.0 * R_tot
                V_old(i-1)    = V_0
!                 T_old(i-1)    = 1e+6        ! inject SNe energy
                call Tsolver(E_SN, V_old(i-1), T_0, T_old(i-1), &
                    ch, a_rad, wm)
                P_old(i-1)    = comp * T_old(i-1) * (1 / V_old(i-1)) &
                    + (.333333 * A_rad * T_old(i-1)**4)
                C_ad(i-1)     = sqrt(gamma * (P_old(i-1) * V_old(i-1)))
                E_old(i-1)    = (cH / wm) * T_old(i-1) &
                    + A_rad * V_old(i-1) * T_old(i-1)**4
                Q_old(i-1)    = 0.9
            else

                U_old(i)     = 0.0
                R_old(i)     = (real(i) / real(zones))**2.0 * R_tot
                V_old(i-1)   = V_0
                T_old(i-1)   = T_0
                P_old(i-1)   = comp * T_old(i-1) * (1 / V_old(i-1)) &
                    + (.333333 * A_rad * T_old(i-1)**4)
                C_ad(i-1)    = sqrt(gamma * (P_old(i-1) * V_old(i-1)))
                E_old(i-1)   = (cH / wm) * T_old(i-1)  &
                    + A_rad * V_old(i-1) * T_old(i-1)**4
                Q_old(i-1)   = 0.9

            end if
        end do

!         No pressure outside boundary.
!         Other quantities reflective across boundary.

        P_old(zones+1) = 0.0
        T_old(zones+1) = T_0
        V_old(zones+1) = V_old(zones-1)
        E_old(zones+1) = (cH / wm) * T_old(zones+1) &
            + A_rad * V_old(zones+1) * T_old(zones+1)**4
        Q_old(zones+1) = Q_old(zones-1)

!        Mass of the core initially starts out at zero
        M_core = 0.0
    else

!        Read in already existing data:
        do i = 0, zones
            read(modfile_unit,*) U_old(i)
            read(modfile_unit,*) R_old(i)
            read(modfile_unit,*) P_old(i+1)
            read(modfile_unit,*) V_old(i+1)
            read(modfile_unit,*) E_old(i+1)
            read(modfile_unit,*) T_old(i+1)
            read(modfile_unit,*) C_ad(i+1)
            read(modfile_unit,*) Q_old(i+1)
        end do
        read(modfile_unit,*) time_total

    endif

!    Calculate the initial starting energy:
!    M_int stands for the mass internal to that zone border
!        e.g. M_int(i+1) includes all Mass(1:i:2) 

    E_grav      = 0.0
    E_int       = 0.0
    E_kin       = 0.0
    Mass(0)     = 0.0
    M_int(0)    = Mass(0)
    time_total  = 0.0

    do i = 2, zones-2, 2
        Mass(i-1)  = (4./3.) * pi * (R_old(i)**3 - R_old(i-2)**3) &
            / V_old(i-1)
        M_int(i)   = M_int(i-2) + Mass(i-1)
        E_grav     = E_grav - (G * (4./3.) * pi * M_int(i) * (R_old(i+2)**3 &
            - R_old(i)**3)) &
            /(((R_old(i) + R_old(i+2)) / 2) * V_old(i+1))
    end do

    if (.NOT.with_gravity) then
            E_grav = 0
    end if

!     Mass(1:zones-1:2) = [((4./3.)*pi*(R_old(i+1)**3-R_old(i-1)**3)/V_old(i), i=1, zones-1,2)]
!     M_int(1:zones+1:2) = [(M_int(i-1) + Mass(i), i=1, zones+1, 2)]
!     if (with_gravity) then
!         E_grav     = sum ([(E_grav-(G*(4./3.)*pi*M_int(i)*(R_old(i+2)**3-R_old(i)**3)) &
!                 /(((R_old(i)+R_old(i+2))/2)*V_old(i+1)), i=2, zones, 2)])
!     endif


!    Add self energy of the central zone
    E_grav = E_grav - (3./5.)*G*Mass(1)**2/R_old(2)
    if (.NOT.with_gravity) then
        E_grav = 0
    end if

!    Add the mass of the outermost zone to the total mass.
    Mass(zones-1)     = 1.333 * pi * (R_old(zones)**3 - R_old(zones-2)**3) &
        / V_old(zones-1)
    M_int(zones)      = M_int(zones-2) + Mass(zones-1)

    do i = 2, zones, 2
        E_int    = E_int + E_old(i-1) * Mass(i-1)
        E_kin    = E_kin + 0.5 * ((U_old(i)**2 + U_old(i-2)**2) / 2) * Mass(i-1)
    end do

    E_tot = E_grav + E_kin + E_int


    write(info_unit, FMT=info_fmt)  0, E_tot, E_grav, E_kin, E_int, &
        M_int(zones), 0.0, time_total


    do i = 0, zones-2, 2
        write(velos_unit, fmt=velos_fmt)  0, i, R_old(i), U_old(i), &
        log10(1/V_old(i+1)), T_old(i+1), Mass(i+1), c_ad(i+1), E_old(i+1), P_old(i+1)
    end do


!    MAIN LOOP
    do k = 1, timesteps    
        c = C_ad(1)
        delta_time = eta * (R_old(2) - R_old(0)) &
        / (sqrt(c**2 + (U_old(2)**2 + U_old(0)**2) / 2)) 
        
!        calculate delta_time to match desired CFL number, eta
!        WARNING - calculates delta_time *before* updating values
!                - if there's a strong shock, the initial + final CFLs
!                    of a cell might differ significantly
        do i = 0, zones-2, 2
            c       = C_ad(i+1)
            delta = eta * (R_old(i+2) - R_old(i)) &
                / (sqrt(c**2 + (U_old(i+2)**2 + U_old(i)**2) / 2))

            if (delta.lt.delta_time)  then
                delta_time = delta
            endif

        end do

           do i = 0,zones,2

            if (i.eq.0) then
                U_new(0) = 0 

            elseif (i.ne.zones) then
                del_P = (P_old(i+1) - P_old(i-1) &
                    +    Q_old(i+1) - Q_old(i-1)) &
                    /  ((R_old(i+2) - R_old(i-2)) / 2)

!                Pressure + Viscous forces
                term_1 = -.5 * (V_old(i+1) + V_old(i-1)) * del_P

!                Use the running total of interior mass.
!                to compute the gravitational force term.

                term_2 = (-G*M_int(i))/(R_old(i)**2)
                if (.NOT.with_gravity) then
                    term_2 = 0
                end if

                U_new(i) = (term_1 + term_2) * delta_time + U_old(i)

!            Decide here whether to fix the outer boundary condition.
            elseif (i.eq.zones) then  
                U_new(zones) = 0
            endif

        end do

        do i = 0, zones,2

            R_new(i) = U_new(i) * delta_time + R_old(i)

        end do
      
        do i = 0, zones-2,2

            V_new(i+1) = 1.333 * pi * (R_new(i+2)**3 - R_new(i)**3) / Mass(i+1)
            if ((U_new(i+2) - U_new(i)).ge.0) then
                Q_new(i+1) = 0
            else
                Q_new(i+1) = ((2*a**2) * (U_new(i+2) - U_new(i))**2) / &
                    (V_new(i+1) + V_old(i+1))
            endif
            E_new(i+1) = -P_old(i+1) * (V_new(i+1) - V_old(i+1)) + E_old(i+1)
            T_a = T_old(i+1)
            T_b = T_new(i+1)
            V_b = V_new(i+1)
            E_b = E_new(i+1)
!             call Tsolver(E_b, V_b, T_a, T_b, ch, a_rad, wm)
!             T_new(i+1) = T_b
!             T_a = T_a
            call Tsolver(E_new(i+1), V_new(i+1), T_old(i+1), T_new(i+1), &
                ch, a_rad, wm)
            P_new(i+1) = comp * T_new(i+1) * (1/V_new(i+1)) &
                + (.33333 * A_rad * T_new(i+1)**4)
            C_ad(i+1)  = sqrt(gamma * P_new(i+1) * V_new(i+1))

        end do





!        Get the half zone radii
        do i = 0, zones-2, 2
            R_new(i+1) = (R_new(i) + R_new(i+2))/2
            if (i.eq.zones-2) then
                R_new(zones+1) = R_new(zones) + (R_new(zones) - R_new(zones-1))
            endif
        end do


!        Get the new boundary conditions:
        P_new(zones+1) = 0
        V_new(zones+1) = V_new(zones-1)
        E_new(zones+1) = E_new(zones-1)

        if (k.ge.0) then

!            Print out quantities for the current timestep.
       
            if ((mod(k,nprint).eq.0).and.(mod(k,1).eq.0)) then

                do i=0, zones-2, 2
                    write(velos_unit,FMT=velos_fmt) k,i,R_new(i),U_new(i),log10(1/V_new(i+1)), &
                        T_old(i+1), Mass(i+1), c_ad(i+1), E_old(i+1), P_old(i+1)
                end do

            end if

        end if
 
!     Update the old quanties in preparation for new timestep.
       
        do i=0, zones, 2
            U_old(i) = U_new(i)
            R_old(i) = R_new(i)
            if(i.ne.0) then
                V_old(i-1) = V_new(i-1)
                P_old(i-1) = P_new(i-1)
                T_old(i-1) = T_new(i-1)
                E_old(i-1) = E_new(i-1)
                Q_old(i-1) = Q_new(i-1)
            endif
        end do

        Q_old(zones+1) = Q_new(zones+1)
        E_old(zones+1) = E_new(zones+1)
        P_old(zones+1) = 0
        V_old(zones+1) = V_new(zones-1)

!        Look at the energy and mass conservation for current timestep

        E_grav = 0
        E_int  = 0
        E_kin  = 0

        do i = 2, zones-2, 2
            E_grav=E_grav-(G*M_int(i)*Mass(i+1))/((R_old(i)+R_old(i+2))/2)
            if (.NOT.with_gravity) then
                E_grav = 0
            end if
        end do

!        Add in the gravitational energy of the central zone.

        E_grav = E_grav - 3./5.*G*Mass(1)**2/R_old(2)
        if (.NOT.with_gravity) then
            E_grav = 0
        end if
!         M_int=M_int+Mass(zones-1)
        do i = 2, zones, 2
            E_int = E_int + E_old(i-1)*Mass(i-1)
            E_kin = E_kin + 0.5*((U_old(i)**2 + U_old(i-2)**2)/2) * Mass(i-1)
        end do

        time_total = time_total+delta_time

        E_tot = E_grav+E_kin+E_int

        if (mod(k,1).eq.0) then
            write(info_unit,fmt=info_fmt) k,E_tot,E_grav,E_kin, &
                E_int,M_int(zones),delta_time,time_total
        end if

    end do
!    END OF MAIN LOOP

!    Store model parameters for a restart:
    do i = 0, zones
        write(modfile_unit,*) U_old(i)
        write(modfile_unit,*) R_old(i)
        write(modfile_unit,*) P_old(i+1)
        write(modfile_unit,*) V_old(i+1)
        write(modfile_unit,*) E_old(i+1)
        write(modfile_unit,*) T_old(i+1)
        write(modfile_unit,*) C_ad(i+1)
        write(modfile_unit,*) Q_old(i+1)
    end do
    write(modfile_unit,*) time_total

end


subroutine Tsolver(E_new,V_new,T_old,T_new,ch,a_rad,wm)
    ! Use the Newton-Raphson method to find the Temperature
    real(kind=8), intent(in)    :: E_new, V_new, T_old, ch, a_rad, wm
    real(kind=8), intent(out)   :: T_new
    real(kind=8)                :: E, dEdT, err
    real(kind=8), parameter     :: tol = 1d-6
    integer                     :: feval=0
    integer, parameter          :: feval_max = 100
!   save        ! save for performance?

    T_new = T_old

    dEdT   = (ch/wm)         + 4 * a_rad * (T_new**3) * V_new
    E      = (ch/wm) * T_new +     a_rad * (T_new**4) * V_new
            
        
    err = abs((E - E_new) / (E_new))

    do while ((err.gt.tol).and.(feval.lt.feval_max))
        T_new  = T_new - ((E - E_new) / dEdT )

        dEdT   = (ch/wm)         + 4 * a_rad * (T_new**3) * V_new
        E      = (ch/wm) * T_new +     a_rad * (T_new**4) * V_new
            
        err   = abs((E - E_new) / (E_new))
        feval = feval + 1
    end do
!        write(*,*) i
!        write(*,*) err
!        write(*,*) T_new
    return
end