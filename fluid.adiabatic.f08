! .....This program calculates the time evolution of a self-gravitating
! .....spherical gas cloud, using the Lagrangian formulation of hydrodynamics.

program fluid
	implicit none
	real*8					:: Mint, Mcore
	real*8 					:: Eint, Ekin, Egrav, Etot
	integer 				:: timesteps
	integer, parameter		:: zones = 50  
	real*8, parameter 		:: Au = 1.496e+13
	real*8, parameter		:: a = 3.0
	real*8, parameter		:: gamma = 7./5
	real*8, parameter		:: V_0 = 2.0944e+15 		! constant density
	real*8, parameter		:: Rj = 1.0e16 				! total size
	real*8, parameter		:: eta = .2 * .1    		! CLF number
	real*8, parameter		:: Rcore = 3.0e+11
	real*8, parameter		:: Msun = 2.0e+33
	real*8, parameter		:: pi =3.14159
	real*8, parameter		:: G = 6.67e-8
	real*8, parameter		:: m_proton = 1.6598e-24
	real*8, parameter		:: wm = 2.					! ?
	real*8, parameter		:: k_boltzmann = 1.38046e-16
	real*8, parameter		:: comp = k_boltzmann / (m_proton * wm)
	real*8, parameter		:: A_rad = 7.56e-15
	real*8, parameter		:: Ch = (5./2.) * 8.317e+7	! ?
	integer, parameter		:: nprint = 5
	real*8, parameter 		:: E_SN = 1e+10 			! set to 1e51 ergs later

!	bounds important to allow 0 indexing
	real*8, dimension(0:zones) 		:: Unew, Uold
	real*8, dimension(0:zones+1)	:: Rnew, Rold
	real*8, dimension(zones+1)		:: Vnew, Vold
	real*8, dimension(zones+1)		:: Pnew, Pold
	real*8, dimension(zones+1)		:: Enew, Eold
	real*8, dimension(zones+1)		:: Qnew, Qold
	real*8, dimension(0:zones+1)	:: Tnew, Told
	real*8, dimension(0:zones+1)	:: C_ad				! adiabatic sound speed
	real*8, dimension(0:zones+1)	:: Mass

	integer 	:: i, k

	integer 	:: info_unit, velos_unit, modfile_unit
	real*8 		:: totalT

	real*8 		:: c, delP, delta, deltat, term1, term2
	real*8 		:: Ta, Tb, Eb, Vb

	logical 	:: interactive 		= .False.
	logical 	:: restart_switch
	logical 	:: gravity_switch

	character(len=*), parameter	::	info_fmt = '(i6,7(1p1e14.5))'

	character(len=*), parameter	:: 	velos_fmt = '(2(i6),7(1p1e12.4))'

	open (newunit=info_unit,file='info',status='unknown')
	open (newunit=velos_unit,file='velos',status='unknown')
	open (newunit=modfile_unit,file='modfile',status='unknown')

	if (interactive) then
		!Currently not working with iPython binding
		write(*,*) 'Enter the number of timesteps:'
			read(*,*) timesteps 

		write(*,*) 'Restart from model file (True) or new simulation (False)?'
			read(*,*) restart_switch

		write(*,*) 'With gravity? Yes (True) or no (False)?'
			read(*,*) gravity_switch
	else
		timesteps 		= 50
		restart_switch 	= .False.
		gravity_switch 	= .False.
	endif


	if (restart_switch) then
		do  i=0,zones,2

			if (i.eq.0) then

				Uold(i) = 0.0
				Rold(i) = 3.0e+11

			else if (i.eq.2) then

				Uold(i) = 0.0
				Rold(i) = (real(i)/real(zones))**2.0*Rj
				Vold(i-1) = V_0
				Told(i-1)= 1.0e+6
				Pold(i-1) = comp*Told(i-1)*(1/Vold(i-1)) &
					+(.333333*A_rad*Told(i-1)**4)
				C_ad(i-1)=sqrt(gamma*(Pold(i-1)/(1/Vold(i-1))))
				Eold(i-1)=(cH/wm)*Told(i-1) &
					+A_rad*Vold(i-1)*Told(i-1)**4
				Qold(i-1) = 0.9
			else

				Uold(i) = 0.0
				Rold(i) = (real(i)/real(zones))**2.0*Rj
				Vold(i-1) = V_0
				Told(i-1)= 50.
				Pold(i-1) = comp*Told(i-1)*(1/Vold(i-1)) &
					+(.333333*A_rad*Told(i-1)**4)
				C_ad(i-1)=sqrt(gamma*(Pold(i-1)/(1/Vold(i-1))))
				Eold(i-1)=(cH/wm)*Told(i-1)  &
					+A_rad*Vold(i-1)*Told(i-1)**4
				Qold(i-1) = 0.9

			end if
		end do

!		No pressure outside boundary.
!		Other quantities reflective across boundary.

		Pold(zones+1)=0.0
		Told(zones+1)=20.0
		Eold(zones+1)=(cH/wm)*Told(zones+1)+A_rad*Vold(zones+1)*Told(zones+1)**4
		Vold(zones+1)=Vold(zones-1)
		Qold(zones+1)=Qold(zones-1)

!		Mass of the core initially starts out at zero
		Mcore = 0.0
	else

!	    Read in already existing data:
		do i=0,zones
			read(modfile_unit,*) Uold(i)
			read(modfile_unit,*) Rold(i)
			read(modfile_unit,*) Pold(i+1)
			read(modfile_unit,*) Vold(i+1)
			read(modfile_unit,*) Eold(i+1)
			read(modfile_unit,*) Told(i+1)
			read(modfile_unit,*) C_ad(i+1)
			read(modfile_unit,*) Qold(i+1)
		end do
		read(modfile_unit,*) totalT

	endif

!	Calculate the initial starting energy:
!	Mint stands for the mass internal to the current radius.

	Mint       = 0.0
	Egrav      = 0.0
	Eint       = 0.0
	Ekin       = 0.0
	Mass(0)    = 0.0
	TotalT     = 0.0
	
	do i=2,zones-2,2
		Mass(i-1) = (4./3.)*pi*(Rold(i)**3-Rold(i-2)**3)/Vold(i-1)
		Mint = Mint + Mass(i-1)
		Egrav = Egrav-(G*(4./3.)*pi*Mint*(Rold(i+2)**3-Rold(i)**3)) &
			/(((Rold(i)+Rold(i+2))/2)*Vold(i+1))
		if (gravity_switch) then
			Egrav = 0
		end if
	end do


!	Add self energy of the central zone
      Egrav = Egrav - (3./5.)*G*Mass(1)**2/Rold(2)
      if (gravity_switch) then
        Egrav = 0
      end if

!	Add the mass of the outermost zone to the total mass.
	Mass(zones-1) = 1.333*pi*(Rold(zones)**3-Rold(zones-2)**3)/Vold(zones-1)
	Mint = Mint + Mass(zones-1)

	do i=2,zones,2
        Eint=Eint+Eold(i-1)*Mass(i-1)
        Ekin=Ekin+0.5*((Uold(i)**2+Uold(i-2)**2)/2)*Mass(i-1)
	end do

	Etot = Egrav+Ekin+Eint


	write(info_unit,FMT=info_fmt ) 0 ,Etot,Egrav,Ekin,Eint,Mint,0.0,totalT


	do i=0,zones-2,2
        write(velos_unit,fmt=velos_fmt) 0,i,Rold(i),Uold(i),log10(1/Vold(i+1)), &
			Told(i+1), Mass(i+1), c_ad(i+1), Eold(i+1)
	end do


!	MAIN LOOP
	do k=1,timesteps    
		c = C_ad(1)
		deltat = eta*(Rold(2)-Rold(0))/(sqrt(c**2+(Uold(2)**2+Uold(0) &
			**2)/2)) 
        
		do i=0,zones-2,2

			c = C_ad(i+1)
			delta=eta*(Rold(i+2)-Rold(i))/(sqrt(c**2+(Uold(i+2)**2 &
				+Uold(i)**2)/2))

			if (delta.lt.deltat)  then
				deltat = delta
			endif

        end do

        Mint = 0

   		do i = 0,zones,2
          
			if (i.eq.0) then

            	Unew(0) = 0 

			elseif (i.ne.zones) then

				delP = (Pold(i+1)-Pold(i-1)+Qold(i+1)-Qold(i-1)) &
					 /((Rold(i+2)-Rold(i-2))/2)

				term1 = -.5*(Vold(i+1)+Vold(i-1))*delP

!				compute running total of interior mass.
!				then compute the gravitational force term.

				Mint = Mint + Mass(i-1)
				term2 = (-G*Mint)/(Rold(i)**2)
				if (gravity_switch) then
					term2 = 0
				end if
				Unew(i) = (term1 + term2)*deltat + Uold(i)

!			Decide here whether to fix the outer boundary condition.

			elseif (i.eq.zones) then  

        	    Unew(zones) = 0

			endif

		end do

		do i = 0,zones,2

			Rnew(i) = Unew(i)*deltat + Rold(i)

		end do
      
		do i = 0,zones-2,2

			Vnew(i+1) = 1.333*pi*(Rnew(i+2)**3-Rnew(i)**3)/Mass(i+1)
			if ((Unew(i+2)-Unew(i)).ge.0) then
				qnew(i+1)=0
			else
				qnew(i+1)=((2*a**2)*(Unew(i+2)-Unew(i))**2)/ &
					(Vnew(i+1)+Vold(i+1))
			endif
			Enew(i+1) = -Pold(i+1)*(Vnew(i+1)-Vold(i+1))+Eold(i+1)
			Ta = Told(i+1)
			Tb = Tnew(i+1)
			Vb = Vnew(i+1)
			Eb = Enew(i+1)
			call rootfinder(Eb,Vb,Ta,Tb,ch,a_rad,wm,zones)
			Tnew(i+1) = Tb
			Pnew(i+1) = comp*Tnew(i+1)*(1/Vnew(i+1)) &
				+(.33333*A_rad*Tnew(i+1)**4)
			C_ad(i+1) = sqrt(gamma*Pnew(i+1)/(1/Vnew(i+1)))

		end do





!		Get the half zone radii
		do i = 0,zones-2,2
			Rnew(i+1) = (Rnew(i) + Rnew(i+2))/2
			if (i.eq.zones-2) then
				Rnew(zones+1) = Rnew(zones) + (Rnew(zones) - Rnew(zones-1))
			endif
		end do


!		Get the new boundary conditions:
		Pnew(zones+1) = 0
		Vnew(zones+1) = Vnew(zones-1)
		Enew(zones+1) = Enew(zones-1)

		if (k.ge.0) then

!			Print out quantities for the current timestep.
       
			if ((mod(k,nprint).eq.0).and.(mod(k,1).eq.0)) then

				do i=0,zones-2,2
					write(velos_unit,FMT=velos_fmt) k,i,Rnew(i),Unew(i),log10(1/Vnew(i+1)), &
						Told(i+1), Mass(i+1), c_ad(i+1), Eold(i+1)
				end do

			end if

		end if
 
!     Update the old quanties in preparation for new timestep.
       
		do i=0,zones,2
			Uold(i) = Unew(i)
			Rold(i) = Rnew(i)
			if(i.ne.0) then
				Vold(i-1) = Vnew(i-1)
				Pold(i-1) = Pnew(i-1)
				Told(i-1) = Tnew(i-1)
				Eold(i-1) = Enew(i-1)
				Qold(i-1) = Qnew(i-1)
			endif
		end do

		Qold(zones+1) = Qnew(zones+1)
		Eold(zones+1) = Enew(zones+1)
		Pold(zones+1) = 0
		Vold(zones+1) = Vnew(zones-1)

!		Look at the energy and mass conservation for current timestep

		Mint = 0
		Egrav = 0
		Eint = 0
		Ekin = 0

		do i=2,zones-2,2
			Mint=Mint+Mass(i-1)
			Egrav=Egrav-(G*Mint*Mass(i+1))/((Rold(i)+Rold(i+2))/2)
			if (gravity_switch) then
				Egrav = 0
			end if
		end do

!		Add in the gravitational energy of the central zone.

		Egrav = Egrav - 3./5.*G*Mass(1)**2/Rold(2)
		if (gravity_switch) then
			Egrav = 0
		end if
		Mint=Mint+Mass(zones-1)
		do i=2,zones,2
			Eint=Eint+Eold(i-1)*Mass(i-1)
			Ekin=Ekin+0.5*((Uold(i)**2+Uold(i-2)**2)/2)*Mass(i-1)
		end do

		totalT = totalT+deltat

		Etot = Egrav+Ekin+Eint

		if (mod(k,1).eq.0) then
			write(info_unit,fmt=info_fmt) k,Etot,Egrav,Ekin, &
				Eint,Mint,deltat,totalt
		end if

	end do
!	END OF MAIN LOOP

!	Store model parameters for a restart:
	do i=0,zones
		write(modfile_unit,*) Uold(i)
		write(modfile_unit,*) Rold(i)
		write(modfile_unit,*) Pold(i+1)
		write(modfile_unit,*) Vold(i+1)
		write(modfile_unit,*) Eold(i+1)
		write(modfile_unit,*) Told(i+1)
		write(modfile_unit,*) C_ad(i+1)
		write(modfile_unit,*) Qold(i+1)
	end do
	write(modfile_unit,*) totalT

end


subroutine Rootfinder(Enew,Vnew,Told,Tnew,ch,a_rad,wm,zones)
	! Use the Newton-Raphson method to find the Temperature
	real*8,  intent(in)		:: Enew, Vnew, Told, ch, a_rad, wm
	integer, intent(in) 	:: zones
	real*8,  intent(out)	:: Tnew
	real*8 					:: T, E, dEdT, err
	real*8, parameter 		:: tol = 1e-6
	integer					:: i
	integer, parameter 		:: i_max =100
! 	save		! save for performance?

	T = Told
	i = 0

	E =  (ch/wm)*T+a_rad*(T**4)*Vnew              
	dEdT = (ch/wm)+4*a_rad*(T**3)*Vnew
	err = abs((E - Enew) / (Enew))

	do while ((err.gt.tol).and.(i.lt.imax))
		T = T - ((E-Enew) / dEdT )

		dEdT = (ch/wm)+4*a_rad*(T**3)*Vnew
		E =  (ch/wm)*T+a_rad*(T**4)*Vnew
            
		err = abs((E - Enew) / (Enew))
		i = i + 1
	end do
!        write(*,*) i
!        write(*,*) err
!        write(*,*) Tnew
	Tnew = T

	return
end