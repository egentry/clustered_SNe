c.....This program calculates the time evolution of a self-gravitating
c.....spherical gas cloud, using the Lagrangian formulation of hydrodynamics.

      program fluid

      implicit real*8 (a-h,o-z), integer*4 (i-n)
      real*8 Msun,Mint,Mass
      integer Tsteps 

c.... J = Number of Lagrangian zones.
c.... gamma = Adiabatic index
c.... Vo = Specific volume.
c.... Rj = Radius of system.
c.... Cinit = Initial sound speed.
c.... Tsteps = Number of timesteps that the program is to run for.
c.... Msun = Initial mass of the protostellar system.


c.....Set J (number of zones) both here and in subroutine RadT
      parameter (Au =1.496e+13,a=3.0,J=200,gamma=1.0,Vo=2.096e+18)
      parameter (Rj=1.0e+17,eta=0.4,Tinit=10.)
      parameter (Rcore=0.0)
      parameter (Msun=2.0e+33,pi=3.14159,G=6.67e-8)
      parameter (proton=1.6598e-24,wm=2.,boltz=1.38046e-16)
      parameter (comp=boltz/(proton*wm))
      parameter (Arad=7.56e-15,Ch=(5./2.)*(8.317e+7))
      parameter (tsteps=50,nprint=1)
      parameter (ESN=1e+30)


      dimension Unew(0:J),Rnew(0:J+1),Vnew(J+1),Pnew(J+1)
      dimension Uold(0:J),Rold(0:J+1),Vold(J+1),Pold(J+1)
      dimension Qnew(J+1),Qold(J+1)
      dimension Told(0:J+1),Cad(0:J+1)
      dimension Mass(0:J+1) 
      

      open (unit=1,file='info',status='unknown')
      open (unit=2,file='velos',status='unknown')

      write(*,*) 'With gravity? Yes (1) or no (0)?'
      read(*,*) gswitch
  
c..   Set up parameters to make a new simulation..

      do i=0,J,2

        if (i.eq.0) then

          Uold(i) = 0.0
          Rold(i) = Rcore

        else

          Uold(i) = 0.0
          Rold(i) = Rcore+(real(i)/real(j))**2.0*(Rj-Rcore)
          Rold(i-1) = Rcore+(real(i-1)/real(j))**2.0*(Rj-Rcore)
          Vold(i-1) = Vo
          Told(i-1)= Tinit
          Pold(i-1) = comp*Told(i-1)*(1/Vold(i-1))
          Cad(i-1)=sqrt(gamma*(Pold(i-1)/(1/Vold(i-1))))
          Qold(i-1) = 0.9

        end if

      end do

c..   Intensive quantities reflective across boundary.

      Pold(J+1)=Pold(J-1)
      Told(J+1)=Told(J-1)
      Vold(J+1)=Vold(J-1)
      Qold(J+1)=Qold(J-1)

c..   Calculate the initial starting energy:
c..   Mint stands for the mass internal to the current radius.

      Mint       = 0.0
      Egrav      = 0.0
      Eint       = 0.0
      Ekin       = 0.0
      Mass(0)    = 0.0
      TotalT     = 0.0
 
      do i=2,J-2,2
        Mass(i-1) = (4./3.)*pi*(Rold(i)**3-Rold(i-2)**3)/Vold(i-1)
        Mint = Mint + Mass(i-1)
        Egrav = Egrav-(G*(4./3.)*pi*Mint*(Rold(i+2)**3-Rold(i)**3))
     +           /(((Rold(i)+Rold(i+2))/2)*Vold(i+1))
        if (gswitch.eq.0) then
          Egrav = 0
        end if

      end do

c..   Inject SN Energy
      Uold(2) = sqrt(2 * ESN / Mass(1))

c..   Add self energy of the central zone
      Egrav = Egrav - (3./5.)*G*Mass(1)**2/Rold(2)
      if (gswitch.eq.0) then
        Egrav = 0
      end if


c..   Add the mass of the outermost zone to the total mass.
      Mass(j-1) = 1.333*pi*(Rold(j)**3-Rold(j-2)**3)/Vold(j-1)
      Mint = Mint + Mass(j-1)
   

      do i=2,J,2
        Eint=Eint+(5./2.)*comp*Told(i-1)*Mass(i-1)
        Ekin=Ekin+0.5*((Uold(i)**2+Uold(i-2)**2)/2)*Mass(i-1)
      end do

      Etot = Egrav+Ekin+Eint
      write(*,*) Etot



c..   MAIN LOOP
      do k=1,Tsteps    
        c = Cad(1)
        deltat = eta*(Rold(2)-Rold(0))/(sqrt(c**2+(Uold(2)**2+Uold(0)
     +           **2)/2)) 
        do i=0,J-2,2

          c = Cad(i+1)
          delta=eta*(Rold(i+2)-Rold(i))/(sqrt(c**2+(Uold(i+2)**2
     +          +Uold(i)**2)/2))

          if (delta.lt.deltat)  then
            deltat = delta
          endif

        end do

        Mint = 0

c..     rigid inner wall boundary conditions

        do i = 0,J,2
          
          if ((i.eq.0).or.(i.eq.J)) then

            Unew(0) = 0. 

          else

            delP = (Pold(i+1)-Pold(i-1)+Qold(i+1)-Qold(i-1))
     +             /(Rold(i+1)-Rold(i-1))

            term1 = -.5*(Vold(i+1)+Vold(i-1))*delP

c..         compute running total of interior mass.
c..         then compute the gravitational force term.

            Mint = Mint + Mass(i-1)
            term2 = (-G*Mint)/(Rold(i)**2)
            if (gswitch.eq.0) then
              term2 =0
            end if

            
            Unew(i) = (term1 + term2)*deltat + Uold(i)


          endif

        end do

        do i = 0,J,2

          Rnew(i) = Unew(i)*deltat + Rold(i)

        end do
      
        do i = 0,J-2,2

          Vnew(i+1) = 1.333*pi*(Rnew(i+2)**3-Rnew(i)**3)/Mass(i+1)
          if ((Unew(i+2)-Unew(i)).ge.0) then
            qnew(i+1)=0
          else
            qnew(i+1)=((2*a**2)*(Unew(i+2)-Unew(i))**2)/
     +                (Vnew(i+1)+Vold(i+1))
          endif
             
       end do

       do i=0,j-2,2

         Pnew(i+1) = comp*Told(i+1)*(1/Vnew(i+1))
         Cad(i+1) = sqrt(gamma*Pnew(i+1)/(1/Vnew(i+1)))

       end do


c..    Get the half zone radii
       do i = 0,J-2,2
         Rnew(i+1) = (Rnew(i) + Rnew(i+2))/2
         if (i.eq.J-2) then
           Rnew(J+1) = Rnew(J) + (Rnew(J) - Rnew(J-1))
         endif
       end do


c..     Get the new boundary conditions:
        Pnew(J+1) = Pnew(J-1)
        Vnew(J+1) = Vnew(J-1)

        if (k.ge.0) then

c..       Print out quantities for the current timestep.
       
          if ((mod(k,nprint).eq.0).and.(mod(k,1).eq.0)) then

            do i=0,J,2
              write(2,107) k,i,Rnew(i)/1.5e+13,Unew(i)/2.0e+4,
     +                      (1./Vnew(i+1))*1.0e+18
            end do

          end if

        end if

107     format (2(i6),4(1p1e12.4))
 
c..     Update the old quanties in preparation for new timestep.
       
        do i=0,J,2
          Uold(i) = Unew(i)
          Rold(i) = Rnew(i)
          if(i.ne.0) then
            Vold(i-1) = Vnew(i-1)
            Pold(i-1) = Pnew(i-1)
            Qold(i-1) = Qnew(i-1)
          endif
        end do

        Qold(J+1) = Qnew(J+1)
        Pold(J+1) = Pnew(J+1)
        Vold(J+1) = Vnew(J+1)

c..     Look at the energy and mass conservation for current timestep

        Mint = 0
        Egrav = 0
        Eint = 0
        Ekin = 0

        do i=2,J-2,2
          Mint=Mint+Mass(i-1)
          Egrav=Egrav-(G*Mint*Mass(i+1))/((Rold(i)+Rold(i+2))/2)
          if (gswitch.eq.0) then
            Egrav = 0
          end if
        end do

c..     Add in the gravitational energy of the central zone.

        Egrav = Egrav - 3./5.*G*Mass(1)**2/Rold(2)
        if (gswitch.eq.0) then
          Egrav = 0
        end if
        Mint=Mint+Mass(j-1)
        do i=2,J,2
          Eint=Eint+(5./2.)*comp*Told(i-1)*Mass(i-1)
          Ekin=Ekin+0.5*((Uold(i)**2+Uold(i-2)**2)/2)*Mass(i-1)
        end do

        totalT = totalT+deltat

        Etot = Egrav+Ekin+Eint

        if (mod(k,1).eq.0) then
          write(1,108) k,Etot,Egrav,Ekin,
     +                Eint,Mint,deltat,totalt
        end if

108     format (i6,7(1p1e14.5))

      end do
c..   END OF MAIN LOOP

      end


