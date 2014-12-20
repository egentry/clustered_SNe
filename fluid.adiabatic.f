c.....This program calculates the time evolution of a self-gravitating
c.....spherical gas cloud, using the Lagrangian formulation of hydrodynamics.

      program fluid

      implicit real*8 (a-h,o-z), integer*4 (i-n)
      real*8 Msun,Mint,Mass,Mcore
      integer Tsteps 

c.... J = Number of Lagrangian zones.
c.... gamma = Adiabatic index
c.... Vo = Specific volume.
c.... Rj = Radius of system.
c.... Cinit = Initial sound speed.
c.... Tsteps = Number of timesteps that the program is to run for.
c.... Msun = Initial mass of the protostellar system.


c.....Set J (number of zones) both here and in subroutine RadT
      parameter (Au =1.496e+13,a=3.0,J=200,gamma=7./5.,Vo=2.0944e+15)
      parameter (Rj=1.0e+16,eta=.20 * .1)
      parameter (Rcore=3.0e+11)
      parameter (Msun=2.0e+33,pi=3.14159,G=6.67e-8)
      parameter (proton=1.6598e-24,wm=2.,boltz=1.38046e-16)
      parameter (comp=boltz/(proton*wm))
      parameter (Arad=7.56e-15,Ch=(5./2.)*(8.317e+7))
      parameter (nprint=5)
      parameter (ESN=1e+10)


      dimension Unew(0:J),Rnew(0:J+1),Vnew(J+1),Pnew(J+1),Enew(J+1) 
      dimension Uold(0:J),Rold(0:J+1),Vold(J+1),Pold(J+1),Eold(J+1)
      dimension Qnew(J+1),Qold(J+1)
      dimension Told(0:J+1),Tnew(0:J+1),Cad(0:J+1)
      dimension Mass(0:J+1) 
 

      open (unit=1,file='info',status='unknown')
      open (unit=2,file='velos',status='unknown')
      open (unit=3,file='modfile',status='unknown')


      write(*,*) 'Enter the number of timesteps:'
      read(*,*) Tsteps 

      write(*,*) 'Restart from model file (1) or new simulation (0)?'
      read(*,*) iswitch

c..      write(*,*) 'With gravity? Yes (1) or no (0)?'
c..      read(*,*) gswitch
      gswitch=0

      if (iswitch.eq.0) then
  
c..     Set up parameters to make a new simulation..
c..     Current setup: Uniform density, temperature startup condition.

        do i=0,J,2

          if (i.eq.0) then

            Uold(i) = 0.0
            Rold(i) = 3.0e+11

          else if (i.eq.2) then

            Uold(i) = 0.0
            Rold(i) = (real(i)/real(j))**2.0*Rj
            Vold(i-1) = Vo
            Told(i-1)= 1.0e+6
            Pold(i-1) = comp*Told(i-1)*(1/Vold(i-1))
     +                  +(.333333*Arad*Told(i-1)**4)
            Cad(i-1)=sqrt(gamma*(Pold(i-1)/(1/Vold(i-1))))
            Eold(i-1)=(cH/wm)*Told(i-1) 
     +                  +arad*Vold(i-1)*Told(i-1)**4
            Qold(i-1) = 0.9


          else

            Uold(i) = 0.0
            Rold(i) = (real(i)/real(j))**2.0*Rj
            Vold(i-1) = Vo
            Told(i-1)= 50.
            Pold(i-1) = comp*Told(i-1)*(1/Vold(i-1))
     +                  +(.333333*Arad*Told(i-1)**4)
            Cad(i-1)=sqrt(gamma*(Pold(i-1)/(1/Vold(i-1))))
            Eold(i-1)=(cH/wm)*Told(i-1) 
     +                  +arad*Vold(i-1)*Told(i-1)**4
            Qold(i-1) = 0.9

          end if

        end do

c..     No pressure outside boundary.
c..     Other quantities reflective across boundary.

        Pold(J+1)=0.0
        Told(J+1)=20.0
        Eold(J+1)=(cH/wm)*Told(J+1)+arad*Vold(J+1)*Told(J+1)**4
        Vold(J+1)=Vold(J-1)
        Qold(J+1)=Qold(J-1)

c..     Mass of the core initially starts out at zero
        Mcore = 0.0

      elseif (iswitch.eq.1) then

c..     Read in already existing data:

        do i=0,J
          read(3,*) Uold(i)
          read(3,*) Rold(i)
          read(3,*) Pold(i+1)
          read(3,*) Vold(i+1)
          read(3,*) Eold(i+1)
          read(3,*) Told(i+1)
          read(3,*) Cad(i+1)
          read(3,*) Qold(i+1)
        end do
        read(3,*) totalT

      end if

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

c..     Incorrect energy injection:
c..      Uold(2) = sqrt(2 * ESN / Mass(1))

c..     Better energy injection (but still doesn't give good T values)
c..      Eold(2) = ESN

c..   Add self energy of the central zone
      Egrav = Egrav - (3./5.)*G*Mass(1)**2/Rold(2)
      if (gswitch.eq.0) then
        Egrav = 0
      end if

c..   Add the mass of the outermost zone to the total mass.
      Mass(j-1) = 1.333*pi*(Rold(j)**3-Rold(j-2)**3)/Vold(j-1)
      Mint = Mint + Mass(j-1)

      do i=2,J,2
        Eint=Eint+Eold(i-1)*Mass(i-1)
        Ekin=Ekin+0.5*((Uold(i)**2+Uold(i-2)**2)/2)*Mass(i-1)
      end do

      Etot = Egrav+Ekin+Eint

      write(1,108) 0.0,Etot,Egrav,Ekin,Eint,Mint,0.0,totalT

      do i=0,J-2,2
        write(2,107) 0,i,Rold(i),Uold(i),log10(1/Vold(i+1)),
     +            Told(i+1), Mass(i+1), cad(i+1), Eold(i+1)
      end do

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

        do i = 0,J,2
          
          if (i.eq.0) then

            Unew(0) = 0 

          elseif (i.ne.J) then

            delP = (Pold(i+1)-Pold(i-1)+Qold(i+1)-Qold(i-1))
     +             /((Rold(i+2)-Rold(i-2))/2)

            term1 = -.5*(Vold(i+1)+Vold(i-1))*delP

c..         compute running total of interior mass.
c..         then compute the gravitational force term.

            Mint = Mint + Mass(i-1)
            term2 = (-G*Mint)/(Rold(i)**2)
            if (gswitch.eq.0) then
              term2 = 0
            end if
            Unew(i) = (term1 + term2)*deltat + Uold(i)

c..         Decide here whether to fix the outer boundary condition.

          elseif (i.eq.J) then  

            Unew(J) = 0

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
          Enew(i+1) = -Pold(i+1)*(Vnew(i+1)-Vold(i+1))+Eold(i+1)
          Ta = Told(i+1)
          Tb = Tnew(i+1)
          Vb = Vnew(i+1)
          Eb = Enew(i+1)
           call rootfinder(Eb,Vb,Ta,Tb,ch,arad,wm)
           Tnew(i+1) = Tb
          Pnew(i+1) = comp*Tnew(i+1)*(1/Vnew(i+1))
     +                 +(.33333*Arad*Tnew(i+1)**4)
          Cad(i+1) = sqrt(gamma*Pnew(i+1)/(1/Vnew(i+1)))

             
       end do





c..     Get the half zone radii
        do i = 0,J-2,2
          Rnew(i+1) = (Rnew(i) + Rnew(i+2))/2
          if (i.eq.J-2) then
            Rnew(J+1) = Rnew(J) + (Rnew(J) - Rnew(J-1))
          endif
        end do


c..     Get the new boundary conditions:
        Pnew(J+1) = 0
        Vnew(J+1) = Vnew(J-1)
        Enew(J+1) = Enew(J-1)

        if (k.ge.0) then

c..       Print out quantities for the current timestep.
       
          if ((mod(k,nprint).eq.0).and.(mod(k,1).eq.0)) then

            do i=0,J-2,2
              write(2,107) k,i,Rnew(i),Unew(i),log10(1/Vnew(i+1)),
     +                     Told(i+1), Mass(i+1), cad(i+1), Eold(i+1)
            end do

          end if

        end if

107     format (2(i6),7(1p1e12.4))
 
c..     Update the old quanties in preparation for new timestep.
       
        do i=0,J,2
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

        Qold(J+1) = Qnew(J+1)
        Eold(J+1) = Enew(J+1)
        Pold(J+1) = 0
        Vold(J+1) = Vnew(J-1)

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
          Eint=Eint+Eold(i-1)*Mass(i-1)
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

c..   Store model parameters for a restart:
      do i=0,J
        write(3,*) Uold(i)
        write(3,*) Rold(i)
        write(3,*) Pold(i+1)
        write(3,*) Vold(i+1)
        write(3,*) Eold(i+1)
        write(3,*) Told(i+1)
        write(3,*) Cad(i+1)
        write(3,*) Qold(i+1)
      end do
      write(3,*) totalT

        end

       subroutine Rootfinder(Enew,Vnew,Told,Tnew,ch,arad,wm)

         implicit real*8 (a-h,o-z), integer*4 (i-n)
         parameter (tol=1e-6)
         parameter (imax = 100)
c..         dimension Enew(J+1),Vnew(J+1),Told(0:J+1)
         save
        T = Told
        i = 0

        E =  (ch/wm)*T+arad*(T**4)*Vnew              
        dEdT = (ch/wm)+4*arad*(T**3)*Vnew
        err = abs((E - Enew) / (Enew))


        do while ((err.gt.tol).and.(i.lt.imax))
            T = T - ((E-Enew) / dEdT )

            dEdT = (ch/wm)+4*arad*(T**3)*Vnew
            E =  (ch/wm)*T+arad*(T**4)*Vnew
            
            err = abs((E - Enew) / (Enew))
            i = i + 1
        end do
c..        write(*,*) i
c..        write(*,*) err
c..        write(*,*) Tnew
        Tnew = T

        return
        end
