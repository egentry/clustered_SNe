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
      dimension Told(0:J+1),Cad(0:J+1)
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

      endif

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
             
       end do

c......Return the temperature array from the freshly calculated energies.
       call Rootfinder(Enew,Vnew,Told,Ch,Arad,wm)

c......Calculate the pressures and the soundspeeds
       do i=0,j-2,2

         Pnew(i+1) = comp*Told(i+1)*(1/Vnew(i+1))
     +                 + (.33333*Arad*Told(i+1)**4)
         Cad(i+1) = sqrt(gamma*Pnew(i+1)/(1/Vnew(i+1)))

       end do
c......2nd iteration
c......For higher accuracy determination of E and P.

       do i=0,j-2,2
         Enew(i+1) = -((Pold(i+1)+Pnew(i+1))/2)*(Vnew(i+1)
     +               -Vold(i+1))+Eold(i+1)
       end do
       do i=0,J-2,2
         Told(i+1) = Enew(i+1)/Ch
         Pnew(i+1) = comp*Told(i+1)*(1/Vnew(i+1))
     +               + (.33333*Arad*Told(i+1)**4)
         Cad(i+1) = sqrt(gamma*Pnew(i+1)/(1/Vnew(i+1)))
       end do

c......Third iteration

       do i=0,j-2,2
         Enew(i+1) = -((Pold(i+1)+Pnew(i+1))/2)*(Vnew(i+1)
     +               -Vold(i+1))+Eold(i+1)
       end do

c......Convert new energies into a temperature array.
       call rootfinder (Enew,Vnew,Told,Ch,Arad,Wm)
 
       do i=0,j-2,2
         Pnew(i+1) = comp*Told(i+1)*(1/Vnew(i+1))
     +                 +(.33333*Arad*Told(i+1)**4)
         Cad(i+1) = sqrt(gamma*Pnew(i+1)/(1/Vnew(i-1)))                
       end do
            
       call RadT(Enew,Vnew,Rnew,deltaT,k)
       call Rootfinder(Enew,Vnew,Told,Ch,Arad,wm)

c       do index=0,J-2,2
c         Eitotnew = Enew(index+1) + Eitotnew
c       end do
c       timescale = ((Eitotold)*deltaT)/(Eitotold-Eitotnew)
c       write(*,*) k,timescale
101     format (1i1,4e11.4)



c..     Get the half zone radii, and find the zones bracketing 1 Au
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
            Eold(i-1) = ((Ch/wm)*Told(index+1))+(Arad
     +                  *(Told(index+1)**4)*(1/Vold(index+1)))
            Qold(i-1) = Qnew(i-1)
          endif
        end do

        Qold(J+1) = Qnew(J+1)
        Eold(J+1) = ((Ch/wm)*Told(J+1))+(Arad*(Told(J+1)**4)
     +              *(1/Vold(J+1))) 
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

       subroutine radt(Enew,Vnew,Rnew,deltaT,k)
       implicit real*8(a-h,o-z), integer*4(i-n)
   
           parameter (J=200)

           dimension Enew(J+1),Vnew(J+1),Rnew(0:J+1)
           dimension Tnew(0:J+1),Told(0:J+1)
           dimension Rho(J+1)
           dimension rKros(J),Fbar(J)
           dimension A(J),B(J),C(J),D(J),E(-1:J),F(-1:J)
           save

c..........define physical constants (Radiation constant and Specific heat.)
           Arad = 7.56e-15
           clight = 3.0e+10
           Ch = (5./2.)*(8.317e+7) 
           wm = 2. 
           Chw = Ch/wm

c..........Use the Rootfinder routine to get the temperature array from the E's
           
           call Rootfinder(Enew,Vnew,Told,Ch,Arad,wm)

c..........Assign the Temperature boundary value:
c..........The reflective boundary condition is handled implicitly.

           Told(J+1) = 20
           Tnew(J+1) = 20
           Touter = Told(J+1)

c..........Get the Rosseland mean opacities

           do i=0,J-2,2
               T=Told(i+1)
               rho(i+1)=1/Vnew(i+1)
               call Opacity (T,rho,rkappa)
               rKros(i+1) = rkappa
           end do

c..........Start the first pass through the temperature array.
      
           do i=0,J-2,2
               Rnew(i+1) = (Rnew(i) + Rnew(i+2))/2
               if (i.eq.(J-2)) then
                   Rnew(J+1) = Rnew(J) + (Rnew(J)-Rnew(J-1))
               end if
            end do

            do i=0,J-2,2

c..............Need to calculate average values for F,Kros,Rho,T
c..............(On the grid points)

               if (i.eq.J-2) then
                   rkros(i+2) = rkros(i+1)
                   rho(i+2) = .5*(1/Vnew(i+1))
                   Told(i+2) = .5*(Told(i+1)+Touter)
               else
                   rKros(i+2) = (rKros(i+1) + rKros(i+3))/2
                   Rho(i+2) = (1/Vnew(i+1) + 1/Vnew(i+3))*.5
                   Told(i+2) = (Told(i+1) + Told(i+3))*.5
               endif
c..............Determine the flux limiter Using Bodenheimer et al 1990 method
                   del_er=arad*(Told(i+3)**4-Told(i+1)**4)
     +                    /(rnew(i+3)-rnew(i+1))
                   Rlim = abs(del_er)/((arad*Told(i+2)**4)
     +                    *rho(i+2)*rkros(i+2))
                   fluxlim = (2+Rlim)/(6+3*Rlim+Rlim**2)
                   fluxlim = 1/fluxlim
c                   fluxlim = 3.

                   Fbar(i+2) = (4*arad*clight*Told(i+2)**3)/
     +                         (fluxlim*rKros(i+2)*Rho(i+2))
           end do

c..........Calculate the A,B,C, and D coefficients:

           do i=0,J-2
               A(i+1) = Rnew(i+2)**2*Fbar(i+2)*Vnew(i+1)/
     +      (Chw*Rnew(i+1)**2*(Rnew(i+2)-Rnew(i))*(Rnew(i+3)-Rnew(i+1)))

               C(i+1) = Rnew(i)**2*Fbar(i)*Vnew(i+1)/
     +      (Chw*Rnew(i+1)**2*(Rnew(i+2)-Rnew(i))*(Rnew(i+1)-Rnew(i-1)))

               B(i+1) = A(i+1) + C(i+1) + (1/deltaT)

               D(i+1) = Told(i+1)/deltaT
           end do

c..........Now get the E and F coefficients by working outwards.

           E(-1) = 1
           F(-1) = 0

           do i=0,J-2
               E(i+1)=A(i+1)/(B(i+1)-C(i+1)*E(i-1))
               F(i+1)=(D(i+1)+C(i+1)*F(i-1))/(B(i+1)-C(i+1)*E(i-1))
           end do
      
c..........Recalculate the changes to the temperature array by working inwards
c..........Also get the energies back.
c           write(*,*) Tnew(41)

           do i=J-2,0,-2
  
               Tnew(i+1) = E(i+1) * Tnew(i+3) + F(i+1)
              if (mod(k,1).eq.0) then
c               write(*,*) i+1,'E ',E(i+1),'F ',F(i+1),'Tnew',Tnew(i+1) 
              endif
               Enew(i+1)=(cH/wm)*Tnew(i+1)+arad*Vnew(i+1)*Tnew(i+1)**4
           end do      
30         format(i3,3e10.3)
           
           return 
           end




           Subroutine Rootfinder(Enew,Vnew,Told,c,a,wm)

               implicit real*8 (a-h,o-z), integer*4 (i-n)
               parameter (J=142)
               dimension Enew(J+1),Vnew(J+1),Told(0:J+1)
               save
 
c..            need to change this to Newton-Raphson or Bisection
c..            if doing serious work!

               do i=0,J-2,2
                   do T=0.,10000.,.05
                      if ((c/wm)*T+a*(t**4)*Vnew(i+1).gt.Enew(i+1)) 
     +                                         goto 10
                   end do
10                 Told(i+1) = T
               end do
         
              Return
            end 

           subroutine opacity(t,rho,k)
           implicit real*8(a-h,o-z), integer*4(i-n)
           real*8 k

c..   a very simple opacity function.


              k=(rho/1.0e-13)+(t/100.)

           return
           end

     

       Subroutine Locate(XX,N,X,J)
c......Given an array XX of length N, and given a value X, this routine returns
c......a value J such that X is between XX(J) and XX(J+1).  J = 0 or J=N is
c......returned to indicate that X is out of range.

       implicit real*8 (a-h,o-z), integer*4(i-n)
       Dimension XX(0:N)
       jl=0
       ju=n+1
10     if (ju-jl.gt.1) then
           jm = (ju+jl)/2
           if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
               jl=jm
           else
               ju=jm
           end if
        go to 10
        end if
        j = jl
        return
        end

        Subroutine Polint(xa,ya,n,x,y,dy)
c.......Given arrays xa and ya, each of length n, and given a value x,  
c.......this routine returns a value y, and an error estimate dy.  If P(x) 
c.......is the polynomial of degree N-1 such that P(xa_i) = ya, i=1,...,n, 
c.......then the returned value y=p(x)

        implicit real*8 (a-h,o-z), integer*4(i-n)
        parameter (nmax=10)
        dimension xa(n),ya(n),c(nmax),d(nmax)
        ns = 1
        dif = abs(x-xa(1))
        do 11 i=1,n
            dift=abs(x-xa(i))
            if (dift.lt.dif) then
                ns = i
                dif = dift
            end if
            c(i)=ya(i)
            d(i)=ya(i)
11      continue
        y=ya(ns)
        ns=ns-1
        do 13 m=1,n-1
            do 12 i=1,n-m
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den = ho - hp
                if (den.eq.0.) pause
                den = w/den
                d(i) = hp*den
                c(i) = ho*den
12          continue
        if (2*ns.lt.n-m) then
            dy = c(ns+1)
        else
            dy = d(ns)
            ns = ns-1
        endif
        y = y+dy
13      continue
        return
        end

