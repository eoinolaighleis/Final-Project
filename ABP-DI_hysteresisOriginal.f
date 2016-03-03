c  ABP-DI model with hysteresis and Verlet scheme
c  (c) Maksym Romenskyy
c  19 August 2014 - latest version
       
       PROGRAM SOCIAL 
         implicit none
	 
       INTEGER BOX,BOX2,L,ST,ip,ipp,np,PIPE
       INTEGER nlist(100000),list(50000,10000),timer,J2,hystt,hystp
       
       LOGICAL WALL,TUBE,xp,yp,TRANSVERSE,MOTOR,LANGEVIN,ONLYBOX
       
       REAL P,Vold,Vnew,pi,Vxx,Vyy,V(3,100000),TEMPERATURE_dpd
       REAL DX,DY,DX2,DY2,dt,rand,TOR,Tcross,Tnew,Force(3,100000)
       REAL Lin,Rrepul,GAMMAl,NOI,Rcut
       REAL cosfall,fall,Rux,Ruy,WD,WR,hysts,Prx,Pry
       REAL Rij,Rijx,Rijy,Field,VCx,VCy,gausx,gausy
       REAL gamma,DotProduct,sigma,Gnoise,Ar,TEMPERATURE
       REAL Fdisx(100000),Fdisy(100000),Fconsx(100000),Fconsy(100000)
       REAL Frandx(100000),Frandy(100000),Ftotx(100000),Ftoty(100000)
       REAL d2,q0,c1,d2q0,DEPOTx,DEPOTy,Flx,Fly  
       REAL R(3,100000),Rb(3,100000),Rtempx,Rtempy,R2(3,100000)
       
       CHARACTER(len=80)::CoordFile,HystFile
       CHARACTER(LEN=4)::NAME
  
!--------------------------Filenames----------------------------------!       
       
       	NAME="test"
       
	write(CoordFile,'(a,a)')name,"_coord.dat"
	write(HystFile,'(a,a)')name,"_hyst.dat"
	    	
!-----------------------Key parameters--------------------------------!
	
	BOX=500            	!       Box size	!
	np=1000           	!  Number of particles	!
	ST=1000            	!    Number of steps	!
	dt=0.01                 !        timsetep       !
	hysts=1.0		!     Field strength    !
	hystp=100               !        Period         !
	
!-----------------------DPD parameters--------------------------------!
	
	TEMPERATURE_dpd=0.0     !    Temperature DPD	!
	GAMMA=1.0        	!     Friction DPD	!
	
	Rcut=3.0          	!     Radius of cut	!
	Rrepul=rcut/2           !    Repulsion radius	!
	Ar=1.0                  !   Repulsion strength	!
	
	TRANSVERSE=.false.

!---------------------Langevin parameters------------------------------!		
	
	GAMMAl=0.3		!     Langevin gamma	!
	TEMPERATURE=0.3      	!  Langevin temperature !
	LANGEVIN=.true.
	MOTOR=.true.
	
	d2=2.0			! Conversion of internal energy into kinetic
	q0=1.0			! Gain of energy from the environment
	c1=0.8			! Energy loss
	
!-------------------------Confinement---------------------------------!

	TUBE=.true.
	PIPE=50     		!      Pipe width	!
	
!----------------------Other parameters-------------------------------!
	
	pi=4.0d0*ATAN(1.0d0)
	P=BOX
	Vold=1.0
	call init_random_seed()
	
!----------------------Pre-installations------------------------------!
	
	hystp=hystp/dt
	BOX2=BOX*0.5
	SIGMA=SQRT(2.0*GAMMA*TEMPERATURE_dpd/dt)
	NOI=SQRT(2.0*GAMMAl*TEMPERATURE/dt)
	d2q0=d2*q0
	
	timer=0
	hystt=0
       
       IF(TRANSVERSE)THEN
	write(*,*)"***TRANSVERSE DPD***"
       END IF
       
       IF(TUBE)THEN
       write(*,*)"INFO) Modelling with walls"
       xp=.true.
       yp=.false.
       ELSE
       xp=.true.
       yp=.true.
       ONLYBOX=.true.
       END IF
       
       write(*,*)" "
	
       OPEN(22,FILE=CoordFile)
       OPEN(82,FILE=HystFile)
       
       write(*,*)" "
!-------------------Position & Velocity Setup-------------------------!       	         		             
       
       DO ip=1,np
       	IF(TUBE)THEN
	R(1,ip)=rand(0)*BOX+1 
	R(2,ip)=rand(0)*((PIPE-2)+1-2)+2
	ELSE
	R(1,ip)=rand(0)*BOX+1 
	R(2,ip)=rand(0)*BOX+1
	END IF
	DX=2*(rand(0)-0.5)                              
       	DY=2*(rand(0)-0.5)              
       	DX2=DX**2
       	DY2=DY**2
       	Vnew=1/SQRT(DX2+DY2)
       	V(1,ip)=DX*Vnew*Vold
	V(2,ip)=DY*Vnew*Vold
	Rb(1,ip)=R(1,ip)
	Rb(2,ip)=R(2,ip)
       END DO
       write(*,*)"Starting integration..."
       do ip=1,np
       R2(1,ip)=R(1,ip)
       R2(2,ip)=R(2,ip)
       enddo
       
!------------------------Main Loop------------------------------------!
       DO L=1,ST
       
       IF(hystt.lt.hystp)THEN
       hystt=hystt+1
       ELSE IF(hystt.eq.hystp)THEN
       hystt=1
       END IF
       Field=hysts*sin((2*pi)/hystp*hystt)
       
	TOR=0.0
       
       DO ip=1,np       
	
!- - - - - - - - - - - - - - WALLS - - - - - - - - - - - - - -!
	
	IF(TUBE)THEN
	WALL=.false.
	
	IF(V(2,ip).gt.0.0)THEN
	Tcross=ABS((PIPE-R(2,ip))/(V(2,ip)))
	IF(Tcross<=dt)THEN
	!WALL=.true.
	Tnew=dt-Tcross
	cosfall=V(2,ip)/sqrt(V(1,ip)**2+V(2,ip)**2)
	fall=pi-2*acos(cosfall)
	IF(V(1,ip).gt.0.0)THEN	
	Vxx=V(1,ip)*cos(fall)+V(2,ip)*sin(fall)
       	Vyy=V(2,ip)*cos(fall)-V(1,ip)*sin(fall)
	ELSE IF(V(1,ip).lt.0.0)THEN
	Vxx=V(1,ip)*cos(fall)-V(2,ip)*sin(fall)
       	Vyy=V(2,ip)*cos(fall)+V(1,ip)*sin(fall)
	END IF
	END IF
	END IF
	
	IF(V(2,ip).lt.0.0)THEN
	Tcross=ABS((0.0-R(2,ip))/(V(2,ip)))
	IF(Tcross<=dt)THEN
	!WALL=.true.
	Tnew=dt-Tcross
	cosfall=(-V(2,ip))/sqrt(V(1,ip)**2+V(2,ip)**2)
	fall=pi-2*acos(cosfall)
	IF(V(1,ip).gt.0.0)THEN
	Vxx=V(1,ip)*cos(fall)-V(2,ip)*sin(fall)
       	Vyy=V(2,ip)*cos(fall)+V(1,ip)*sin(fall)
	ELSE IF(V(1,ip).lt.0.0)THEN
	Vxx=V(1,ip)*cos(fall)+V(2,ip)*sin(fall)
       	Vyy=V(2,ip)*cos(fall)-V(1,ip)*sin(fall)
	END IF
	END IF
	END IF
	
	IF(WALL)THEN
	IF(V(1,ip).eq.0.0)THEN
	V(2,ip)=-V(2,ip)
	Rtempx=R(1,ip)+V(1,ip)*Tnew
	Rtempy=R(2,ip)+V(2,ip)*Tnew
	!Rtempx=2*R(1,ip)-Rb(1,ip)+Force(1,ip)*Tnew**2
       	!Rtempy=2*R(2,ip)-Rb(2,ip)+Force(2,ip)*Tnew**2
       	!V(1,ip)=(Rtempx-R(1,ip))/Tnew
       	!V(2,ip)=(Rtempy-R(2,ip))/Tnew
       	Rb(1,ip)=R(1,ip)
       	Rb(2,ip)=R(2,ip)
       	R(1,ip)=Rtempx
       	R(2,ip)=Rtempy
	R2(1,ip)=MODULO(R(1,ip),P)
	R2(2,ip)=R(2,ip)
	
	ELSE
	
	V(1,ip)=Vxx
	V(2,ip)=Vyy
	Rtempx=R(1,ip)+V(1,ip)*Tnew
	Rtempy=R(2,ip)+V(2,ip)*Tnew
	!Rtempx=2*R(1,ip)-Rb(1,ip)+Vxx*Tnew**2
       	!Rtempy=2*R(2,ip)-Rb(2,ip)+Vyy*Tnew**2
       	!V(1,ip)=(Rtempx-R(1,ip))/Tnew
       	!V(2,ip)=(Rtempy-R(2,ip))/Tnew
       	Rb(1,ip)=R(1,ip)
       	Rb(2,ip)=R(2,ip)
       	R(1,ip)=Rtempx
       	R(2,ip)=Rtempy
       	R2(1,ip)=MODULO(R(1,ip),P)
	R2(2,ip)=R(2,ip)
	END IF
	
	ELSE!(no WALL)
	
	Rtempx=2*R(1,ip)-Rb(1,ip)+Force(1,ip)*dt**2
       	Rtempy=2*R(2,ip)-Rb(2,ip)+Force(2,ip)*dt**2
	
	!if((Rtempy.gt.pipe).or.(Rtempy.lt.0.0))then
	!Rtempy=Rb(2,ip)
	!end if
       	V(1,ip)=(Rtempx-R(1,ip))/dt
       	V(2,ip)=(Rtempy-R(2,ip))/dt
       	Rb(1,ip)=R(1,ip)
       	Rb(2,ip)=R(2,ip)
       	R(1,ip)=Rtempx
       	R(2,ip)=Rtempy
	R2(1,ip)=MODULO(R(1,ip),P)
	R2(2,ip)=MODULO(R(2,ip),50.0)
	END IF
	END IF
	
! - - - End wall
	
	WRITE(22,*)L,ip,R(1,ip),R(2,ip)

       END DO
       
    !------------------DPD-------------------!
       timer=timer+1
       if(timer==8) timer=0
       if(timer==1)then
       call new_vlist(R2,np,list,nlist,box,box2,xp,yp)
       end if
       
       DO ip=1,np
       Fdisx(ip)=0
       Fdisy(ip)=0
       Fconsx(ip)=0
       Fconsy(ip)=0
       Frandx(ip)=0
       Frandy(ip)=0
       Ftotx(ip)=0
       Ftoty(ip)=0
       END DO
              
       DO ip=1,np
       DO ipp=1,nlist(ip)

! ----> Calculate distance      
       
       J2=list(ip,ipp)
       IF(J2.le.ip)CYCLE
        Rijx=R2(1,ip)-R2(1,J2)
        Rijy=R2(2,ip)-R2(2,J2)
       
       IF(xp)THEN
       if(Rijx.gt.box2)then
	Rijx=Rijx-box
       else if(Rijx.lt.-box2)then
	Rijx=Rijx+box
       end if
       END IF
       
       IF(yp)THEN
       if(Rijy.gt.box2)then
	Rijy=Rijy-box
       else if(Rijy.lt.-box2)then
	Rijy=Rijy+box
       end if
       END IF
       
        Rij=SQRT(Rijx**2+Rijy**2)
       
       IF(Rij.lt.Rcut)THEN
       
! ----> Define unit vectors
      	
	IF(Rijx.ne.0.0)THEN
	Rux=Rijx/Rij
	ELSE
	Rux=0.0
	END IF
	IF(Rijy.ne.0.0)THEN
      	Ruy=Rijy/Rij
	ELSE
	Ruy=0.0
	END IF

! ----> Define weight functions
	
	WR=1-Rij/Rcut
	WD=WR**2
		
! ----> Define forces	
	
! ----> Dissipative

	IF(TRANSVERSE)THEN ! Transverse DPD
	VCx=((1-Rux*Rux)*(V(1,ip)-V(1,J2))-Rux*Ruy*(V(2,ip)-V(2,J2)))
	VCy=(-(Ruy*Rux)*(V(1,ip)-V(1,J2))+(1-Ruy*Ruy)*(V(2,ip)-V(2,J2)))
	Fdisx(ip)=-gamma*WD*VCx
     	Fdisy(ip)=-gamma*WD*VCy
	Fdisx(J2)=-Fdisx(ip)
	Fdisy(J2)=-Fdisy(ip)
	ELSE		   ! Traditional DPD
      	DotProduct=(V(1,ip)-V(1,J2))*Rux+(V(2,ip)-V(2,J2))*Ruy
       	Fdisx(ip)=-gamma*WD*DotProduct*Rux
	Fdisy(ip)=-gamma*WD*DotProduct*Ruy
	Fdisx(J2)=-Fdisx(ip)
	Fdisy(J2)=-Fdisy(ip)
	END IF
	
! ----> Conservative

	IF(Rij.lt.Rrepul)THEN
	Lin=1-Rij/Rrepul
	Fconsx(ip)=Ar*Lin*Rux
	Fconsy(ip)=Ar*Lin*Ruy
	Fconsx(J2)=-Fconsx(ip)
	Fconsy(J2)=-Fconsy(ip)
	ELSE
	Fconsx(ip)=0.0
	Fconsy(ip)=0.0
	Fconsx(J2)=0.0
	Fconsy(J2)=0.0
	END IF

! ----> Random
	
	IF(TEMPERATURE_dpd.ne.0.0)THEN
	gausx=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
	gausy=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
	IF(TRANSVERSE)THEN
	Prx=(1-Rux*Rux)-Rux*Ruy
	Pry=-Ruy*Rux+(1-Ruy*Ruy)
	Frandx(ip)=sigma*WR*gausx*Prx
	Frandy(ip)=sigma*WR*gausy*Pry
	Frandx(J2)=-Frandx(ip)
	Frandy(J2)=-Frandy(ip)
	ELSE
	Frandx(ip)=sigma*WR*gausx*Rux
	Frandy(ip)=sigma*WR*gausy*Ruy
	Frandx(J2)=-Frandx(ip)
	Frandy(J2)=-Frandy(ip)
	END IF
	END IF

! ----> Total "collective" force
	
	Ftotx(ip)=Ftotx(ip)+(Fdisx(ip)+Fconsx(ip)+Frandx(ip))
	Ftoty(ip)=Ftoty(ip)+(Fdisy(ip)+Fconsy(ip)+Frandy(ip))
	Ftotx(J2)=Ftotx(J2)+(Fdisx(J2)+Fconsx(J2)+Frandx(J2))
	Ftoty(J2)=Ftoty(J2)+(Fdisy(J2)+Fconsy(J2)+Frandy(J2))
       END IF
       END DO
       END DO
       
       DO ip=1,np

! ----> Internal energy depot
       
       IF(MOTOR)THEN
       	DEPOTx=d2q0/(c1+d2*(V(1,ip)**2+V(2,ip)**2))
       	DEPOTy=d2q0/(c1+d2*(V(1,ip)**2+V(2,ip)**2))
       ELSE
       	DEPOTx=0.0
       	DEPOTy=0.0
       END IF
       
       IF(LANGEVIN)THEN
       Gnoise=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
       Flx=-(GAMMAl-DEPOTx)*V(1,ip)+NOI*Gnoise
       Gnoise=sqrt(-2.0*log(rand(0)))*cos(2*pi*rand(0))
       Fly=-(GAMMAl-DEPOTy)*V(2,ip)+NOI*Gnoise
       ELSE
       Flx=0.0
       Fly=0.0
       END IF  

! ----> Adding "collective" force to ABP force
       
       Force(1,ip)=Ftotx(ip)+Flx
       Force(2,ip)=Ftoty(ip)+Fly
       
       END DO
              
       DO ip=1,np
       
! - - - Adding acceleration along x-axis
       
       Force(1,ip)=Force(1,ip)+Field

! - - - Integrator for periodic box only
       
	IF(ONLYBOX)THEN
	Rtempx=2*R(1,ip)-Rb(1,ip)+Force(1,ip)*dt**2
       	Rtempy=2*R(2,ip)-Rb(2,ip)+Force(2,ip)*dt**2
       	V(1,ip)=(Rtempx-R(1,ip))/dt
       	V(2,ip)=(Rtempy-R(2,ip))/dt
       	Rb(1,ip)=R(1,ip)
       	Rb(2,ip)=R(2,ip)
       	R(1,ip)=Rtempx
       	R(2,ip)=Rtempy
       	R2(1,ip)=MODULO((R(1,ip)),P)
       	R2(2,ip)=MODULO((R(2,ip)),P)
	END IF
	
	TOR=TOR+V(1,ip)
       
       END DO                             
       
	TOR=TOR/np
	IF(L.gt.hystp)THEN  
	write(82,*)Field,TOR
	END IF
       
       
      END DO                                             !--END OF THE MAIN LOOP--!
       
       write(*,*)"Integration done"
       
       CLOSE(22)
       CLOSE(82)
       
       WRITE(*,*) "DONE:",name  
           
      END PROGRAM SOCIAL
      
!- - - - - - - - - - - - VERLET LIST - - - - - - - - - - - -!      
      
      SUBROUTINE new_vlist(R2,npart,list,nlist,box,box2,xp,yp)               
      Implicit none
      
      INTEGER ii,jj,nlist(100000),list(50000,10000),npart
      INTEGER box,box2,VERLET
      REAL Rx,Ry,xr,R2(3,100000)
      LOGICAL xp,yp
      
      VERLET=5**2
      
      do ii=1,npart
      	nlist(ii)=0      	
      end do
      do ii=1,npart-1
      do jj=ii+1,npart
      	Rx=R2(1,ii)-R2(1,jj)
	Ry=R2(2,ii)-R2(2,jj)
      if(xp)then
       if(Rx.gt.box2)then
	Rx=Rx-box
       else if(Rx.lt.-box2)then
	Rx=Rx+box
       end if
       end if
       if(yp)then
       if(Ry.gt.box2)then
	Ry=Ry-box
       else if(Ry.lt.-box2)then
	Ry=Ry+box
       end if
       end if
	xr=Rx**2+Ry**2	      	
      if(abs(xr).lt.VERLET)then
      	nlist(ii)=nlist(ii)+1
      	nlist(jj)=nlist(jj)+1
      	list(ii,nlist(ii))=jj
      	list(jj,nlist(jj))=ii
      end if
      end do
      end do
      return
      END

!- - - - - - - - - RANDOM SEED INITIALISATION - - - - - - - - -!
     
      SUBROUTINE init_random_seed()
            INTEGER*8 :: i, n, clock,seedi
            INTEGER*8, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed=abs((clock+37*(/ (i - 1, i = 1, n) /))**2)    
            CALL RANDOM_SEED(PUT = seed)
          	seedi=rand(seed)
            DEALLOCATE(seed)
      END SUBROUTINE
