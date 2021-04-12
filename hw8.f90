	program PipeFlow
c Taha Yaşar Demir / 1881978
c CE-580 HomeWrok #8
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err
 	real Tol
 	call Init
 	call Prior
 	Tol =  1.
c 	do while(Tol.gt.err)
	do k=1,10
 		call Solution
 		call Update
	enddo
	open(11,file='test.dat')
	write(11,*) 0,y(1)
	do i=1,N
		write(11,*) u(i),yc(i)
	enddo
	close(11)

 	stop
 	end
c-----------------------------------------------------------------------
	subroutine Init
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err

 	Um  = 4.   ! m/s
 	vis = 1e-6 ! m^2/s
 	G   = 0.95   ! Grid Ratio 1 - 0.95 - 0.9 - 0.86
 	Rad = 0.1  ! m
 	rho = 1000.! kg/m^3
 	N   = 20   ! Grid points
 	err = 1e-5 ! Error criteria
 	Ap  = 26.
 	Pi  = 22./7.
 	rkp = 0.41
 	beta= 5.3

 	call MakeGrid
 	call Analytic

 	return
 	end

c-----------------------------------------------------------------------
	subroutine MakeGrid
 	parameter(mx=100)
 	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
      real sum

	sum = 0.0
	do i=0,N-1
	 	sum = sum + G**i
	enddo
	dy(N+1) = Rad/sum
	y(N+1)  = Rad

	do i=N+1,2,-1
	  	y(i-1)  = y(i)-dy(i)
	  	dy(i-1) = G*dy(i)
	  	yc(i-1) = (y(i)+y(i-1))/2.
		r(i)    = Rad - y(i)
		rc(i-1) = Rad - yc(i-1)
c		print*, 'grid', y(i-1),dy(i-1),yc(i-1),r(i),rc(i-1)
	enddo
	do i=2,N
		dm(i) = (dy(i)+dy(i-1))/2
	enddo
	r(1) = Rad
	y(1) = 0.

 	return
 	end

c-----------------------------------------------------------------------
	subroutine Analytic
	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)

	do i=1,N
	   u(i) = Um*(yc(i)/Rad)**(1./7.)
	enddo


 	return
 	end
c-----------------------------------------------------------------------
 	subroutine Prior
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err

c 	Disc = 0.
c 	do i=1,N
c	   Disc = Disc + 2.*pi*rc(i)*u(i)*dy(i)
c	enddo
c	V_ave = Disc/(pi*Rad**2)
	V_ave = 0.85*Um
	Re    = 2.*V_ave*Rad/vis
	fm    = 0.25/((log10(5.74/(Re**0.9)))**2)
	fd    = fm 
	tau_w = (rho*fd*V_ave**2)/8.
c	tau_w = 27.27
	us    = sqrt(tau_w/rho)
	do i=1,N
		yp(i)   = yc(i)*us/vis
		fu(i)   = 1-exp(-yp(i)/Ap)
		fml(i)  = Rad*(0.14-0.08*(1-(yc(i)/Rad))**2
     &                         -0.06*(1-(yc(i)/Rad))**4)*fu(i)
	enddo
	do i=2,N
		vist(i) = (fml(i)**2)*abs((u(i)-u(i-1))/dm(i))
		vise(i) = (vis + vist(i))
	enddo

	return
	end
c-----------------------------------------------------------------------
	subroutine Solution
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err

 	Cp = -2.*(tau_w/Rad)/rho
 	Cw = tau_w/u(1)
 	do i=1,N
 		cd(i) = dy(i)*rc(i)
 	enddo
 	do i=2,N
 		ca(i) = r(i)/dm(i)
 	enddo

c 	rtau(N+1) = 0.
c 	do i=N,2,-1
c 		rtau(i) = rtau(i+1) - dy(i)*rc(i)*Cp 
c 	enddo
c 	u(N) = 4.
c 	do i=N,3,-1
c 		u(i-1)= u(i) - (rtau(i)*dm(i))/(r(i)*vise(i))
c 		print*, 'deneme', u(i-1),rtau(i),vise(i-1)
c 	enddo
c 	u(1) = ( (u(2) - (Cp*cd(1)+Rad*tau_w)/(ca(2)*vise(2))))
c 	print*,'u1', u(1) ,'tau_w',tau_w,'Cp',Cp

 	do i=1,N
	if(i.eq.1) then
 		a(i) = 0.
 		c(i) =-ca(i+1)*vise(i+1)
 		b(i) =-c(i) + Cw*Rad
 		d(i) =-cd(i)*Cp
 	elseif(i.eq.N) then
		a(i) =-ca(i)*vise(i)
 		c(i) = 0.
 		b(i) =-a(i)
 		d(i) =-cd(i)*Cp 
 	else
 		a(i) =-ca(i)*vise(i)
 		c(i) =-ca(i+1)*vise(i+1)
 		b(i) =-(a(i)+c(i))
 		d(i) =-cd(i)*Cp
	endif
	enddo
c      call TRID (a,b,c,d,u,N,0,0,1,4.)
	call THOMAS(1,N,a,b,c,d)
	do i=1,N
		u(i) = d(i)
		print*, 'solution', u(i),tau_w,Cp,fml(i)
	enddo
	return
 	end

c-----------------------------------------------------------------------
	subroutine Update
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err
 	real law, us_old

c 	law  = (1/rkp)*log(yp(1))+beta
c 	us   = u(1)/(law)
	us_old = us
	call ustar
	us_old = us
 	tau_w = (us**2)*rho
 	tol = abs(tau_old-tau_w)/tau_old
 	print*, 'Error',tol,tau_old,us
	do i=1,N
		yp(i)   = yc(i)*us/vis
		fu(i)   = 1-exp(-yp(i)/Ap)
		fml(i)  = Rad*(0.14-0.08*(1-(yc(i)/Rad))**2
     &                         -0.06*(1-(yc(i)/Rad))**4)*fu(i)
	enddo
	do i=2,N
c		vist(i) =(vist(i)+((fml(i)**2)*abs((u(i)-u(i-1))/dm(i))))/2.
		vist(i) =(((fml(i)**2)*abs((u(i)-u(i-1))/dm(i))))
		vise(i) = vis + vist(i)
c		print*, 'eff', vise(i)
	enddo

	Disc = 0.
 	do i=1,N
	   Disc = Disc + u(i)*dy(i)
	enddo
	V_ave = Disc/(Rad)
	Re    = 2.*V_ave*Rad/vis
	fm    = 0.25/((log10(5.74/(Re**0.9)))**2)
	fd    = (8.*tau_w)/(rho*(V_ave**2))
c	print*, fm,fd,tau_w,V_ave,Pi,u(1)

 	return
 	end
c-----------------------------------------------------------------------
	subroutine ustar
	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err
 	real test, yplus
 	yplus = us*yc(1)/vis
 	test  = us*(1/rkp)*log(yplus) + us*beta
 	if(test.lt.u(1)) then 
 		do while(test.lt.u(1))
 			us   = 1.0001*us
 			yplus = us*yc(1)/vis
 			test = us*(1/rkp)*log(yplus) + beta*us
 		enddo
 	else
		do while(test.gt.u(1))
 			us = 0.9999*us
 			yplus = us*yc(1)/vis
 			test  = us*(1/rkp)*log(yplus) + beta*us
 		enddo
 	endif

	return
	end

c-----------------------------------------------------------------------
      subroutine THOMAS(il,iu,aa,bb,cc,ff)
c............................................................
c  Solution of a tridiagonal system of n equations of the form
c  A(i)*x(i-1) + B(i)*x(i) + C(i)*x(i+1) = R(i)  for i=il,iu
c  the solution X(i) is stored in F(i)
c  A(il-1) and C(iu+1) are not used
c  A,Bb,C,R are arrays to bbe provided bby the user
c............................................................
      parameter (mx=100)
      dimension aa(mx),bb(mx),cc(mx),ff(mx),tmp(mx)

      tmp(il)=cc(il)/bb(il)
      ff(il)=ff(il)/bb(il)
      ilp1 = il+1
      do i=ilp1,iu
         z=1./(bb(i)-aa(i)*tmp(i-1))
         tmp(i)=cc(i)*Z
         ff(i)=(ff(i)-aa(i)*ff(i-1))*z
      enddo
      iupil=iu+il
      do ii=ilp1,iu
         i=iupil-ii
         ff(i)=ff(i)-tmp(i)*ff(i+1)
      enddo
      return
      end

      SUBROUTINE TRID (aa,bb,cc,dd,S,M,L1,B1,LM,BM)
      DIMENSION E(5001),F(5001),S(M),aa(M),bb(M),cc(M),dd(M)
      IF(L1.EQ.1) E(1)=0.
      IF(L1.EQ.1) F(1)=B1
      IF(L1.EQ.2) E(1)=1.
      IF(L1.EQ.2) F(1)=-B1
      MM1=M-1
      DO I=2,MM1
    		DEN=bb(I)-cc(I)*E(I-1)
    		E(I)=aa(I)/DEN
    		F(I)=(dd(I)+cc(I)*F(I-1))/DEN
	END DO
	IF(LM.EQ.1) S(M)=BM
    	IF(LM.EQ.2) S(M)=(F(MM1)+BM)/(1.-E(MM1))
    	DO I=MM1,1,-1
    	S(I)=E(I)*S(I+1)+F(I)
	END DO
    	RETURN
    	END