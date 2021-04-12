	program PipeFlow
c Taha Yaşar Demir / 1881978
c CE-580 HomeWrok #8
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol
 	real Tol
 	call Init
 	call Prior
 	err =  1.
 	do while(err.gt.Tol)
c	do k=1,10
 		call Solution
 		call Update
	enddo
c	open(11,file='test.dat')
c	write(11,*) 0,y(1)
c	do i=1,N
c		write(11,*) u(i),yc(i)
c	enddo
c	close(11)

 	stop
 	end
c-----------------------------------------------------------------------
	subroutine Init
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol

 	Um  = 4.   ! m/s
 	vis = 1e-6 ! m^2/s
 	G   = 1.   ! Grid Ratio 1 - 0.95 - 0.9 - 0.86
 	Rad = 0.1  ! m
 	rho = 1000.! kg/m^3
 	N   = 20   ! Grid points
 	Tol = 1e-5 ! Error criteria
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
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol

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
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol

 	Cp = -2.*(tau_w/Rad)/rho
 	Cw = tau_w/u(1)
 	do i=1,N
 		cd(i) = dy(i)*rc(i)
 	enddo
 	do i=2,N
 		ca(i) = r(i)/dm(i)
 	enddo

 	do i=1,N
	if(i.eq.1) then
 		c(i) = 0.
 		a(i) = ca(i+1)*vise(i+1)
 		b(i) = a(i) + Cw*Rad
 		d(i) =-cd(i)*Cp
 	elseif(i.eq.N) then
		c(i) = ca(i)*vise(i)
 		a(i) = 0.
 		b(i) = c(i)
 		d(i) =-cd(i)*Cp 
 	else
 		c(i) = ca(i)*vise(i)
 		a(i) = ca(i+1)*vise(i+1)
 		b(i) = (a(i)+c(i))
 		d(i) =-cd(i)*Cp
	endif
	enddo
	! L1 and B1 is zero, LM , BM would be 2,0. but that does not work
      call TRID (a,b,c,d,u,N,0,0.,1,4.)
	
	do i=1,N
		print*, 'solution', u(i)
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
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol
 	real law, us_old
      ! Update the mixing-length parameters for next calculation
	us_old = us
	call ustar
 	tau_w = (us**2)*rho
 	err = abs(us_old - us)/us_old
 	print*, 'Error',err,us_old,us
	do i=1,N
		yp(i)   = yc(i)*us/vis
		fu(i)   = 1-exp(-yp(i)/Ap)
		fml(i)  = Rad*(0.14-0.08*(1-(yc(i)/Rad))**2
     &                         -0.06*(1-(yc(i)/Rad))**4)*fu(i)
	enddo
	do i=2,N
		vist(i) =(((fml(i)**2)*abs((u(i)-u(i-1))/dm(i))))
		vise(i) = vis + vist(i)
	enddo

	Disc = 0.
 	do i=1,N
	   Disc = Disc + u(i)*dy(i)
	enddo
	V_ave = Disc/(Rad)
	Re    = 2.*V_ave*Rad/vis
	fm    = 0.25/((log10(5.74/(Re**0.9)))**2)
	fd    = (8.*tau_w)/(rho*(V_ave**2))

 	return
 	end
c-----------------------------------------------------------------------
	subroutine ustar
	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol
 	real test, yplus
 	! compares loglaw left and right sides and adjusts us for next calculations
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