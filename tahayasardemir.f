	program PipeFlow
c Taha Yaşar Demir / 1881978
c CE-580 HomeWrok #8
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count
 
 	call Init
 	call Prior
 	open(11,file="yplus.dat")
 	open(12,file="output.dat")
 	err =  1.
 	count= 0.
 	do while(err.gt.Tol)
c	do k=1,100
 		call Solution
 		count = count +1
 		call Update
 		print*, 'iteration number', count , 'Error' , err
 		call Output
	enddo
	close(11)
	close(12)

 	stop
 	end
c-----------------------------------------------------------------------
	subroutine Init
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count

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
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count

 	Disc = 0.
 	do i=1,N
	   Disc = Disc + 2.*pi*rc(i)*u(i)*dy(i)
	enddo
	V_ave = Disc/(pi*Rad**2)
c	V_ave = 0.85*Um
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
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count

 	Cp = -2.*(tau_w/Rad)/rho
 	Cw = tau_w/u(1)
 	do i=1,N
 		cd(i) = dy(i)*rc(i)
 	enddo
 	do i=2,N
 		ca(i) = r(i)/dm(i)
 	enddo

 	rtau(N+1) = 0.
 	do i=N,2,-1
 		rtau(i) = rtau(i+1) - dy(i)*rc(i)*Cp 
 	enddo
c 	u(N) = 4.
 	do i=N,3,-1
 		u(i-1)= u(i) - (rtau(i)*dm(i))/(r(i)*vise(i))
 	enddo
 	u(1) = u(2) - (rtau(2)*dm(2))/(r(2)*vise(2))/r(2) - tau_w/rho
 	print*,'u1', u(1) ,'tau_w',tau_w,'Cp',Cp, 'uN', u(N)

	return
 	end

c-----------------------------------------------------------------------
	subroutine Update
 	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count
 	real law, us_old

	us_old = us
	call ustar
	tau_w = (us**2)*rho
 	err = abs(us_old-us)/us_old
 	print*, 'error',err,us
	do i=1,N
		yp(i)   = yc(i)*us/vis
		fu(i)   = 1-exp(-yp(i)/Ap)
		fml(i)  = Rad*(0.14-0.08*(1-(yc(i)/Rad))**2
     &                         -0.06*(1-(yc(i)/Rad))**4)*fu(i)
	enddo
	do i=2,N
		! Take average of two consequent viscous stress to eliminate oscilation
		vist(i) =(vist(i)+((fml(i)**2)*abs((u(i)-u(i-1))/dm(i))))/2.
c		vist(i) =(((fml(i)**2)*abs((u(i)-u(i-1))/dm(i))))
		vise(i) = vis + vist(i)
	enddo

 	return
 	end
c-----------------------------------------------------------------------
	subroutine ustar
	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count
 	real test, yplus
 	yplus = us*yc(1)/vis
 	test  = us*(1/rkp)*log(yplus) + us*beta
 	if(test.lt.u(1)) then 
 			us   = us+ 0.005*us/(10*count) 
 	else

 			us = (0.995*us+us)/2
 	endif

	return
	end
c-----------------------------------------------------------------------
	subroutine Output
	parameter(mx=100)
	common/grid/ Rad,r(mx),rc(mx),y(mx),yc(mx),dy(mx),dm(mx),N,G
 	common/var/  tau_w,u(mx),vis,vist(mx),vise(mx),rho,Um,rtau(mx)
 	common/const/Cp,rkp,beta,Cw,ca(mx),cd(mx),a(mx),b(mx),c(mx),d(mx)
 	common/turb/ fml(mx),fu(mx),yp(mx),Ap,us,up(mx)
 	common/out/  fm,fd,V_ave,Disc,Re,Pi,err,Tol,count
 	real loglaw,relative
 	do i=1,N
	   Disc = Disc + 2.*pi*rc(i)*u(i)*dy(i)
	enddo
	V_ave = Disc/(pi*Rad**2)
c	V_ave = 0.85*Um
	Re    = 2.*V_ave*Rad/vis
	fm    = 0.25/((log10(5.74/(Re**0.9)))**2)
	fd    = (8*tau_w)/(rho*(V_ave**2))
	if(err.lt.Tol) then
	do i=1,N
		yp(i) = us*yc(i)/vis
		up    = u(i)/us
		loglaw= (1./rkp)*log(yp(i))+beta
		write(11,*) yp(i),up(i),loglaw
	enddo
	else
	relative = 100*abs(fd-fm)/fm
	write(12,*) G,yp(1),fd,relative,count
	print*,'Ratio',' yp1',' fc',' difference',' Iteration'
	print*, G,yp(1),fd,relative,count
	endif
	return
	end