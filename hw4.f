	  program TurbulentPipeFlow
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  common/error/u_old(mx),rel_err,iteration


	  open(11,file='yplus.dat')
	  open(22,file='error.dat')
	  open(33,file='difno.dat')
	  open(44,file='turbv.dat')
	  open(55,file='loglaw.dat')

	  iteration = 10000000

	  call init
	  u(1) = 0.
	  u(N) = u_max
	  do i=1,iteration
	  	call stress
	  	call coefficients
	  	call update
	  	if (i.gt.1) then
	  		call error_calc
	  		call output(i)
	  	endif
	  enddo
	  do i=2,N
	  print*, vise(i),tau,Cp,u(i),y(i),us,sml(i)
	  enddo

	  close(11)
	  close(22)
	  close(33)
	  close(44)
	  close(55)

	  stop
	  end

c-----------------------------------------------------------------------
	  subroutine init
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  common/error/u_old(mx),rel_err,iteration

	  u_max = 4.
	  vis   = 1E-6
	  Beta  = 0.82
	  rho   = 1000.
	  Rad   = 0.1
	  Beta  = 0.82
	  N     = mx
	  dT    = 1E-4

	  call makegrid 
	  call analytic


	  return
	  end
c-----------------------------------------------------------------------
	  subroutine makegrid
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  common/error/u_old(mx),rel_err,iteration


      real sum

	  sum = 0.0
	  do i=0,N-2
	  	sum = sum + Beta**i
	  enddo
	  dy(N) = Rad/sum
	  y(N)  = Rad

	  do i=N,2,-1
	  	y(i-1)  = y(i)-dy(i)
	  	dy(i-1) = Beta*dy(i)
	  	yc(i)   = (y(i)+y(i-1))/2.
	    r(i)    = Rad - y(i)
	  	rc(i)   = Rad - yc(i)
	  enddo
	  r(1) = Rad
	  y(1) = 0.


	  return
	  end
c-----------------------------------------------------------------------
	  subroutine analytic
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)

	  do i=2,N-1
	  	u(i) = u_max*(y(i)/Rad)**(1./7.)
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine stress
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau

	  tau = rho*vis*(u(2)-u(1))/dy(2)
	  us  = sqrt(tau/rho)
	  Ap  = 26.
	  do i=2,N
	  	yp(i)  = yc(i)*us/vis
	  	fu(i)  = 1-exp(-yp(i)/Ap)
	    sml(i) = Rad*(0.14-0.08*(1-(yc(i)/Rad))**2
     +           -0.06*(1-(yc(i)/Rad))**4)*fu(i)
	    vist(i)= (sml(i)**2.)*((u(i)-u(i-1))/dy(i))
	    vise(i)= vist(i) + vis
	  enddo

	  Cp  =-2.*(tau/Rad)

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine coefficients
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  real term

	  do i=2,N-1
	  	term   = (dy(i) + dy(i+1))
	  	A(i)   = ( 2.*rc(i+1)) / (r(i)*term*dy(i+1) )
	  	C(i)   = ( 2.*rc(i)) / (r(i)*term*dy(i) )
	  	eps(i) = dT*( -(Cp/rho) + A(i)*vise(i+1)*(u(i+1)-u(i))
     +                -C(i)*vise(i)*(u(i)-u(i-1))  )


	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine update
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/coef/ A(mx), C(mx), eps(mx), dT

	  do i=2,N-1
	  	u(i) = u(i) + eps(i)
	  enddo

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine error_calc
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  common/error/u_old(mx),rel_err,iteration
	  real ttest
	  rel_err = 0.0

	  do i=2,N-1
	  	rel_err = rel_err + ((1./(N-2)/u_max)*abs(eps(i)))
	  enddo

	  ttest = rel_err

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine output(iter)
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  common/error/u_old(mx),rel_err,iteration
	  real U_p(mx),law(mx),Dis,V_ave,Re,fm,fd
	  integer iter

	  write(22,*) iter,rel_err
	  if (iter.eq.iteration) then
	  	do k=2,N
	  		yp(k)  = y(k)*us/vis
	  		U_p(k) = u(k)/us
	  		write(11,*) yp(k),U_p(k)
	  		law(k) = (1./0.41)*log(yp(k))+5.1
	  		write(55,*) yp(k),law(k)
	  		if (k.lt.N) then
	  			df(k)  = (0.5*(vise(k)+vise(k+1))*dT)/(dy(k)**2)
	  			write(33,*) yp(k),df(k)
	  		endif
	  		write(44,*) vist(k)/vis, y(k)/Rad
	  	enddo
	  	call output_par(Dis,V_ave,Re,fm,fd)
	  	print*, "Discharge:Average Vel:Reynolds:fm:fd"
	  	print*,  Dis,V_ave,Re,fm,fd
	  endif

	  return
	  end
c-----------------------------------------------------------------------
	  subroutine output_par(discha,ave,reynolds,fr_m,fr_d)
	  parameter (mx=31)
	  common/flow/ rho,Rad,Cp,vis,vist(mx),vise(mx),u(mx),u_max
	  common/grid/ Beta, y(mx),yc(mx),dy(mx),r(mx),rc(mx),N,df(mx)
	  common/turb/ sml(mx), fu(mx), yp(mx), Ap, us, tau
	  common/coef/ A(mx), C(mx), eps(mx), dT
	  common/error/u_old(mx),rel_err,iteration
	  real discha,ave,reynolds,fr_m,fr_d,pi

	  discha = 0.
	  pi     = 22./7.
	  do i=2,N
	  	discha = discha + 2.*pi*rc(i)*((u(i)+u(i-1))/2.)*dy(i)
	  enddo
	  ave = discha/(pi*Rad**2)
	  Reynolds  = 2.*ave*Rad/vis
	  fr_m= 0.25/((log10(5.74/(Reynolds**0.9)))**2)
	  fr_d= (8.*tau)/(rho*ave**2)

	  return
	  end