! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Apply_filter

  use sep
  use adj_mod
  use sfft_mod
  use dottest
  use symmetry_mod
  use quantile_calc_mod

  implicit none

  
  real, dimension(2) :: dot1, dot2

  real,       dimension(:,:),   allocatable :: in,gin  ! input data
  complex,    dimension(:,:),   allocatable :: inw,ginw! input Fourier domain
  real,       dimension(:,:),   allocatable :: gt      ! gain function
  real,       dimension(:,:),   allocatable :: qt
  real,       dimension(:,:),   allocatable :: rt,rtmp
  real,       dimension(:),     allocatable :: regt
  real,       dimension(:),     allocatable :: reg2t
  real,       dimension(:),     allocatable :: reg3t
  real,       dimension(:),     allocatable :: reg4t
  
  real,             dimension(:), allocatable :: ut,utp ! filter time domain

  complex,        dimension(:), allocatable :: uw
  complex,        dimension(:), allocatable :: tmp,tmp2

  integer,    dimension(:),     allocatable :: n
  real,       dimension(:),     allocatable :: o, d
  integer(8), dimension(:),     allocatable :: plan1,plan2
  logical,    parameter                     :: FORWARD=.FALSE.,BACKWARD=.TRUE.

  integer :: ntf,nwf
  real :: owf,dwf

  character(len=1024) :: label
  integer :: i,j,it,n_dim,nt,nt2,n23,nw,nq,max_iter,neval,count,lagp,lagn,verb,stat
  real    :: ow,dw,dt,ot,oky,dky,percentile,t,t2,maxin,quant,weight,eps2,qbar,eps3,max1
  logical :: zero_neg,new_init,phase_only,overwrite_quant
  double complex, dimension(:), allocatable :: euw

  !LBFGS arrays
  double precision, dimension(:), allocatable :: wd
  double precision, dimension(:), allocatable :: diagd

  real             :: f
  double precision :: fd

  integer,dimension(2) :: iprint
  double precision     :: EPS,XTOL
  
  logical              :: cont,inv
  integer              :: iter,info,nmax1,norm
  integer              :: stat2,iflag,myinfo,lag,maxlag
  
  DOUBLE PRECISION GTOL,STPMIN,STPMAX
  INTEGER MP,LP,NDIM,NWORK,MSAVE

  external :: LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

  call sep_init(SOURCE)
  n_dim=sep_dimension()
  allocate(n(n_dim),o(n_dim),d(n_dim))

  do i=1,n_dim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do

  call sep_get_data_axis_par("filter",1,nwf,owf,dwf,label)
  call from_aux("filter",'nt',ntf)
  call from_param('inv',inv,.false.)
  call from_param('lag',lag,ntf/20)
  call from_param('maxlag',maxlag,0)
  
  nt=n(1); dt=d(1); ot=o(1)

  if (n_dim.ge.2) then
     n23=product(n(2:))
  else
     n23=1
  end if
  nw=nt/2+1
  nt2=nt/2
 
  allocate(in(nt,n23))
  
  write(0,*) 'FFT input data'
  allocate(plan1(2),plan2(2))
  write(0,*) 'INFO: Initialize FFT plans'
  call init_fft_dimension_1d_r2c_c2r(nt,dt,dw,ow)
  call initialize_fft1d_r2c_c2r(plan1)
  allocate(inw(nw,n23))
  call sreed('in',in,4*nt*n23)
  write(0,*) 'INFO: FFT input data'
  call fft1d_r2c_c2r_n23(in,inw,n23,plan1,FORWARD)
  call destroy_fft(plan1)

  allocate(uw(nwf),tmp(nw))
  allocate(ut(ntf))
  allocate(rt(nt,n23))

  uw=0.

  write(0,*) 'INFO: Read filter in frequency and FFT it to time' 
  call sreed('filter',uw,8*nwf)
  call init_fft_dimension_1d_r2c_c2r(ntf,dt,dwf,owf)
  call initialize_fft1d_r2c_c2r(plan1)
  call fft1d_r2c_c2r_n23(ut,uw,1,plan1,BACKWARD)
  call destroy_fft(plan1)
 
  if (maxlag.ne.0) then
     write(0,*) 'INFO: maxlag=',maxlag
     write(0,*) 'INFO:    lag=',lag
     do it=maxlag+1,maxlag+1+lag
        weight=cos(.5*3.14159*(it-maxlag-1)/(lag))**2
        ut(it)=weight*ut(it)
     end do
     ut(maxlag+2+lag:ntf/2)=0.
     call init_fft_dimension_1d_r2c_c2r(ntf,dt,dwf,owf)
     call initialize_fft1d_r2c_c2r(plan1)
     call fft1d_r2c_c2r_n23(ut,uw,1,plan1,FORWARD)
     call srite('flt_out',uw,8*nwf)
     call to_history('n1',nwf,'flt_out')
     call to_history('d1',dwf,'flt_out')
     call to_history('o1',owf,'flt_out')
     call to_history('nt',ntf,'flt_out')
     call to_history('dt',dt,'flt_out')
     call to_history('ot',ot,'flt_out')
     call to_history('esize',8,'flt_out')
     uw=0.
  end if

  call init_fft_dimension_1d_r2c_c2r(nt,dt,dw,ow)
  call initialize_fft1d_r2c_c2r(plan1)
  
  write(0,*) 'INFO: Copy time filter to bigger filter array'
  allocate(utp(nt))
  utp=0.
  utp(1:ntf/2)=ut(1:ntf/2)
  utp(nt-ntf/2+1:nt)=ut(ntf/2+1:)
  deallocate(uw)
  allocate(uw(nw))

  call srite('time_filter',ut,4*ntf)
  call to_history('n1',ntf,'time_filter')
  call srite('time_filter2',utp,4*nt)
  call to_history('n1',nt,'time_filter2')

  write(0,*) 'INFO: FFT time filter'
  call fft1d_r2c_c2r_n23(utp,uw,1,plan1,FORWARD)


  if (inv) then    
     write(0,*) 'INFO: Apply inverse filter to data: add wavelet'
     do j=1,n23
        tmp=inw(:,j)*exp(-uw)
        call fft1d_r2c_c2r_n23(rt(:,j),tmp,1,plan1,BACKWARD) ! r  = IFT(D*exp(U))
     end do
  else
     write(0,*) 'INFO: Apply filter to data: remove wavelet'
     do j=1,n23
        tmp=inw(:,j)*exp(uw)
        call fft1d_r2c_c2r_n23(rt(:,j),tmp,1,plan1,BACKWARD) ! r  = IFT(D*exp(U))
     end do
  end if

  call destroy_fft(plan1)
  call srite('out',rt,4*nt*n23)
   
     
  call sep_close()

end program Apply_filter


subroutine weight_lags(wt,nt,nmax1,lagp,lagn)
  real, dimension(nt)::wt
  integer            ::   nt,nmax1,lagp,lagn

  integer :: it

  wt=0.
  do it=1,lagn
     wt(nt-it+1)=sin(.5*3.14159*(it-1)/(lagn-1))**2
  end do
  do it=1,lagp
     wt(nmax1-it+1)=cos(.5*3.14159*(it-1)/(lagp-1))**2
  end do
  wt(nmax1+1:nt-lagn)=1.

end subroutine weight_lags
