
program Sparse_Decon_LBFGS

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
  real,       dimension(:,:),   allocatable :: rt,rtmp,rt_save
  real,       dimension(:),     allocatable :: regt
  real,       dimension(:),     allocatable :: reg2t
  real,       dimension(:),     allocatable :: reg3t
  real,       dimension(:),     allocatable :: reg4t
  
  real,             dimension(:), allocatable :: ut, dut, uts ! filter time domain
  double precision, dimension(:), allocatable :: ut_dp,dut_dp,ut_save
  real,             dimension(:), allocatable :: tmp1,utb,regw

  double complex, dimension(:), allocatable :: euw
  complex,        dimension(:), allocatable :: uw, duw
  complex,        dimension(:), allocatable :: tmp,tmp2

  integer,    dimension(:),     allocatable :: n
  real,       dimension(:),     allocatable :: o, d
  integer(8), dimension(:),     allocatable :: plan1,plan2
  logical,    parameter                     :: FORWARD=.FALSE.,BACKWARD=.TRUE.

  character(len=1024) :: label
  integer :: i,j,it,n_dim,nt,nt2,n23,nw,nq,max_iter,neval,count,lagp,lagn,verb,stat
  real    :: ow,dw,dt,ot,oky,dky,percentile,t,t2,maxin,quant,weight,eps2,qbar,eps3,max1
  logical :: zero_neg,new_init,phase_only

  !LBFGS arrays
  double precision, dimension(:), allocatable :: wd
  double precision, dimension(:), allocatable :: diagd

  real             :: f
  double precision :: fd

  integer,dimension(2) :: iprint
  double precision     :: EPS,XTOL
  
  logical              :: cont
  integer              :: iter,info,nmax1,norm
  integer              :: stat2,iflag,myinfo
  
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

  call from_param('t2',t2,0.)

  nt=n(1); dt=d(1); ot=o(1)

  !LBFGS INIT 
  NDIM=nt
  MSAVE=5
  NWORK=NDIM*(2*MSAVE+1)+2*MSAVE

  call from_param('verb',verb,0)

  iprint(1)=verb
  iprint(2)=0
  XTOL=epsilon(XTOL)
  EPS=1.0D-20

  if (n_dim.ge.2) then
     n23=product(n(2:))
  else
     n23=1
  end if
  nw=nt/2+1
  nt2=nt/2

  call from_param('lag',lagn,12)
  lagp=lagn
  call from_param('max1',max1,nt2*d(1)/3)
  call from_param('eps2',eps2,1.)
  call from_param('eps3',eps3,1.)
  call from_param('norm',norm,12)
  call from_param('zero_neg',zero_neg,.false.)
  call from_param('phase_only',phase_only,.false.)
  call from_param('new_init',new_init,.false.)
  call from_param('percentile',percentile,90.)
  if ((percentile.lt.0).or.(percentile.gt.100.)) then
     call erexit('ERROR: percentile must be between 0 and 100')
  end if

  nmax1=max1/d(1)
  if (nmax1.gt.nt2) then
     write(0,*) 'INFO: max1 is too big>',nt2,'change max1 to maxnt'
     nmax1=nt2
  end if
 
  allocate(in(nt,n23),gt(nt,n23))
  gt=1.
  if (exist_file('gain')) then
     call sreed('gain',gt,4*nt*n23)
  end if

  do i=1,nt
     t=(i-1)*dt+ot
     gt(i,:)=t**t2*gt(i,:)
  end do

  if (exist_file('gt')) then
     call srite('gt',gt,4*nt*n23)
     call to_history('n1',nt,'gt')
     call to_history('d1',dt,'gt')
     call to_history('o1',ot,'gt')
     call to_history('n2',n23,'gt')
     call to_history('d2',d(2),'gt')
     call to_history('o2',o(2),'gt')
  end if

  allocate(wd(NWORK))  
  allocate(diagd(NDIM))
  allocate(dut_dp(NDIM))! double precision
  allocate(dut(NDIM))   ! single precision
  allocate(ut_dp(NDIM)) ! double precision
  allocate(ut(NDIM))    ! single precision
  allocate(uts(NDIM))    ! single precision
  
  iter=0
  iflag=0
  cont=.True.
  myinfo=0
  info=0
  count=0 

  call from_param('niter',max_iter,10)
  neval=max_iter

  allocate(plan1(2),plan2(2))
  write(0,*) 'INFO: Initialize FFT plans'
  call init_fft_dimension_1d_r2c_c2r(nt,dt,dw,ow)
  call initialize_fft1d_r2c_c2r(plan1)

  allocate(inw(nw,n23))
  call sreed('in',in,4*nt*n23)

  write(0,*) 'INFO: FFT input data'
  call fft1d_r2c_c2r_n23(in,inw,n23,plan1,FORWARD)

  allocate(uw(nw),duw(nw),euw(nw),tmp(nw),tmp2(nw))
  allocate(rt(nt,n23),rt_save(nt,n23),rtmp(nt,n23),qt(nt,n23),regt(nt),reg2t(nt),reg3t(nt),reg4t(nt),regw(nt))

  ut=0.
  uw=0.
  if (new_init) then
     allocate(ginw(nw,n23),gin(nt,n23))
     gin=gt*in
     call fft1d_r2c_c2r_n23(gin,ginw,n23,plan1,FORWARD)     
     allocate(utb(nt))
     ! Forming the regularization array utb
     tmp=0.
     do j=1,n23
        tmp=tmp+ginw(:,j)*conjg(ginw(:,j))
     end do
     tmp=sqrt(real(tmp))
     uw=clog(tmp)
     uw=uw-sum(uw)/nw
     call fft1d_r2c_c2r_n23(utb,uw,1,plan1,BACKWARD)
     call mostly_causal_amplitude(utb,nt,lagn)
     call fft1d_r2c_c2r_n23(utb,uw,1,plan1,FORWARD)
     ut=-utb
     uw=-uw
     deallocate(gin,ginw,utb)
     
  end if

  if (exist_file('utb')) then
     call srite('utb',uw,8*nw)
     call to_history('n1',nw,'utb')
     call to_history('d1',dw,'utb')
     call to_history('o1',ot,'utb')
     call to_history('n2',1,'utb')
     call to_history('n3',1,'utb')
     call to_history('nt',nt,'utb')
     call to_history('ot',ot,'utb')
     call to_history('dt',dt,'utb')
     call to_history('esize',8,'utb')
  end if

  allocate(tmp1(nt*n23))

  do j=1,n23
     tmp=inw(:,j)*exp(uw)
     call fft1d_r2c_c2r_n23(rt(:,j),tmp,1,plan1,BACKWARD) ! r  = IFT(D*exp(U))
  end do

  do i=1,n23
     do j=1,nt
        tmp1(j+(i-1)*nt)=abs(rt(j,i)*gt(j,i))
     end do
  end do
  nq=size(in)*percentile/100.
  quant=quantile_calc(nq,tmp1); deallocate(tmp1)

  allocate(tmp1(nt))

  if (quant.eq.0.) then
     write(0,*) 'INFO: quant=',quant
     write(0,*) 'INFO: replace quant with maxval'
     quant=maxval(abs(rt*gt))
  end if
  write(0,*) 'INFO:',percentile,'th percentile gained data=',quant

  gt=gt/quant
  qt=gt*rt

  duw=0.
  euw=0.

  call weight_lags(regw,nt,nmax1,lagp,lagn)
  if (exist_file('weight_reg')) call srite('weight_reg',regw,4*nt)

  call symmetry_init(lagn,nt2,nt)

  call dot_test(symmetry_lop,nt,nt,dot1,dot2)

  do while (cont)

     if(info.eq.1) then
        ut_save=ut_dp
        rt_save=rt
        uts=sngl(ut_dp)
     end if

     fd=0.
     duw=0.
     dut=0.

     ! Compute gradient and function
     call fft1d_r2c_c2r_n23(ut,uw,1,plan1,FORWARD)
  
     stat=symmetry_lop(.false.,.false.,ut,regt)   ! r_u2=W_2Ju
     stat=symmetry_lop(.true.,.false.,reg3t,regt) ! Du_2=J'W_2r_u2
     reg2t=regw*ut                                ! r_u1=W_1u

     fd=fd+eps2*sum(dprod(reg2t,reg2t))*0.5+eps3*sum(dprod(regt,regt))*0.5

     euw=exp(uw)
     do j=1,n23
        tmp2=inw(:,j)*euw
        call fft1d_r2c_c2r_n23(rt(:,j),tmp2,1,plan1,BACKWARD)! r  = IFT(D*exp(U))
        qt(:,j)=gt(:,j)*rt(:,j)                              ! q  = g*r
	if (norm.eq.2) then
		fd = fd+sum(qt(:,j)**2)/2
                tmp1=gt(:,j)*qt(:,j)	
	else if (norm.eq.1) then
                fd = fd+sum(abs(qt(:,j)))
                where (qt(:,j).eq.0) 
                  tmp1=0
                else where(qt(:,j).lt.0)
                  tmp1=-gt(:,j)
                else where(qt(:,j).gt.0)
                  tmp1=gt(:,j)
                end where
	else if (norm.eq.12) then
       		fd = fd+sum(sqrt(1.d0+qt(:,j)**2)-1.d0)       ! f  = f+Sum(sqrt(1+q^2)-1)
                tmp1=gt(:,j)*qt(:,j)/sqrt(qt(:,j)**2+1)       !tmp1= g*softclip(q)
	end if
        tmp1=gt(:,j)*qt(:,j)/sqrt(qt(:,j)**2+1)              !tmp1= g*softclip(q)
        call fft1d_r2c_c2r_n23(tmp1,tmp,1,plan1,FORWARD)     !tmp = FT(tmp1)
        call fft1d_r2c_c2r_n23(rt(:,j),tmp2,1,plan1,FORWARD) !tmp2= FT(r) 
        duw=duw+conjg(tmp2)*tmp
        
     end do
     call fft1d_r2c_c2r_n23(dut,duw,1,plan1,BACKWARD)
     
     dut=dut+eps3*reg3t+eps2*regw*reg2t

     if (zero_neg) dut=(1-regw)*dut

     call fft1d_r2c_c2r_n23(dut,duw,1,plan1,FORWARD)
     if (phase_only) duw=cmplx(0.,aimag(duw))
     duw=(duw-sum(duw)/nw)
!     if (phase_only) duw=cmplx(0.,aimag(duw))
     call fft1d_r2c_c2r_n23(dut,duw,1,plan1,BACKWARD)

     dut_dp=dble(dut)
     ut_dp=dble(ut)

     call LBFGSS(NDIM,MSAVE,ut_dp,fd,dut_dp,&
     &          .False.,diagd,iprint,EPS,&
     &          XTOL,wd,iflag,myinfo)

     ut=sngl(ut_dp)

     info=myinfo
     if (myinfo.eq.1) then
          myinfo=0
          count=count+1
          if (count.ge.max_iter) cont=.false.
!          if (exist_file('movie'))        call srite('movie',rt,4*nt*n23)
!          if (exist_file('filter_movie')) then
!             call fft1d_r2c_c2r_n23(ut,uw,1,plan1,FORWARD)
!             call srite('filter_movie',uw,8*nw)
!          end if
!          if (exist_file('fct'))          call srite('fct',sngl(fd),4)
!          if(exist_file('qbar')) then 
!             qbar=0
!             do j=1,n23
!                qbar=qbar+sum(sqrt(1.+qt(:,j)**2))
!             end do
!             qbar=sqrt((qbar/(nt*n23))**2-1)
!             call srite('qbar',qbar,4)
!          end if
!          if (exist_file('qt_movie')) then
!             qbar=0
!             do j=1,n23
!                qbar=qbar+sum(sqrt(1.+qt(:,j)**2))
!             end do
!             qbar=sqrt((qbar/(nt*n23))**2-1)
!             rtmp=qt
!             do j=1,50
!                rtmp(1:50,j)=qbar
!             end do
!             do j=n23-50,n23
!                rtmp(1:50,j)=-qbar
!             end do
!             call srite('qt_movie',rtmp,4*nt*n23)
!          end if
       end if

       if(iflag.le.0) cont=.false.
       iter=iter+1
       if(iter.ge.neval) cont=.false.

    end do

    call srite('out',rt_save,4*nt*n23)
    if (exist_file('qt')) call srite('qt',qt,4*nt*n23)

  if (exist_file('filter')) then
     call fft1d_r2c_c2r_n23(uts,uw,1,plan1,FORWARD)
     call srite('filter',uw,8*nw)
     call to_history('n1',nw,'filter')
     call to_history('d1',dw,'filter')
     call to_history('n2',1,'filter')
     call to_history('n3',1,'filter')
     call to_history('nt',nt,'filter')
     call to_history('ot',ot,'filter')
     call to_history('dt',dt,'filter')
     call to_history('esize',8,'filter')
  end if
  if (exist_file('qbar')) then
     call to_history('n1',count,'qbar')
     call to_history('d1',1,'qbar')
     call to_history('o1',1,'qbar')
  end if
  if (exist_file('filter_movie')) then
     call to_history('n1',nw,'filter_movie')
     call to_history('d1',dw,'filter_movie')
     call to_history('o1',ot,'filter_movie')
     call to_history('n2',1,'filter_movie')
     call to_history('n3',count,'filter_movie')
     call to_history('nt',nt,'filter_movie')
     call to_history('ot',ot,'filter_movie')
     call to_history('dt',dt,'filter_movie')
     call to_history('esize',8,'filter_movie')
  end if
  if (exist_file('qt_movie')) then
     call to_history('n3',count,'qt_movie')
     call to_history('n2',n23,'qt_movie')
     call to_history('n1',n(1),'qt_movie')
     call to_history('d1',d(1),'qt_movie')
     call to_history('d2',d(2),'qt_movie')
     call to_history('o2',o(2),'qt_movie')
  end if
  if (exist_file('qt')) then
     call to_history('n3',1,'qt')
     call to_history('n2',n23,'qt')
     call to_history('n1',n(1),'qt')
     call to_history('o1',o(1),'qt')
     call to_history('d1',d(1),'qt')
     call to_history('d2',d(2),'qt')
     call to_history('o2',o(2),'qt')
  end if
  if (exist_file('movie')) then
     call to_history('n3',count,'movie')
     call to_history('n2',n23,'movie')
     call to_history('n1',n(1),'movie')
     call to_history('d1',d(1),'movie')
     call to_history('d2',d(2),'movie')
     call to_history('o2',o(2),'movie')
  end if
  if (exist_file('fct')) then
     call to_history('n1',count,'fct')
     call to_history('d1',1,'fct')
     call to_history('o1',1,'fct')
  end if
  if (exist_file('weight_reg')) then
     call to_history('n1',nt,'weight_reg')
     call to_history('d1',dt,'weight_reg')
     call to_history('o1',0.,'weight_reg')
  end if

     
  call sep_close()

end program Sparse_Decon_LBFGS

subroutine mostly_causal_amplitude(ut,nt,lag)
  integer            ::               nt,lag
  real, dimension(nt)::            ut

  integer :: it,nt2
  real    :: weight, asym

  nt2=nt/2
  do it=2,nt2
     ut(it)         = 2.*ut(it)
     ut(2*nt2-it+2) = 0.
  end do

  do it=2,lag
     asym=(ut(it)-ut(2*nt2-it+2))/2.
     weight=cos(.5*3.14159*(it-1)/(lag-1.))**2
     ut(      it)   = ut(      it)   - weight*asym
     ut(2*nt2-it+2) = ut(2*nt2-it+2) + weight*asym
  end do
   
end subroutine mostly_causal_amplitude

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
