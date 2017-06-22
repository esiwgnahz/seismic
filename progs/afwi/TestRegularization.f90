! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program TestRegularization

  use sep
  use helix
  use print
  use helicon_mod
  use helicon2_mod
  use testreg_mod
  use createhelixmod

  real, dimension(:), allocatable :: d,m,w ! Here, d and m are the same
  integer          :: n1, n2, niter ,neval
  real             :: epsx,epszp,epszm,threshp,threshm
  type(filter)     :: xx,zz
  integer, dimension(:), allocatable :: n0,x,y,z,center,gap,npef
  integer :: stat,regtypep,regtypem,nthreads,task,counter

  call sep_init()

  call from_history('n1',n1)
  call from_history('n2',n2)

!  allocate(d(n1*n2),m(n1*n2),g(n1*n2))
  allocate(d(n1*n2),m(n1*n2))
  call sreed('in',d,4*n1*n2)
  allocate(w(n1*n2))
  if(exist_file('weight')) then
     call sreed('weight',w,4*n1*n2)
  else
     w=1.
  end if

  m=0.
  g=0.

  allocate(x(2),y(2),z(2),center(2),gap(2),n0(2),npef(2))
  center=1; gap=0; npef=(/n1,n2/); n0=npef
  z=(/2,1/)
  x=(/1,2/) 
  
  zz=createhelix(npef,center,gap,z)
  xx=createhelix(npef,center,gap,x)            
  zz%flt=-1
  xx%flt=-1

  write(0,*) 'INFO: Filter z:'
  call printn(n0,center,z,zz)
  write(0,*) 'INFO: Filter x:' 
  call printn(n0,center,x,xx)

  call from_param('epszp',epszp,1.0)
  call from_param('epszm',epszm,epszp)
  call from_param('task',task,0)
  call from_param('epsx',epsx,epszp)
  call from_param('regtypep',regtypep,3)
  call from_param('regtypem',regtypem,2)
  call from_param('threshp',threshp,1.0)
  call from_param('threshm',threshm,threshp)
  call from_param('niter',niter,10)
  call from_param('neval',neval,10)
  call from_param('num_threads',nthreads,4)

  write(0,*) 'INFO: threshp=',threshp,'threshm=',threshm
  write(0,*) 'INFO: regtypep=',regtypep,'regtypem=',regtypem
  write(0,*) 'INFO: n1=',n1,'n2=',n2
  write(0,*) 'INFO: nthreads=',nthreads
  write(0,*) 'INFO: epszm=',epszm
  write(0,*) 'INFO: epszp=',epszp,' niter=',niter
  write(0,*) 'INFO: epsx=',epsx,' neval=',neval
  write(0,*) 'INFO: task=',task

  call omp_set_num_threads(nthreads)

  call testreg_init(regtypep,regtypem,epszp,epszm,epsx,threshp,threshm,w)
  call helicon_mod_init(zz)
  call helicon2_mod_init(xx)

  if (task.eq.0) then
     call l_bfgs_smp(simple_test,d,m,n1*n2,niter,neval)
  else
     call l_bfgs_smp(simple_test_hinge,d,m,n1*n2,niter,neval)
  end if

  call testreg_count(counter)
!  stat=simple_test(n1*n2,m,g,f,eps2)

  call srite('out',m,4*n1*n2)

  call to_history('n1',n1,'grad')
  call to_history('n2',n2,'grad')
  call to_history('n3',counter,'grad')

end program TestRegularization

subroutine l_bfgs_smp(fctgdt,d,x,sizex,niter,neval)

  interface
     function fctgdt(sizex,d,x,g,f) result (stat)
       integer                ::  stat
       integer                :: sizex
       real, dimension(sizex)     ::  d
       real, dimension(sizex)     ::  x
       real, dimension(sizex)     ::  g
       double precision       ::  f
    end function fctgdt
  end interface

  integer :: sizex
  real, dimension(sizex)          :: x,d
  integer                     :: count

  real, dimension(:),     allocatable :: g


  double precision, dimension(:), allocatable :: gd
  double precision, dimension(:), allocatable :: xd
  double precision, dimension(:), allocatable :: xsave
  double precision, dimension(:), allocatable :: xdsave
  double precision, dimension(:), allocatable :: wd
  double precision, dimension(:), allocatable :: diagd

  real             :: f,g0
  double precision :: fd,f0

  integer,dimension(2) :: iprint
  double precision     :: EPS,XTOL

  logical              :: cont,found
  integer              :: i,iter,info,niter,neval
  integer              :: stat2,iflag,myinfo

  DOUBLE PRECISION GTOL,STPMIN,STPMAX,XMIN,XMAX
  INTEGER MP,LP,NDIM,NWORK,MSAVE

  external :: LB2S
  COMMON /LB3S/MP,LP,GTOL,STPMIN,STPMAX

  NDIM=size(x)
  write(0,*) 'NDIM',NDIM
  MSAVE=5
  NWORK=NDIM*(2*MSAVE+1)+2*MSAVE

  iprint(1)=1
  iprint(2)=0
  XTOL=epsilon(XTOL)
  EPS=1.0D-20

  allocate(wd(NWORK))

  allocate(diagd(NDIM))
  allocate(gd(NDIM))
  allocate( g(NDIM))
  allocate(xd(NDIM))
  allocate(xsave(NDIM))
  allocate(xdsave(NDIM))

  fd=0.
  g=0.
  gd=0.

  xd=0.
  xsave=0.
  iter=0
  eval=0
  iflag=0
  myinfo=0
  info=0
  count=0
  cont=.True.
  found=.False.

  ! We keep the starting guess in memory
  xsave = dble(x)
  xdsave=xsave

  do while (cont.and.(.not.found))

     stat2=fctgdt(sizex=sizex,d=d,x=x,g=g,f=fd)

     xd=dble(x)
     gd=dble(g)
     xdsave=xd

     call LBFGSS(NDIM,MSAVE,xd,fd,gd,&
     &          .False.,diagd,iprint,EPS,&
     &          XTOL,wd,iflag,myinfo)

     x=sngl(xd)
     info=myinfo

     if (myinfo.eq.1) then
        xsave=xdsave
        myinfo=0
        count=count+1
        if(count.ge.niter) then
           found=.true.
           write(0,*) 'INFO: iter reached max_iter - Iterations will stop'
        end if
     end if

     if(iflag.le.0) then
        cont=.false.
        write(0,*) 'INFO: LINE SEARCH STOPPED IFLAG = ',iflag
        if (iflag.eq.-1) then
           if (info.eq.0) write(0,*) 'INFO:       IMPROPER INPUT PARAMETERS'
           if (info.eq.2) write(0,*) 'INFO:       INTERVAL OF UNCERTAINTY ISSUE'
           if (info.eq.3) write(0,*) 'INFO:       MORE THAN 20 FCT EVALS'
           if (info.eq.4) write(0,*) 'INFO:       STEP TOO SMALL'
           if (info.eq.4) write(0,*) 'INFO:       STEP TOO LARGE'
           if (info.eq.6) write(0,*) 'INFO:       ROUNDING ERRORS'
        end if
     end if

     eval=eval+1
     if(eval.ge.neval) then
        write(0,*) 'INFO: iter reached neval - Iterations will stop'
        cont=.false.
     end if
  end do

  ! We didn't really find a new model vector, the line-search
  ! was still going. We take the previous best solution.
  if (.not.cont) then
     x=sngl(xsave)
  end if
  deallocate(xd,xsave,xdsave,g,wd,gd,diagd)

end subroutine l_bfgs_smp
