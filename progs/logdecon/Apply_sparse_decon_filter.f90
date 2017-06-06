! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Apply_sparse_decon_filter

  use sep
  use sfft_mod
 

  implicit none

  real,       dimension(:),     allocatable :: in,ou   
  real,       dimension(:),     allocatable :: rt

  complex,    dimension(:),     allocatable :: uw, tmp, inw

  double complex, dimension(:), allocatable :: euw

  real,       dimension(:),     allocatable :: ut


  integer,    dimension(:),     allocatable :: n
  real,       dimension(:),     allocatable :: o, d
  integer(8), dimension(:),     allocatable :: plan1,plan2
  logical,    parameter                     :: FORWARD=.FALSE.,BACKWARD=.TRUE.

  character(len=1024) :: label

  integer :: i,j,k
  integer :: ndim,nt,n2,nshots,nw,nb
  real    :: ow,dw,dt,ot

  logical :: div

  call sep_init(SOURCE)

  ndim=sep_dimension()
  allocate(n(ndim),o(ndim),d(ndim))

  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do

  call from_param('div',div,.false.)

  nt=n(1); dt=d(1); ot=o(1)
  nw=nt/2+1

  if (ndim.ge.3) then
     nshots=product(n(3:))
     n2=n(2)
  else if (ndim.eq.2) then
     nshots=1 
     n2=n(2)
  else if (ndim.eq.1) then
     n2=1
     nshots=1
  end if

  allocate(in(nt),ou(nt),inw(nw),tmp(nw),ut(nt),uw(nw),euw(nw))

  allocate(plan1(2))
  write(0,*) 'INFO: Initialize FFT plans'
  call init_fft_dimension_1d_r2c_c2r(nt,dt,dw,ow)
  call initialize_fft1d_r2c_c2r(plan1)

  if (.not.exist_file('filter')) call erexit('Error: need filter file')

  call from_aux('filter','n1',nb)
  if (nb.ne.nt) call erexit('Error: filter should have same length as input file')
  call sreed('filter',ut,4*nt)

!  call fft1d_r2c_c2r_n23(cshift(ut,-nt/2),uw,1,plan1,FORWARD)
  call fft1d_r2c_c2r_n23(ut,uw,1,plan1,FORWARD)
  if (div) uw=1/uw

  SHOTS:do k=1,nshots
     write(0,*) 'INFO: processing shot',k,'/',nshots


     TRACES: do j=1,n2

        call sreed('in',in,4*nt)
        call fft1d_r2c_c2r_n23(in,inw,1,plan1,FORWARD)
        tmp=inw*uw
        call fft1d_r2c_c2r_n23(ou,tmp,1,plan1,BACKWARD)
        call srite('out',ou,4*nt)        

     end do TRACES

  end do SHOTS
  call sep_close()

end program Apply_sparse_decon_filter
