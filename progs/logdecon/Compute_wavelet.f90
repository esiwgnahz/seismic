program Compute_wavelet

  use sep  
  use sfft_mod

  implicit none
 
  complex,        dimension(:), allocatable :: uw, tmp
  double complex, dimension(:), allocatable :: euw
  real,           dimension(:), allocatable :: ut

  integer :: ndim,nt,n3,i,nw,j
  real    :: ot,dt,ow,dw,o3,d3
  logical :: verb

  character(len=1024) :: label

  integer(8), dimension(:),     allocatable :: plan1,plan2
  logical,    parameter                     :: FORWARD=.FALSE.,BACKWARD=.TRUE.

  call sep_init(SOURCE)

  ndim=sep_dimension()

  call sep_get_data_axis_par("in",1,nw,ow,dw,label)
  call sep_get_data_axis_par("in",3,n3,o3,d3,label)
  call from_param('verb',verb,.false.)

  write(0,*) 'INFO: There are ',n3,' filters'

  call from_history('nt',nt)
  call from_history('ot',ot)
  call from_history('dt',dt)
  allocate(plan1(2))
  write(0,*) 'INFO: Initialize FFT plans'
  call init_fft_dimension_1d_r2c_c2r(nt,dt,dw,ow)
  call initialize_fft1d_r2c_c2r(plan1)

  allocate(uw(nw),euw(nw),tmp(nw),ut(nt))

  do i=1,n3

     if (verb) write(0,*) 'INFO:    processing filter',i,'/',n3
     call sreed('in',uw,8*nw)
     euw=cexp(uw)
     tmp=1./euw
          
     do j=2,nw,2
        tmp(j)=-tmp(j)
     end do
        
     call fft1d_r2c_c2r_n23(ut,tmp,1,plan1,BACKWARD)
     
     call srite('out',ut,4*nt)

  end do

  write(0,*) 'INFO: done with processing'
  call to_history('o1',-dt*(nt-1)/2)
  call to_history('d1',dt)
  call to_history('n1',nt)
  call to_history('esize',4)
  call sep_close()
  
end program Compute_wavelet

