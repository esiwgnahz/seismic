! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program IlluScale

  use sep

  implicit none

  real,    dimension(:), allocatable:: in,ou

  integer, dimension(:), allocatable:: n
  real,    dimension(:), allocatable:: o
  real,    dimension(:), allocatable:: d
  character(len=1024)               :: label

  integer :: ndim,i,j,index
  real    :: denom
  logical :: found

  call sep_init(SOURCE)
  ndim=sep_dimension()
  if (ndim.eq.2) ndim=ndim+1
  allocate(n(ndim), o(ndim), d(ndim))
  n=0.
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do
  n(3)=max(n(3),1)
  allocate(in(product(n)),ou(product(n)))

  call from_param('denom',denom,100.)

  write(0,*) 'INFO: Scaling and invert illumination file - (I+max(I)/denom)/sqrt(sum(I''I)', denom
  write(0,*) 'INFO: n=',n

  call sep_read(in) 
  
  ou=(in+maxval(in)/denom)/sqrt(sum(dprod(in,in))/size(in))
  ou=1./ou

  call sep_write(ou)
        
  do i=1,ndim
     call sep_put_data_axis_par("ou",i,n(i),o(i),d(i),label)
  end do

end program IlluScale
