program Change

  use sep

  implicit none

  real,    dimension(:,:), allocatable:: in,ou

  integer, dimension(:), allocatable:: n
  real,    dimension(:), allocatable:: o
  real,    dimension(:), allocatable:: d
  character(len=1024)               :: label

  integer :: ndim,i,j,index
  real    :: value
  logical :: found

  call sep_init(SOURCE)
  ndim=sep_dimension()
  allocate(n(ndim), o(ndim), d(ndim))
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do
  
  allocate(in(n(1),product(n(2:))),ou(n(1),product(n(2:))))
  call from_param('value',value,0.)
  call sep_read(in)

  where(in.le.value) 
     ou=0.
  elsewhere
     ou=1.
  end where

  call sep_write(ou)
        
end program Change
