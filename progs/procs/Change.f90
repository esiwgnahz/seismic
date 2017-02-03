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
  if (ndim.eq.2) ndim=ndim+1
  allocate(n(ndim), o(ndim), d(ndim))
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do
  n(3)=1
  allocate(in(n(1),n(2)),ou(n(1),n(2)))
  call from_param('value',value,0.)

  write(0,*) 'INFO: value below which everything is set to zero= ',value
  write(0,*) 'INFO: n=',n
  do j=1,n(3)

     in=0.
     ou=0.

     call sep_read(in)
     
     where(in.le.value) 
        ou=0.
     elsewhere
        ou=1.
     end where
     
     call sep_write(ou)
        
  end do

  do i=1,ndim
     call sep_put_data_axis_par("ou",i,n(i),o(i),d(i),label)
  end do

end program Change
