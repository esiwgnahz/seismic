Program Cliptomax
  use sep

  implicit none

  real,   dimension(:),allocatable :: in
  integer,dimension(:),allocatable :: n
  real    :: o,d
  integer :: ndim,i,ntotal
  real    :: max_val
  character(len=1024):: label

  call sep_init(SOURCE)
  ndim=sep_dimension()
  allocate(n(ndim))
 
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o,d,label)
  end do
  
  ntotal=product(n)
  write(0,*) "INFO n=",n," ntotal=",ntotal

  allocate(in(ntotal))
  call sreed('in',in,4*ntotal)
  
  max_val=maxval(in)
  write(0,*) "INFO maxval in=",max_val

  call srite('out',in/max_val,4*ntotal)
  call sep_close()

end Program Cliptomax
