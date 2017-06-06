! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Convolution

  use tcai1
  use sep

  real, dimension(:),      pointer :: bb
  real, dimension(:), allocatable  :: yy,xx

  integer,    dimension(:),     allocatable :: n
  real,       dimension(:),     allocatable :: o, d

  character(len=1024) :: label
  integer             :: ndim,nb,i,j,stat,nt,ny

  call sep_init(SOURCE)

  ndim=sep_dimension()
  allocate(n(ndim),o(ndim),d(ndim))

  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do

  nt=n(1); dt=d(1)
  if (ndim.ge.2) then     
     n23=product(n(2:))
  else
     n23=1
  end if

  allocate(xx(nt))
  if (.not.exist_file('filter')) call erexit('Error: need filter file')

  call from_aux('filter','n1',nb)
  ny=nt+nb-1
  allocate(bb(nb),yy(ny))
  call sreed('filter',bb,4*nb)

  call tcai1_init(bb)
  do i=1,n23
     xx=0.
     yy=0.
     call sreed('in',xx,4*nt)
     stat=tcai1_lop(.false.,.false.,xx,yy)
     do j=nt+1,ny ! Folds the end of the convolution back to origin 
        yy(j-nt)=yy(j-nt)+yy(j)
     end do
     call srite('out',yy,4*nt)
  end do

  call sep_close()

end program Convolution
