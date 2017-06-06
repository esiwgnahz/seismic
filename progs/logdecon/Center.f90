! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Center

  use sep
  implicit none

  integer kk,n1,n2,n3,cent
  real, allocatable, dimension(:)  :: in, ou
  

  call sep_init()
  call from_par('n1',n1)
  call from_par('n2',n2)
  call from_par('n3',n3)

  allocate( in(n1*n2*n3), ou(n1*n2*n3) )
  
  ! Read data on Helix
  call sep_read(in)

  do kk=1,n3;
     cent = (n2/2-1)*n1+n1/2-1
     if (kk==1) then
        ou(              1+cent :               n1*n2 ) = in(              1      :              n1*n2-cent )
     else
        ou((kk-1)*n1*n2+ 1      : (kk-1)*n1*n2+ n1*n2 ) = in( (kk-1)*n1*n2+1-cent : (kk-1)*n1*n2+n1*n2-cent )
     end if
  end do
  
  call sep_write(ou)

end program Center
