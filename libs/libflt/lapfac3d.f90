! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module lapfac3_mod

  use wilson_mod
  use helix

  implicit none

contains
  
  function lapfac3d_mod( eps, n1, n2, na) result (aa)

    type(filter)        :: aa, lap
    real,    intent(in) :: eps
    integer, intent(in) :: n1, n2, na
    integer             :: i,j
    real                :: a0, lap0
    
    call allocatehelix(lap, 3)    ! laplacian filter
    lap0 = 6. + eps               ! zero lag coeff.
    lap%lag = (/ 1, n1, n1*n2 /)  ! lag(1)=1;lag(2)=n1;lag(3)=n1*n2
    lap%flt = -1                  
    
    call allocatehelix(aa, 3*na)
    aa%flt = 0.

    do i=1,na 
       aa%lag(i) = i
       aa%lag(na+i) = n1 + i - na
       aa%lag(2*na+i) = n1*n2 + i - na
    end do

    call wilson_mod_init(10*n1*n2)
    call wilson_mod_factor(30, lap0, lap, a0, aa)
    call wilson_mod_close()
    call deallocatehelix(lap)

  end function lapfac3d_mod

end module lapfac3_mod

