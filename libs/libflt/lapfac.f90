! Copyright 1999 Stanford University
!All rights reserved
!
!Author: Sergey Fomel, Jon Claerbout, and Paul Sava
module lapfac_mod  
! Factor 2-D Laplacian.
  use wilson_mod
  implicit none
  contains
  function lapfac2_mod( eps, n1, na)   result (aa)
    type( filter)                       :: aa, lap
    real,                   intent( in) :: eps
    integer,                intent( in) :: n1, na
    integer                             :: i
    real                                :: a0, lap0
    call allocatehelix( lap, 2)         ! laplacian filter
    lap0 = 4. + eps                     ! zero lag coeff.
    lap%lag = (/ 1, n1 /)		! lag(1)= 1; lag(2)=n1  # one side only
    lap%flt = -1.			! flt(1)=-1; flt(2)=-1
    call allocatehelix( aa, 2*na)       ! laplacian derivative
    aa%flt = 0.
! probably done already in allocation.
    do i = 1, na  
      aa%lag( i   ) =      i		! early lags (first row)
      aa%lag( na+i) = n1 + i - na	! late lags (second row)
    end do 
    call wilson_mod_init( 10 * n1 )
    call wilson_mod_factor( 20, lap0, lap, a0, aa)
    call wilson_mod_close()
    call deallocatehelix( lap)
  end function 
end module 
