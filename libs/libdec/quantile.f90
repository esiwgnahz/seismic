! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module quantile_calc_mod

  implicit none

contains

  recursive function quantile_calc( k, a) result( value)
    integer,             intent (in)  :: k          
    ! position in array
    real, dimension (:), intent (in)  :: a
    real                              :: value      
    ! output value of quantile
    integer                           :: j
    real                              :: ak
    value=0.
    ak = a( k)
    j = count( a < ak) 
    ! how many a(:) < ak
    if ( j >= k) then
       value = quantile_calc( k, pack( a, a < ak))
    else
       j = count( a > ak) + k - size( a)
       if ( j > 0) then
          value = quantile_calc( j, pack( a, a > ak))
       else
          value = ak
       end if
    end if
  end function quantile_calc

end module quantile_calc_mod
