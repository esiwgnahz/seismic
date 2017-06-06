! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Logistic_mod

  implicit none

contains

  function logistic_comp(rd,k) result (stat)
    real, dimension(:) ::rd
    real               ::   k
    double precision, dimension(:), allocatable :: rdd
    integer            ::stat
    double precision   :: one

    allocate(rdd(size(rd)))

    one=1.0D+0
    rdd=dble(rd)

    rdd=one/(one+exp(-dble(k)*rdd))

    rd=sngl(rdd)
    deallocate(rdd)
    stat=0
  end function logistic_comp

  function Heaviside_comp(rd) result (stat)
    real, dimension(:) :: rd
    integer            ::             stat

    where(rd.le.0)
       rd=0.
    elsewhere
       rd=1
    end where
    
    stat=0
  end function Heaviside_comp

  function logistic_der(rd,k) result (stat)
    real, dimension(:) ::rd
    real               ::   k
    integer            :: stat, stat1

    stat1=logistic_comp(rd,k)
    rd=k*rd*(1-rd)
    stat=0
  end function logistic_der

end module logistic_mod
