module FD_types

  implicit none 

  type FDcoefs
     real :: c0
     real :: c1
     real :: c2
     real :: c3
     real :: c4
  end type FDcoefs

  type ScaledFDcoef
     real :: c0x
     real :: c1x
     real :: c2x
     real :: c3x
     real :: c4x
     real :: c0y
     real :: c1y
     real :: c2y
     real :: c3y
     real :: c4y
     real :: c0z
     real :: c1z
     real :: c2z
     real :: c3z
     real :: c4z
  end type ScaledFDcoef

  type FDbounds
     real :: min1
     real :: max1
     real :: min2
     real :: max2
     real :: min3
     real :: max3
  end type FDbounds

end module FD_types
