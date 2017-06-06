! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module FD_types

  use GeneralParam_types

  implicit none 

  type UnscaledFDcoefs
     real :: c0
     real :: c1
     real :: c2
     real :: c3
     real :: c4
  end type UnscaledFDcoefs

  type ScaledFDcoefs
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
  end type ScaledFDcoefs

  type FDbounds
     integer :: nmin1
     integer :: nmax1
     integer :: nmin2
     integer :: nmax2
     integer :: nmin3
     integer :: nmax3
     
  end type FDbounds

contains

  subroutine FD_types_assign_scaled_coefs(coefin,coefou,genpar)
    type(ScaledFDcoefs)    ::                    coefou
    type(UnscaledFDcoefs)  ::             coefin
    type(GeneralParam)     ::                           genpar
    real                   :: dxi,dyi,dzi

    ! genpar%coefpower=2 for scalar   wave equation (v only)
    ! genpar%coefpower=1 for acoustic wave equation (v and rho)
    dxi=1./genpar%delta(2)
    dyi=1./genpar%delta(3)
    dzi=1./genpar%delta(1)
    
    coefou%c0x=coefin%c0*dxi**genpar%coefpower
    coefou%c1x=coefin%c1*dxi**genpar%coefpower
    coefou%c2x=coefin%c2*dxi**genpar%coefpower
    coefou%c3x=coefin%c3*dxi**genpar%coefpower
    coefou%c4x=coefin%c4*dxi**genpar%coefpower

    coefou%c0y=coefin%c0*dyi**genpar%coefpower
    coefou%c1y=coefin%c1*dyi**genpar%coefpower
    coefou%c2y=coefin%c2*dyi**genpar%coefpower
    coefou%c3y=coefin%c3*dyi**genpar%coefpower
    coefou%c4y=coefin%c4*dyi**genpar%coefpower

    coefou%c0z=coefin%c0*dzi**genpar%coefpower
    coefou%c1z=coefin%c1*dzi**genpar%coefpower
    coefou%c2z=coefin%c2*dzi**genpar%coefpower
    coefou%c3z=coefin%c3*dzi**genpar%coefpower
    coefou%c4z=coefin%c4*dzi**genpar%coefpower

  end subroutine FD_types_assign_scaled_coefs

end module FD_types
