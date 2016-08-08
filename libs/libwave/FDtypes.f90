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
    type(ScaledFDcoefs):: coefou
    type(UnscaledFDcoefs)  :: coefin
    type(GeneralParam) :: genpar
    real               :: dxi,dyi,dzi

    dxi=1./genpar%dx
    dyi=1./genpar%dy
    dzi=1./genpar%dz
    
    coefou%c0x=coefin%c0*dxi*dxi
    coefou%c1x=coefin%c1*dxi*dxi
    coefou%c2x=coefin%c2*dxi*dxi
    coefou%c3x=coefin%c3*dxi*dxi
    coefou%c4x=coefin%c4*dxi*dxi

    coefou%c0y=coefin%c0*dyi*dyi
    coefou%c1y=coefin%c1*dyi*dyi
    coefou%c2y=coefin%c2*dyi*dyi
    coefou%c3y=coefin%c3*dyi*dyi
    coefou%c4y=coefin%c4*dyi*dyi

    coefou%c0z=coefin%c0*dzi*dzi
    coefou%c1z=coefin%c1*dzi*dzi
    coefou%c2z=coefin%c2*dzi*dzi
    coefou%c3z=coefin%c3*dzi*dzi
    coefou%c4z=coefin%c4*dzi*dzi

  end subroutine FD_types_assign_scaled_coefs

end module FD_types
