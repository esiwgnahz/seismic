! Copyright 1999 Stanford University
!All rights reserved
!
!Author: Sergey Fomel
!
module createnhelixmod_smp 
  ! Create non-stationary helix filter lags and mis

  use createhelixmod
  use nhelix
  use helix
  use nbound_mod

  implicit none

contains

  function createnhelix_smp( nd, center, gap, na, pch, ncoef, nmis) result (nsaa)
    type( nfilter)                     :: nsaa   ! needed by nhelicon	
    type( filter)                      :: aa

    integer, dimension(:)  :: nd, na ! data and filter axes
    integer, dimension(:)  :: center ! normally (na1/2,na2/2,...,1)
    integer, dimension(:)  :: gap    ! normally ( 0,   0,  0,...,0)
    integer, dimension(:)  :: pch	 ! (prod(nd)) patching
    logical, dimension(:), optional    :: nmis	

    integer                            :: n123, np, ip, ncoef
    integer, dimension(:), allocatable :: nh

    aa = createhelix( nd, center, gap, na)

    n123 = product( nd)
    np = maxval(pch)

    allocate (nh (np))
    nh = size (aa%lag)

    call nallocate( nsaa, nh, pch)

    ncoef=size(aa%lag)

    do ip = 1, np 
       nsaa%hlx( ip)%lag = aa%lag
    end do

    call deallocatehelix (aa)

    if(present(nmis)) then
       allocate(nsaa%mis(n123))
       nsaa%mis=nmis
    else
       call nboundn_mod(1, nd, na, nsaa)
    end if
   
    deallocate(nh)

  end function createnhelix_smp

end module createnhelixmod_smp
