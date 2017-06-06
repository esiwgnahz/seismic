! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module ncnbound  
  ! mark helix filter outputs where input is off data.
  use nhelix
  use ncbound
  use Filter_types
  implicit none
contains 
  subroutine ncnboundn ( ip, nd, na, aa)
    integer, dimension(:), intent( in) :: nd, na ! (ndim)
    type(NSfilter)                     :: aa
    integer, intent (in)               :: ip
    allocate( aa%nmatch%mis( product (nd)))
    call ncboundn (nd, nd, na, aa%nmatch%hlx( ip))
    aa%nmatch%mis = aa%nmatch%hlx( ip)%mis
    deallocate( aa%nmatch%hlx( ip)%mis)
    nullify( aa%nmatch%hlx( ip)%mis)
  end subroutine ncnboundn
end module ncnbound
