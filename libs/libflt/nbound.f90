module nbound_mod  
! mark helix filter outputs where input is off data.
  use nhelix
  use bound_mod
  implicit none
  contains 
  subroutine nboundn_mod ( ip, nd, na, aa)
    integer, dimension(:), intent( in) :: nd, na ! (ndim)
    type( nfilter)                     :: aa
    integer, intent (in)               :: ip
    allocate( aa%mis( product (nd)))
    call boundn_mod (nd, nd, na, aa%hlx( ip))
    aa%mis = aa%hlx( ip)%mis
    deallocate( aa%hlx( ip)%mis)
    nullify( aa%hlx( ip)%mis)
  end subroutine 
end module 
