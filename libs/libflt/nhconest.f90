! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module nhconest_mod
  ! Nonstationary convolution, adjoint is the filter.
  use adj_mod
  use nhelix

  use omp_lib

  implicit none

  real ,private,  dimension(:), pointer          :: x 
  type( nfilter), private :: aa 
  integer, private :: nhmax 
  integer,   private                             :: np 

  !$OMP THREADPRIVATE(x,aa)

contains
  subroutine nhconest_mod_init ( x_in,aa_in,nhmax_in )
    real ,   dimension(:), target                    :: x_in 
    type( nfilter)    :: aa_in 
    integer    :: nhmax_in 
    x => x_in 
    aa = aa_in 
    nhmax = nhmax_in 
    np = size( aa%hlx)
  end subroutine nhconest_mod_init
  function nhconest_mod_lop ( adj, add, a, y) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: a,y 
    call adjnull (adj,add,a,y )
    call nhconest_mod_lop2(adj,a,y )
    stat=0
  end function nhconest_mod_lop
  subroutine nhconest_mod_lop2(adj,a,y)
    logical,intent(in) :: adj
    real, dimension ( nhmax, np)  :: a 
    real, dimension (:)  :: y 
    integer                        :: ia, ix, iy, ip
    integer, dimension(:), pointer :: lag

    do iy = 1, size( y)
       if ( aa%mis( iy)) then
          cycle
       end if
       ip = aa%pch( iy)
       lag => aa%hlx( ip)%lag
       do ia = 1, size( lag)
          ix = iy - lag( ia)
          if ( ix < 1) then
             cycle
          end if
          if ( adj) then
             a( ia, ip) =    a( ia, ip) +  y( iy) * x( ix)
          else
             y( iy) =         y( iy) +  a( ia, ip) * x( ix)
          end if
       end do
    end do
  end subroutine nhconest_mod_lop2
  subroutine nhconest_mod_close()
    if (associated(x)) nullify(x)
  end subroutine nhconest_mod_close
end module nhconest_mod
