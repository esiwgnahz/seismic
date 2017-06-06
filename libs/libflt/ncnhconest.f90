! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module ncnhconest  
! non-causal non-stationary convolution, adjoint is the filter.
  use adj_mod
  use nhelix
  use omp_lib
  implicit none
  real ,private,  dimension(:), pointer                  :: x 
  type( nfilter), private :: aa 
  integer, private :: nhmax 
  integer,   private                             :: np

  contains
  subroutine ncnhconest_init ( x_in,aa_in,nhmax_in )
    real ,   dimension(:),target                :: x_in 
    type( nfilter)    :: aa_in 
    integer    :: nhmax_in 
    x => x_in 
    aa = aa_in 
    nhmax = nhmax_in 
    np = size( aa%hlx)
  end subroutine
  function ncnhconest_lop ( adj, add, a, y) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: a,y 
    call adjnull (adj,add,a,y )
    call ncnhconest_lop3(adj,a,y )
!    call ncnhconest_lop2(adj,a,y )
    stat=0
  end function 
  subroutine ncnhconest_lop2(adj,a,y)
    logical,intent(in) :: adj 
    real, dimension ( nhmax, np)  :: a 
    real, dimension (:)  :: y 
    integer                        :: ia, ix, iy, ip, nx
    integer, dimension(:), pointer :: lag
    nx=size(x)

    !$OMP PARALLEL DO PRIVATE(iy,ip,lag,ia,ix)
    do iy = 1, size( y)
       if (associated(aa%mis)) then
          if ( aa%mis( iy)) then
             cycle
          end if
       end if
       ip = aa%pch( iy)
       lag => aa%hlx( ip)%lag

       do ia = 1, size( lag)
          ix = iy - lag( ia)
          if ((ix<1).or.(ix>nx)) then
             cycle
          end if
          if ( adj) then
             a( ia, ip) =    a( ia, ip) +  y( iy) * x( ix)
          else
             y( iy) =         y( iy) +  a( ia, ip) * x( ix)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine ncnhconest_lop2
  subroutine ncnhconest_lop3(adj,a,y)
    logical,intent(in) :: adj 
    real, dimension ( nhmax, np)  :: a 
    real, dimension (:)  :: y 
    integer                        :: ia, ix, iy, ip, nx
    integer, dimension(:), pointer :: lag
    nx=size(x)

    !$OMP PARALLEL DO PRIVATE(iy,ip,lag,ia,ix)
    do iy = 1, size( y)
       if (associated(aa%mis)) then
          if ( aa%mis( iy)) then
             cycle
          end if
       end if
       ip = aa%pch( iy)
       lag => aa%hlx( ip)%lag

       do ia = 1, size( lag)
          ix = iy - lag( ia)
          if ((ix<1).or.(ix>nx)) then
             cycle
          end if
          if ( adj) then
             a( ia, ip) =    a( ia, ip) +  y( ix) * x( iy)
          else
             y( ix) =         y( ix) +  a( ia, ip) * x( iy)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine ncnhconest_lop3
  subroutine ncnhconest_close()
    if (associated(x)) nullify(x)
  end subroutine ncnhconest_close
end module ncnhconest
