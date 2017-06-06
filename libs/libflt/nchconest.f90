! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module nchconest  
  use omp_lib
! non-causal masked helix convolution, adjoint is the filter.
  use adj_mod
  use helix
  implicit none
  real ,private,  dimension (:), pointer  :: x 
  type( filter), private :: aa

  !$OMP THREADPRIVATE(x,aa)
  contains
  subroutine nchconest_init ( x_in,aa_in )
    real ,   dimension (:), target :: x_in 
    type( filter)    :: aa_in 
    x => x_in 
    aa = aa_in
  end subroutine
  function nchconest_lop ( adj, add, a, y) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: a,y 
    call adjnull (adj,add,a,y )
    call nchconest_lop2(adj,a,y )
    stat=0
  end function 
  subroutine nchconest_lop2(adj,a,y)
    logical,intent(in) :: adj 
    real, dimension (:)  :: a 
    real, dimension (:)  :: y 
    integer  ia, ix, iy
    do ia = 1, size( a)
      do iy = 1+max(aa%lag(ia),0), size(y)+min(aa%lag(ia),0)
        if ( associated(aa%mis)) then
           if ( aa%mis( iy)) then
              cycle  
           end if
        end if
        ix = iy - aa%lag( ia)
        if ( adj) then
          a( ia) =    a( ia) +  y( iy) * x( ix)
        else
          y( iy) =    y( iy) +  a( ia) * x( ix)
        end if
      end do
    end do
  end subroutine 
  subroutine nchconest_close()
    if (associated(x)) nullify(x)
  end subroutine 
end module 
