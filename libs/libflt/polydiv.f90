! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module polydiv_omp  
  ! Helix polynomial division
  use adj_mod
  use helix

  use omp_lib


  implicit none
  
  integer, private :: nd 
  type( filter), private :: aa 
  real ,private,  dimension (:), allocatable  :: tt 

  !$omp threadprivate(aa,tt)

contains
  subroutine polydiv_omp_init ( nd_in,aa_in )
    integer    :: nd_in 
    type( filter)    :: aa_in 
    nd = nd_in 
    aa = aa_in 
    if (.not. allocated(tt)) then
       allocate(tt ( nd)) 
    end if
  end subroutine polydiv_omp_init
  function polydiv_omp_lop ( adj, add, xx, yy) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: xx,yy 
    call adjnull (adj,add,xx,yy )
    call polydiv_omp_lop2(adj,add,xx,yy )
    stat=0
  end function polydiv_omp_lop
  subroutine polydiv_omp_lop2(adj,add,xx,yy)
    logical,intent(in) :: adj,add 
    real, dimension (:)  :: xx 
    real, dimension (:)  :: yy 
    integer  ia, ix, iy
    tt = 0.
    if ( adj) then
       do ix= nd, 1, -1  
          tt( ix) = yy( ix)
          do ia = 1, size( aa%lag)
             iy = ix + aa%lag( ia)
             if ( iy > nd) then
                cycle
             end if
             tt( ix) =                        tt( ix) -  aa%flt( ia) * tt( iy)
          end do
       end do
       xx =        xx + tt
    else
       do iy= 1, nd  
          tt( iy) = xx( iy)
          do ia = 1, size( aa%lag)
             ix = iy - aa%lag( ia)
             if ( ix < 1) then
                cycle
             end if
             tt( iy) =                        tt( iy) -  aa%flt( ia) * tt( ix)
          end do
       end do
       yy =        yy + tt
    end if
  end subroutine polydiv_omp_lop2
  subroutine polydiv_omp_close()
    deallocate( tt )
  end subroutine polydiv_omp_close
end module polydiv_omp
