module nchelicon  
! Non-causal convolution, no inverse.
!  Do not confuse with ordinary helicon
  use adj_mod
  use helix
  use omp_lib
  implicit none
  type( filter), private :: aa
  !$OMP THREADPRIVATE(aa)
 
contains
  subroutine nchelicon_init ( aa_in )
    type( filter)    :: aa_in 
    aa = aa_in
  end subroutine nchelicon_init
  function nchelicon_lop ( adj, add, xx, yy) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: xx,yy 
    call adjnull (adj,add,xx,yy )
    call nchelicon_lop2(adj,xx,yy )
    stat=0
  end function nchelicon_lop
  subroutine nchelicon_lop2(adj,xx,yy)
    logical,intent(in) :: adj 
    real, dimension (:)  :: xx 
    real, dimension (:)  :: yy 
    integer iy, ix, ia
    do ia = 1, size( aa%lag)
      do iy = 1  + max(aa%lag( ia),0),size( yy)+min(aa%lag( ia),0)
        if ( associated( aa%mis)) then
          if ( aa%mis( iy)) then
            cycle
          end if
        end if
        ix = iy - aa%lag( ia)
        if ( adj) then
          xx(ix) =                        xx(ix) + yy(iy) * aa%flt(ia)
        else
          yy(iy) =                        yy(iy) + xx(ix) * aa%flt(ia)
        end if
      end do
    end do
  end subroutine nchelicon_lop2
  subroutine nchelicon_close()
  end subroutine nchelicon_close
end module nchelicon

