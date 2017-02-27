module helicon_mod
  ! Convolution, inverse to deconvolution.
  use adj_mod
  !             Requires the filter be causal with an implicit "1." at the onset.

  use omp_lib
  use helix

  implicit none
  type( filter), private :: aa 

  !$OMP THREADPRIVATE(aa)

contains
  subroutine helicon_mod_init ( aa_in )
    type( filter)    :: aa_in 
    aa = aa_in
  end subroutine helicon_mod_init
  function helicon_mod_lop ( adj, add, xx, yy) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: xx,yy 
    call adjnull (adj,add,xx,yy )
    call helicon_mod_lop2(adj,add,xx,yy )
    stat=0
  end function helicon_mod_lop
  subroutine helicon_mod_lop2(adj,add,xx,yy)
    logical,intent(in) :: adj,add 
    real, dimension (:)  :: xx 
    real, dimension (:)  :: yy 
    integer iy, ix, ia
    if ( adj) then
       ! zero lag
       xx =        xx + yy
    else
       yy =        yy + xx
    end if

    !$OMP PARALLEL DO PRIVATE(ia,iy,ix)
    do ia = 1, size( aa%lag)
       do iy = 1  + aa%lag( ia), size( yy)
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
    !$OMP END PARALLEL DO

  end subroutine helicon_mod_lop2
  subroutine helicon_mod_close()
  end subroutine helicon_mod_close
end module helicon_mod
