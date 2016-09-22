module ncnhelicon  
! Nonstationary convolution, inverse to deconvolution.
  use adj_mod
  use nhelix
  use omp_lib
  implicit none
  type( nfilter), private :: aa
  !$OMP THREADPRIVATE(aa)
  contains
  subroutine ncnhelicon_init ( aa_in )
    type( nfilter)    :: aa_in 
    aa = aa_in
  end subroutine
  function ncnhelicon_lop ( adj, add, xx, yy) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: xx,yy 
    call adjnull (adj,add,xx,yy )
    call ncnhelicon_lop2(adj,xx,yy )
    stat=0
  end function 
  subroutine ncnhelicon_lop2(adj,xx,yy)
    logical,intent(in) :: adj 
    real, dimension (:)  :: xx 
    real, dimension (:)  :: yy 
    integer                         :: iy, ix, ia, ip, nx
    integer, dimension( :), pointer :: lag
    real,    dimension( :), pointer :: flt
    nx=size(xx)
    do iy = 1, size( yy)
      if ( associated( aa%mis)) then
        if ( aa%mis( iy)) then
          cycle
        end if
      end if
      ip = aa%pch( iy)
      lag => aa%hlx( ip)%lag
      flt => aa%hlx( ip)%flt
      do ia = 1, size( lag)
        ix = iy - lag( ia)
        if ((ix < 1).or.(ix>nx)) then
          cycle
        end if
        if ( adj) then
          xx(ix) =                        xx(ix) + yy(iy) * flt( ia)
        else
          yy(iy) =                        yy(iy) + xx(ix) * flt( ia)
        end if
      end do
    end do
  end subroutine 
  subroutine ncnhelicon_close()
  end subroutine 
end module 
