! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module nhelicon2_mod  
  ! Nonstationary convolution, inverse to deconvolution.
  use adj_mod
  !             Requires the filter be causal with an implicit "1." at the onset.
  use nhelix

  use omp_lib

  implicit none
  type( nfilter), private :: aa 

  !$OMP THREADPRIVATE(aa)

contains

  subroutine nhelicon2_mod_init ( aa_in )
    type( nfilter)    :: aa_in 
    aa = aa_in
  end subroutine nhelicon2_mod_init
  function nhelicon2_mod_lop ( adj, add, xx, yy) result(stat)
    integer            :: stat 
    logical,intent(in) :: adj,add 
    real,dimension(:)  :: xx,yy 
    call adjnull (adj,add,xx,yy )
    call nhelicon2_mod_lop2(adj,add,xx,yy )
    stat=0
  end function nhelicon2_mod_lop
  subroutine nhelicon2_mod_lop2(adj,add,xx,yy)
    logical,intent(in) :: adj,add 
    real, dimension (:)  :: xx 
    real, dimension (:)  :: yy 
    integer                         :: iy, ix, ia, ip
    integer, dimension( :), pointer :: lag
    real,    dimension( :), pointer :: flt
    if ( adj) then
       ! zero lag
       xx =        xx + yy
    else
       yy =        yy + xx
    end if
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
          if ( ix < 1) then
             cycle
          end if
          if ( adj) then
             xx(ix) =                        xx(ix) + yy(iy) * flt( ia)
          else
             yy(iy) =                        yy(iy) + xx(ix) * flt( ia)
          end if
       end do
    end do
  end subroutine nhelicon2_mod_lop2
  subroutine nhelicon2_mod_close()
  end subroutine nhelicon2_mod_close
end module nhelicon2_mod
