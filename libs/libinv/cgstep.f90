! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
!!$
!!$=head1 NAME
!!$
!!$cgstep  - on step of conjugate gradien step
!!$
!!$=head1 SYNOPSIS
!!$
!!$OPERATOR:C<ierr= cgstep(forget,x,g,rr,gg)>
!!$
!!$C<CLOSE:  call cgstep_close()>
!!$
!!$=head1 PARAMETERS
!!$
!!$=over 4
!!$
!!$=item forget - logical
!!$
!!$      Wheter or not to forget previous step
!!$
!!$=item x      - C<real(:)>
!!$
!!$      Model
!!$
!!$=item g      - C<real(:)>
!!$
!!$      Gradient
!!$
!!$=item rr      - C<real(:)>
!!$
!!$=item gg      - C<real(:)>
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$One step of conjugate gradient method
!!$
!!$=head1 SEE ALSO
!!$
!!$L<cgmeth>, L<cdstep>,L<solver_reg>,L<conjgrad>,L<partan>
!!$
!!$=head1 LIBRARY
!!$
!!$B<geef90>
!!$
!!$=cut
!!$
module cgstep_mod2   
  use sep
  use omp_lib
  implicit none
  real, dimension (:), allocatable, private   :: s, ss
  contains
  integer function cgstep( forget, x, g, rr, gg)
    real, dimension (:)  :: x, g, rr, gg
    logical              :: forget 
    double precision     :: sds, gdg, gds, determ, gdr, sdr, alfa,&
      & beta
    if ( .not. allocated (s)) then
      forget = .true.
      allocate ( s (size ( x))) 
      allocate (ss (size (rr)))
    end if
    if ( forget) then
      s = 0.
      ss = 0.
      beta = 0.d0    ! steepest descent
      if ( dot_product(gg, gg) .eq. 0 ) then
        call erexit('cgstep: grad vanishes identically')
      end if
      alfa = - sum( dprod( gg, rr)) / sum( dprod( gg, gg))
    else
      gdg = sum( dprod( gg, gg))       
      ! search plane by solving 2-by-2
      sds = sum( dprod( ss, ss))       
      !  G . (R - G*alfa - S*beta) = 0
      gds = sum( dprod( gg, ss))       
      !  S . (R - G*alfa - S*beta) = 0
      if ( gdg.eq.0. .or. sds.eq.0.) then
        cgstep = 1
        return
      end if
      determ = gdg * sds * max (1.d0 - (gds/gdg)*(gds/sds), 1.d-12)
      gdr = - sum( dprod( gg, rr))
      sdr = - sum( dprod( ss, rr))
      alfa = ( sds * gdr - gds * sdr ) / determ
      beta = (-gds * gdr + gdg * sdr ) / determ
    end if 
    s  = alfa *  g + beta *  s            ! update solution step
    ss = alfa * gg + beta * ss            ! update residual step
    x  =  x +  s                          ! update solution
    rr = rr + ss                          ! update residual
    forget = .false.
    cgstep = 0
  end function 
  subroutine cgstep_close ( )
    if ( allocated( s)) then
      deallocate( s, ss)
    end if
  end subroutine 
end module 
