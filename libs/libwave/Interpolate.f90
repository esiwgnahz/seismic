! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Interpolate_mod

  use FD_types
  use ModelSpace_types

  implicit none

contains

  subroutine mksinc(sinc,lsinc,shift)
    !----------------------------------------------------------------
    implicit none
    integer :: n, lsinc
    real :: pi, shift, pos, fact, fsum
    real :: sinc(lsinc)
    !
    pi = 3.141592654
    fsum = 0.
    do n=1,lsinc
       if (n.le.int(float(lsinc/2+1)+shift)) then
          fact = (float(n-1)-shift)/float(lsinc/2)
       else
          fact = (float(lsinc-n)+shift)/float(lsinc/2)
       endif
       if (fact.le.0.) fact=0.
       pos = pi * (float(n-(lsinc/2+1)) - shift)
       if (abs(pos).gt.0.005) then
          sinc(n) = fact*sin(pos)/pos
       else
          sinc(n) = fact*(1.-pos*pos/6.)
       endif
       fsum = fsum + sinc(n)
    end do
    do n=1,lsinc
       sinc(n) = sinc(n) / fsum
    end do
  end subroutine mksinc

  subroutine Interpolate(model,bounds,genpar)
   
    type(GeneralParam)::              genpar
    type(ModelSpace) ::  model
    type(FDbounds)   ::        bounds
    integer :: i, j, k, l, lsinc, nmin, nmax
    real :: fdum
    parameter (lsinc=9)
    real, allocatable :: sinc(:), trace(:)
    !
    allocate(sinc(lsinc))
    nmin = min(bounds%nmin1-4,bounds%nmin2-4,bounds%nmin3-genpar%nbound) - lsinc
    nmax = max(bounds%nmax1+4,bounds%nmax2+4,bounds%nmax3+genpar%nbound) + lsinc
    allocate(trace(nmin:nmax))
    !
    ! z-axis
    call mksinc(sinc,lsinc,0.5)
    do k=bounds%nmin3-genpar%nbound,bounds%nmax3+genpar%nbound
       do j=bounds%nmin2-4,bounds%nmax2+4
          do i=bounds%nmin1-4,bounds%nmax1+4
             trace(i) = model%rho(i,j,k)
          end do
          do i=bounds%nmin1-4-lsinc,bounds%nmin1-4-1
             trace(i) = trace(bounds%nmin1-4)
          end do
          do i=bounds%nmax1+4+1,bounds%nmax1+4+lsinc
             trace(i) = trace(bounds%nmax1+4)
          end do
          do i=bounds%nmin1-4,bounds%nmax1+4
             fdum  = 0.
             do l=-lsinc/2,lsinc/2
                fdum = fdum + sinc(lsinc/2+1+l)*trace(i+l)
             end do
             model%rho2(i,j,k) = fdum
          end do
       end do
    end do
    !
    ! x-axis
    call mksinc(sinc,lsinc,0.5)
    do k=bounds%nmin3-genpar%nbound,bounds%nmax3+genpar%nbound
       do i=bounds%nmin1-4,bounds%nmax1+4
          do j=bounds%nmin2-4,bounds%nmax2+4
             trace(j) = model%rho2(i,j,k)
          end do
          do j=bounds%nmin2-4-lsinc,bounds%nmin2-4-1
             trace(j) = trace(bounds%nmin2-4)
          end do
          do j=bounds%nmax2+4+1,bounds%nmax2+4+lsinc
             trace(j) = trace(bounds%nmax2+4)
          end do
          do j=bounds%nmin2-4,bounds%nmax2+4
             fdum  = 0.
             do l=-lsinc/2,lsinc/2
                fdum = fdum + sinc(lsinc/2+1+l)*trace(j+l)
             end do
             model%rho2(i,j,k) = fdum
          end do
       end do
    end do
    !
    ! y-axis
    if ((bounds%nmax3+genpar%nbound).gt.1) then
       call mksinc(sinc,lsinc,0.5)
       do j=bounds%nmin2-4,bounds%nmax2+4
          do i=bounds%nmin1-4,bounds%nmax1+4
             do k=bounds%nmin3-genpar%nbound,bounds%nmax3+genpar%nbound
                trace(k) = model%rho2(i,j,k)
             end do
             do k=bounds%nmin3-genpar%nbound-lsinc,bounds%nmin3-genpar%nbound-1
                trace(k) = trace(bounds%nmin3-genpar%nbound)
             end do
             do k=bounds%nmax3+genpar%nbound+1,bounds%nmax3+genpar%nbound+lsinc
                trace(k) = trace(bounds%nmax3+genpar%nbound)
             end do
             do k=bounds%nmin3-genpar%nbound,bounds%nmax3+genpar%nbound
                fdum  = 0.
                do l=-lsinc/2,lsinc/2
                   fdum = fdum + sinc(lsinc/2+1+l)*trace(k+l)
                end do
                model%rho2(i,j,k) = fdum
             end do
          end do
       end do
    endif
    deallocate(trace,sinc)
  end subroutine Interpolate

end module Interpolate_mod
