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

  subroutine Interpolate(model,bounds)
   
    type(ModelSpace) ::  model
    type(FDbounds)   ::        bounds
    integer :: i, j, k, l, lsinc, nmin, nmax
    real :: fdum
    parameter (lsinc=9)
    real, allocatable :: sinc(:), trace(:)
    !
    allocate(sinc(lsinc))
    nmin = min(bounds%nmin1,bounds%nmin2,bounds%nmin3) - lsinc
    nmax = max(bounds%nmax1,bounds%nmax2,bounds%nmax3) + lsinc
    allocate(trace(nmin:nmax))
    !
    ! z-axis
    call mksinc(sinc,lsinc,0.5)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             trace(i) = model%rho(i,j,k)
          end do
          do i=bounds%nmin1-lsinc,bounds%nmin1-1
             trace(i) = trace(bounds%nmin1)
          end do
          do i=bounds%nmax1+1,bounds%nmax1+lsinc
             trace(i) = trace(bounds%nmax1)
          end do
          do i=bounds%nmin1,bounds%nmax1
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
    do k=bounds%nmin3,bounds%nmax3
       do i=bounds%nmin1,bounds%nmax1
          do j=bounds%nmin2,bounds%nmax2
             trace(j) = model%rho2(i,j,k)
          end do
          do j=bounds%nmin2-lsinc,bounds%nmin2-1
             trace(j) = trace(bounds%nmin2)
          end do
          do j=bounds%nmax2+1,bounds%nmax2+lsinc
             trace(j) = trace(bounds%nmax2)
          end do
          do j=bounds%nmin2,bounds%nmax2
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
    if (bounds%nmax3.gt.1) then
       call mksinc(sinc,lsinc,0.5)
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             do k=bounds%nmin3,bounds%nmax3
                trace(k) = model%rho2(i,j,k)
             end do
             do k=bounds%nmin3-lsinc,bounds%nmin3-1
                trace(k) = trace(bounds%nmin3)
             end do
             do k=bounds%nmax3+1,bounds%nmax3+lsinc
                trace(k) = trace(bounds%nmax3)
             end do
             do k=bounds%nmin3,bounds%nmax3
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
