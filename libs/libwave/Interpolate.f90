module Interpolate_mod

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

end module Interpolate_mod
