module Futterman

  implicit none
  
  integer, private :: nt
  real,    private :: dt,ot,dtw

  contains

subroutine futterman_init(nt_init,dt_init,ot_init,dtw_init)
  integer :: nt_init
  real    :: dt_init, ot_init, dtw_init
  nt=nt_init
  dt=dt_init
  ot=ot_init
  dtw=dtw_init
end subroutine futterman_init

subroutine makefuts( Q, futs, power, phase, precision)
  real Q, futs(nt,nt)
  real power, phase, precision, tmute
  complex cx(nt)
  real omega(nt),weight
  integer i,j
  do i=1, nt/2+1
     omega(i) = 3.14159265 * (i-1.)/nt
  end do
  do i=1, nt/2
     omega(nt-i+1) = omega(i+1) ! Notice FFT domain symmetry.
  end do
  do j=1,nt
     cx=exp(-abs(omega*j/Q)) ! exp(-|omega| t/Q )

     call kolmogoroff( cx, power, phase, precision)

     do i=1, nt ! delay
        cx(i) = cx(i) * cexp(cmplx(0., 2.*3.14159265*(i-1.)/nt * (j-1)))
     end do

     call ftu( -1., nt, cx) ! Have time domain wavelet
     !     futs(:,j)=cx
     do i=1, nt
        weight=(i-1)*dt-ot+dtw ! optimal weight t-te+dt (see Jon2 in sep149)
                               ! but here te is always zero by construction
!        weight=(i-1)+.1 ! optimal weight t-te+dt (see Jon2 in sep149)
                               ! but here te is always zero by construction
!        futs(i,j) = cx(i) * sqrt(weight)**power ! Need to see if we need sqrt or not
        futs(i,j) = cx(i) * sqrt(weight)**power ! Need to see if we need sqrt or not
     end do
     futs(1,j) = 0.
  end do
  return
end subroutine makefuts

subroutine kolmogoroff( cx, power, phase, precision)
  ! Spectral factorization.
  integer i
  ! input: cx = amplitude spectrum
  complex cx(nt)
  ! output: cx = FT of min phase wavelet
  real power, phase, precision, gain, minlog, logamp
  real rr(nt)

  cx=clog(cx+cmplx(1e-37,0))

  call ftu( -1., nt, cx)

  do i= 2, nt/2
     ! Make it causal changing only the odd part.
     cx(i) = cx(i) * 2.
     cx(nt-i+2) = 0.
  end do

  call ftu( +1., nt, cx)

  if ( precision > 1. ) then
     ! Limit spectral range
     rr=-real(cx)
     minlog = minval(rr)
 
     where(rr.gt.(minlog+alog(precision))) 
        rr=minlog+alog(precision)
     endwhere
     cx=cmplx(rr,phase*aimag(cx))

  else
     cx=cmplx(power*real(cx),phase*aimag(cx))
  end if

  cx=cexp(cx)

  return
end subroutine kolmogoroff

subroutine ftu( signi, nx, cx )
  ! complex fourier transform with traditional scaling (FGDP)
  !
  ! 1 nx signi*2*pi*i*(j-1)*(k-1)/nx
  ! cx(k) = -------- * sum cx(j) * e
  ! scale j=1 for k=1,2,...,nx=2**integer
  !
  ! scale=1 for forward transform signi=1, otherwise scale=1/nx
  integer nx, i, j, k, m, istep
  real signi, arg
  complex cx(nx), cmplx, cw, cdel, ct
  i=1
  do while ( i<nx)
     i=2*i
  end do
  if ( i .ne. nx ) then
     call erexit('ftu: nx not a power of 2')
  end if
  do i= 1, nx
     if ( signi<0.) then
        cx(i) = cx(i) / nx
     end if
  end do
  j = 1
  k = 1
  do i= 1, nx
     if (i<=j) then
        ct = cx(j)
        cx(j) = cx(i)
        cx(i) = ct
     end if
     m = nx/2
     do while (j>m .and. m>1)
        j = j-m
        m = m/2
     end do
     ! "&&" means .AND.
     j = j+m
  end do
  do
     istep = 2*k
     cw = 1.
     arg = signi*3.14159265/k
     cdel = cmplx( cos(arg), sin(arg))
     do m= 1, k
        do i= m, nx, istep
           ct=cw*cx(i+k)
           cx(i+k)=cx(i)-ct
           cx(i)=cx(i)+ct
        end do
        cw = cw * cdel
     end do
     k = istep
     if (k>=nx) then
        exit
     end if
  end do
  return
end subroutine ftu

end module Futterman
