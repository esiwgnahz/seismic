! <in.H Rickdecon debubl=.06 ricker=.06 tresol=.01 shot=shot.H >out.H
program zyxabc
  use sep
  use omp_lib
  implicit none
  real, allocatable, dimension (:,:) :: data
  real, allocatable, dimension (:,:) :: datain
  real, allocatable, dimension (:) :: shot
  ! Inputs: d1=dt,debubl,ricker,tresol must have identical units of time (seconds)
  ! debubl=.06 =0 means no debubble
  ! ricker=.06 =0 means no Ricker compliance
  ! tresol=.01 time resolution, =0 for whitening decon
  ! Outputs: out.H is deconvolved data
  ! Outputs: shot(n) centered at (n+1)/2. Note n=2^^N > n1
  !
  ! subroutine rickdecon( shot,n, data,n1,n2, d1,debubl,ricker,tresol)
  !real dt, debubl, ricker, tresol
  integer n1,n2,n3,m1,i,j
  integer getch,hetch,putch,num_threads
  real d1, debubl, ricker, tresol, pink

  call initpar()
  call doc('/work/seismic/progs/logdecon/Rickerdecon.rst')
  if (0>= hetch('n1','d',n1 )) then
     call erexit('Trouble reading  n1 from history  ')
  end if
  if (0>= hetch('n2','d',n2 )) then
     n2 = 1
  end if
  if (0>= hetch('n3','d',n3 )) then
     n3 = 1  
  end if
  m1=1
  do while ( m1<n1)
     m1= 2*m1
  end do
  !m1 = 2*m1

  call from_param('num_threads',num_threads,omp_get_num_threads())

  if (0>= hetch('d1','f',d1 )) then
    call erexit('Trouble reading  d1 from history  ')
  end if
  if (0>= getch('debubl','f',debubl )) then
    debubl = .06
  end if
  if (0.ne.putch('Read  from param: #debubl ','f',debubl)) then
    call erexit('Trouble writing debubl to history')
  end if
  if (0>= getch('ricker','f',ricker )) then
    ricker = .06
  end if
  if (0.ne.putch('Read  from param: #ricker ','f',ricker)) then
    call erexit('Trouble writing ricker to history')
  end if
  if (0>= getch('tresol','f',tresol )) then
    tresol = .01
  end if
  if (0.ne.putch('Read  from param: #tresol ','f',tresol)) then
    call erexit('Trouble writing tresol to history')
  end if
  if (0>= getch('pink','f',pink )) then
    pink = 0.
  end if
  if (0.ne.putch('Read  from param: #pink ','f',pink)) then
    call erexit('Trouble writing pink to history')
  end if

  call omp_set_num_threads(num_threads)
  write(0,*) num_threads

  SHOTLOOP:do j=1,n3
     allocate (datain(n1,n2)) 
     write(0,*) 'INFO: starting reading panel',j
     call sreed( 'in', datain, 4*n1*n2)
     write(0,*) 'INFO: done reading panel',j
     !$OMP PARALLEL DO PRIVATE(i,shot)
     TRACELOOP:do i=1,n2
        !if (mod(i,n2/100).eq.0) write(0,*) 'Rickerdecon processing',i,'/',n2
        allocate (shot(m1))
        call rickdecon( shot,m1, datain(:,i),n1,1, d1,debubl,ricker,tresol,pink)
        deallocate(shot)
     end do TRACELOOP
    !$OMP END PARALLEL DO
     write(0,*) 'INFO: starting writing panel',j
     do i=1,n2
        call srite('out',datain(:,i),4*n1)
     end do
     write(0,*) 'INFO: done writing panel',j
     deallocate(datain)
  end do SHOTLOOP

  call hclose()
end program zyxabc

subroutine rickdecon(shot,m1, data,n1,n2, dt,debubl,ricker,tresol,pink)
integer m1, n1,n2,i1,i2
real shot(m1), data(n1,n2), dt,debubl,ricker,tresol&
  &,pink
complex cdata(m1)
complex cspec(m1)
complex Z
do i1=1,m1
  cspec(i1) = 0.0
end do
do i2=1,n2
  do i1=1,m1
    cdata(i1) = 0.
  end do
  do i1=1,n1
    cdata(i1) = data(i1,i2)
  end do
  call ftu(1., m1, cdata)
  do i1=1,m1
    cspec(i1) = cspec(i1) + cdata(i1) * conjg(cdata(i1))
  end do
end do
do i1=1,m1
  cspec(i1) = sqrt( real( cspec(i1)))
end do
call kolmogoroff( m1, cspec, dt, debubl,ricker,tresol)
! Spectral factorization.
do i2=1,n2
  do i1=1,m1
    cdata(i1) = 0.
  end do
  do i1=1,n1
    cdata(i1) = data(i1,i2)
  end do
  call ftu(1., m1, cdata)
  do i1=1,m1
    cdata(i1) = cdata(i1)/cspec(i1)
  end do
  call ftu(-1., m1, cdata)
  do i1=1,n1
    data(i1,i2) = cdata(i1)
  end do
end do
do i1=2,m1,2
  cspec(i1) = - cspec(i1)
  ! Moves the origin to the middle of the trace.
end do
call ftu(-1.,m1, cspec)
do i1=1,m1
  shot(i1) = cspec(i1)
  ! real values, big only near middle of vector
end do
return
end
subroutine kolmogoroff( n, cx, dt, debubl,ricker,tresol)
! Spectral factorization.
real dt, debubl,ricker,tresol, weight, tau,&
  & evn, odd
! Adapted from PVI, converted energy-->amplitude
integer i, n
! input: cx = amplitude spectrum
complex cx(n)
! output: cx = FT of min phase wavelet
average=0
do i= 1, n
! For precision sake, get logs surrounding zero.
  average = average + cx(i)
end do
average = average/n
do i= 1, n
  cx(i) = cx(i)/average
  ! Amplitude spectrum avg is unity, ready for logs
end do
do i= 1, n
  cx(i) = clog( cx(i) )
end do
call ftu( -1., n, cx)
do i= 2, n/2
! Make it causal changing only the odd part.
  cx(i) = cx(i) * 2.
  cx(n-i+2) = 0.
end do
!# BEGIN stuff to fiddle with inner lags of the odd part.
tau = dt
i=2
do while ( tau < debubl)
  weight = sin( .5 * 3.14159265 * tau/(debubl+1.e-20))**2
  cx(i) = cx(i) * weight
  cx(n-i+2) = cx(n-i+2) * weight
  i = i+1
  tau = tau + dt
end do
tau = dt
i=2
do while ( tau < tresol)
  weight = sin( .5 * 3.14159265 * tau/(tresol+1.e-20))**2
  cx(i) = cx(i) * weight
  cx(n-i+2) = cx(n-i+2) * weight
  i = i+1
  tau = tau + dt
end do
tau = dt
i=2
do while ( tau < ricker)
  weight = sin( .5 * 3.14159265 * tau/(ricker+1.e-20))**2
  evn = (cx(i) + cx(n-i+2))/2.
  odd = (cx(i) - cx(n-i+2))/2.
  odd = odd * weight
  cx(i) = evn + odd
  cx(n-i+2) = evn - odd
  i = i+1
  tau = tau + dt
end do
! END stuff to fiddle with inner lags of the odd part
call ftu( +1., n, cx)
do i= 1, n
  cx(i) = cexp( cx(i))
end do
return
end
subroutine ftu( signi, nx, cx )
! complex fourier transform with traditional scaling (FGDP)
!
! 1 nx signi*2*pi*i*(j-1)*(k-1)/nx
! cx(k) = -------- * sum cx(j) * e
! scale j=1 for k=1,2,...,nx=2**integer
!
! scale=1 for forward transform signi=1, otherwise scale=1/nx
integer nx, i, j, k, m, istep, pad2
real signi, arg
complex cx(nx), cmplx, cw, cdel, ct
if ( nx .ne. pad2(nx) ) then
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
end
