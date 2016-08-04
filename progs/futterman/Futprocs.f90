!! Futprocs.x nt=128 Q=10 > out.H
program zyxabc
  ! nt must be a power of 2
  ! Auxiliary output futmatrix.H contains real(nt,nt) futterman wavelets at increasing delays.
  ! Pleasing values of Q for tutorial purposes are about nt/10.
  ! UNITS: Since units cancel from omega*t, I'm choosing sample units.
  use lapack95
  use f95_precision
  implicit none
  integer n, n1, nt
  real Q,precision, Qinv
  integer getch,putch
  call initpar()
  call doc('/work/FWI/WI/Futterman/futterman/Futprocs.rst')
  if (0>= getch('nt','d',nt )) then
     nt = 256
  end if
  if (0.ne.putch('Read  from param: #nt ','d',nt)) then
     call erexit('Trouble writing nt to history')
  end if
  if (0>= getch('Q','f',Q )) then
     Q = nt/10.
  end if
  if (0>= getch('Qinv','f',Qinv )) then
     Qinv = Q
  end if
  if (0>= getch('precision','f',precision )) then
     precision=40.
  end if
  if (0.ne.putch('Read  from param: #Q ','f',Q)) then
     call erexit('Trouble writing Q to history')
  end if
  if (0.ne.putch('Read  from param: #Precision ','f',precision)) then
     call erexit('Trouble writing precision to history')
  end if
  if (0.ne. putch('Q','f',Q )) then
     call erexit('Trouble writing  Q to history  ')
  end if
  n = nt
  call doit ( n, Q, Qinv, precision)
end program zyxabc
subroutine doit( n, Q, Qinv, precision)
  integer n, i, j, iseed
  real Q, sign, Qinv
  real power, phase, precision, prec_in, rcond, norm, norm_inv
  integer putch,info,ipiv(n),iwork(n),ipiv_inv(n),iwork_inv(n)
  real futs(n,n), invfuts(n,n), rwork(5*n), rwork_inv(5*n),idenfuts(n,n), cond, cond_inv ! Futterman operator
  double precision dfuts(n,n),dinvfuts(n,n),didenfuts(n,n)
  real modl(n) ! given model
  real data(n) ! data
  real data_inv(n) ! data
  real hmod(n) ! estimated model
  if (0 .ne. putch ( 'n1', 'i', n)) then
     call erexit('trouble writing to file ')
  end if
  if (0 .ne. putch ( 'd1', 'f', 0.004)) then
     call erexit('trouble writing to file ')
  end if
  if (0 .ne. putch ( 'n2', 'i', 5)) then
     call erexit('trouble writing to file ')
  end if
  if (0 .ne. putch ( 'n3', 'i', 1)) then
     call erexit('trouble writing to file ')
  end if
  call hclose()
  power=1.
  phase=1.
  prec_in=0.
  call makefuts(n, Q, futs, power, phase, prec_in)
  ! make Futterman Operator matrix for Figure
  do i=1,n
     do j=1,n
        futs(i,j) = futs(i,j) * sqrt( (i-1) +.1)
     end do
  end do
  call auxputch('d1','f', 1., 'futmatrix.H')
  call auxputch('o1','f', 0., 'futmatrix.H')
  call auxputch('d2','f', 1., 'futmatrix.H')
  call auxputch('o2','f', 0., 'futmatrix.H')
  call auxputch('esize','i', 4, 'futmatrix.H')
  call auxputch('n1','i', n, 'futmatrix.H')
  call auxputch('n2','i', n, 'futmatrix.H')
  call auxputch('label1','s','time','futmatrix.H')
  call auxputch('label2','s','traveltime distance','futmatrix.H')
  call srite( 'futmatrix.H', futs, n*n*4) ! Where does this get made?
  sign=1.
  do i=1,n
     modl(i) = 0.
  end do
  do i= n/10,n, n/10
     modl(i) = sign
     sign = -sign
  end do
  !iseed = 2015
  !do i=1, n {
  ! modl(i) = (rand01(iseed)-.5)**9
  ! }
  call srite('out', modl, 4*n) ! See the model
  ! Here I will explicitly change makefuts to make other operators, then use them.
  ! This will not ouput the new operators but will use them here.
  ! Plan: Give makefuts a definition of adjoint, AND change adj mult below to matmult.
  power=1.
  phase=1
  prec_in=0. ! Data
  call makefuts(n, Q, futs, power, phase, prec_in)
  ! make Futterman type Operators again!
  do j=1,n
     data(j) = 0.
  end do
  do j=1,n
     do i=1,n
        data(j) = data(j) + futs(j,i) * modl(i)
     end do
     ! Output the data
  end do
  call srite('out', data, 4*n) ! See the synthetic data.
  power=1.
  phase=-1.
  prec_in=0. ! Adjoint
  call makefuts(n, Qinv, futs, power, phase, prec_in)
  do i=1,n
     hmod(i) = 0.
  end do
  do j=1,n
     do i=1,n
        hmod(j) = hmod(j) + futs(j,i) * data(i)
     end do
  end do
  call srite('out', hmod, 4*n)
  power=0.
  phase=-1.
  prec_in=0. ! Unitary
  call makefuts(n, Qinv, futs, power, phase, prec_in)
  do i=1,n
     hmod(i) = 0.
  end do
  do j=1,n
     do i=1,n
        hmod(j) = hmod(j) + futs(j,i) * data(i)
     end do
  end do
  call srite('out', hmod, 4*n)
  power= -1
  phase=-1.
  prec_in=precision ! PseudoBest
  call makefuts(n, Qinv, invfuts, power, phase, prec_in)
!  power= 1
!  phase=1.
!  prec_in=,arw0 ! Forward
!  call makefuts(n, Qinv, futs, power, phase, prec_in)
!!  do j=1,n
!!     do i=j+1,n
!!        futs(j,i)=0
!!     end do
!!  end do
!!
!!  do j=1,n
!!     do i=j+1,n
!!        invfuts(j,i)=0
!!     end do
!!  end do
!
!  futs(1,1)=futs(2,2)
!  invfuts(1,1)=invfuts(2,2)
!
!  idenfuts=matmul(invfuts,futs)
!  data_inv=matmul(invfuts,data)
!
!  call srite('out', data_inv, 4*n)
!!  write(0,*) 'data_inv',data_inv
!
!  call srite( 'futmatrices.H', futs, n*n*4) ! Where does this get made?
!  call srite( 'futmatrices.H', invfuts, n*n*4) ! Where does this get made?
!  call srite( 'futmatrices.H', idenfuts, n*n*4) ! Where does this get made?
!
!  idenfuts(1,1)=idenfuts(2,2)
!
!!  call strcon('1','L','N',n,idenfuts,n,cond_inv,rwork_inv,iwork_inv,info)
!!  write(0,*) 'info: condition number BA=',cond_inv
!!
!!  call strcon('1','L','N',n,futs,n,cond,rwork,iwork,info)
!!  write(0,*) 'info: condition number A=',cond
!!
!  norm_inv=slange('1',n,n,idenfuts,n,rwork_inv)
!  write(0,*) 'L1 norm matrix BA=',norm_inv
!
!  call sgetrf(n,n,idenfuts,n,ipiv_inv,info)
!
!  call sgecon('1',n,idenfuts,n,norm_inv,cond_inv,rwork_inv,iwork_inv,info)
!  write(0,*) 'info: condition number BA=',cond_inv
!
!  norm=slange('1',n,n,futs,n,rwork)
!  write(0,*) 'L1 norm matrix A=',norm
!
!  call sgetrf(n,n,futs,n,ipiv,info)
!
!  call sgecon('1',n,futs,n,norm,cond,rwork,iwork,info)
!  write(0,*) 'info: condition number A=',cond
!
!!  call srite( 'futmatrices.H', idenfuts, n*n*4) ! Where does this get made?
!
!!  idenfuts=0.
!!  do i=1,n
!!     idenfuts(i,i)=1
!!  end do
!!  call sgetrs('N',n,n,futs,n,ipiv,idenfuts,n,info)
!!  idenfuts=sngl(didenfuts)
!!  call srite( 'futmatrices.H', idenfuts, n*n*4) ! Where does this get made?
!
!  call sgetrs('N',n,1,futs,n,ipiv,data,n,info)
!  call sgetrs('N',n,1,idenfuts,n,ipiv_inv,data_inv,n,info)
!
!  write(0,*) 'INFO',info,data
  do i=1,n
     hmod(i) = 0.
  end do
  do j=1,n
     do i=1,n
        hmod(j) = hmod(j) + invfuts(j,i) * data(i)
     end do
  end do
  call srite('out', hmod, 4*n)
!  call srite('out', data_inv, 4*n)
!
!  idenfuts=0
!  do i=1,n
!     idenfuts(i,i)=1
!  end do
!  call sgetrs('N',n,n,futs,n,ipiv,idenfuts,n,info)
!  call srite( 'futmatrices.H', idenfuts, n*n*4) ! Where does this get made?
!
!  call auxputch('d1','f', 1., 'futmatrices.H')
!  call auxputch('o1','f', 0., 'futmatrices.H')
!  call auxputch('d2','f', 1., 'futmatrices.H')
!  call auxputch('o2','f', 0., 'futmatrices.H')
!  call auxputch('esize','i', 4, 'futmatrices.H')
!  call auxputch('n1','i', n, 'futmatrices.H')
!  call auxputch('n2','i', n, 'futmatrices.H')
!  call auxputch('n3','i', 5, 'futmatrices.H')  
!  call auxputch('label1','s','time','futmatrices.H')
!  call auxputch('label2','s','traveltime distance','futmatrices.H')
  return
end subroutine doit
subroutine makefuts( n, Q, futs, power, phase, precision)
  integer n
  real Q, futs(n,n)
  real power, phase, precision
  complex cx(n)
  real omega(n)
  do i=1, n/2+1
     omega(i) = 3.14159265 * (i-1.)/n
  end do
  do i=1, n/2
     omega(n-i+1) = omega(i+1) ! Notice FFT domain symmetry.
  end do

  do j=1,n
!     write(0,*) -abs(omega(1:10)*j/Q)
     cx=exp(-abs(omega*j/Q)) ! exp(-|omega| t/Q )
     call kolmogoroff( n, cx, power, phase, precision)
     
    do i=1, n ! delay
        cx(i) = cx(i) * cexp(cmplx(0., 2.*3.14159265*(i-1.)/n * (j-1)))
     end do
     call ftu( -1., n, cx) ! Have time domain wavelet
     do i=1, n
        futs(i,j) = cx(i) * sqrt( (i-1)+.1)**power
     end do
     futs(1,j) = 0.
  end do
  return
end subroutine makefuts
subroutine kolmogoroff( n, cx, power, phase, precision)
  ! Spectral factorization.
  integer i, n
  ! input: cx = amplitude spectrum
  complex cx(n)
  ! output: cx = FT of min phase wavelet
  real power, phase, precision, gain, minlog, logamp
  real rr(n),cxin(n)
  cxin=cx
  do i= 1, n
     cx(i) = clog( cx(i)+ cmplx(1e-37,0.))
  end do

  call ftu( -1., n, cx)

  do i= 2, n/2
     ! Make it causal changing only the odd part.
     cx(i) = cx(i) * 2.
     cx(n-i+2) = 0.
  end do
  call ftu( +1., n, cx)

  if ( precision > 1. ) then
     ! Limit spectral range
     rr=-real(cx)
     minlog = minval(rr)
 
     where(rr.gt.(minlog+alog(precision))) 
        rr=minlog+alog(precision)
     endwhere

     call snap('rr.H', n, 1, rr)
     cx=cmplx(rr,phase*aimag(cx))

  else
     cx=cmplx(power*real(cx),phase*aimag(cx))
  end if
!  write(0,*) 'cx before cexp',cx
  cx=cexp(cx)
!  write(0,*) 'cx after cexp',cx
!  do i= 1, n
!     cx(i) = cexp( cx(i))
!  end do
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
