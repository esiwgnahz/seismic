! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Source_mod

  implicit none

contains

  subroutine comp2real(forw,carray,array,n)
    complex,dimension(:) :: carray
    real   ,dimension(:) :: array
    integer              :: i,j,n
    logical              :: forw

    j=1

    if (forw) then
       array=0
    else
       carray=0
    end if

    if (forw) then
       do i=1,2*n,2
          array(i)=real(carray(j))
          array(i+1)=aimag(carray(j))
          j=j+1
       end do
    else
       do i=1,2*n,2
          carray(j)=cmplx(array(i),array(i+1))
          j=j+1
       end do
    end if

  end subroutine comp2real

  subroutine SourceFct(wave,nt,dt,nhx)
    !
    integer :: nt, ntdata, nfft, ifft, it, j, nhx,k,i
    real :: wave(nt,nhx), dt
    real :: tpi, dw
    complex :: cdum
    complex, allocatable :: trace(:)
    real,    allocatable :: tracer(:)
    !
    tpi=6.283185307
    ntdata = 2**(int(log(float(nt+100))/log(2.0))+1)
    nfft = ntdata/2+1
    allocate(trace(ntdata))
    allocate(tracer(2*ntdata))

    dw = tpi/float(ntdata)/dt
    do j=1,nhx
       do it=1,ntdata
          if (it.le.nt) then
             trace(it) = cmplx(wave(it,j),0.)
          else
             trace(it) = cmplx(0.,0.)
          endif
       end do

       call comp2real(.True.,trace,tracer,ntdata)   
       call FOUR1(tracer,ntdata,+1)
       call comp2real(.False.,trace,tracer,ntdata)

       do it=1,nfft
          cdum = trace(it)*csqrt(cmplx(0.,-float(it-1)*dw))
          trace(it) = cdum
          ifft=ntdata-it+2
          if (ifft.gt.nfft .and. ifft.le.ntdata) &
          &   trace(ifft) = conjg(cdum)
       end do
       call comp2real(.True.,trace,tracer,ntdata) 
       call FOUR1(tracer,ntdata,-1)
       call comp2real(.False.,trace,tracer,ntdata) 

       do it=1,nt
          wave(it,j) = real(trace(it))/float(ntdata)
       end do
    end do
    deallocate(trace,tracer)
  end subroutine SourceFct  

  subroutine SourceRandom(wave,nt,dt,nhx,phase)
    !
    integer :: nt, ntdata, nfft, ifft, it, j, nhx,k,i
    real :: wave(nt,nhx), dt
    real :: tpi, dw, phase
    complex :: cdum
    complex, allocatable :: trace(:)
    real,    allocatable :: tracer(:)
    !
    tpi=6.28318530717959
    ntdata = 2**(int(log(float(nt+100))/log(2.0))+1)
    nfft = ntdata/2+1
    allocate(trace(ntdata))
    allocate(tracer(2*ntdata))

    dw = tpi/float(ntdata)/dt
    do j=1,nhx
       do it=1,ntdata
          if (it.le.nt) then
             trace(it) = cmplx(wave(it,j),0.)
          else
             trace(it) = cmplx(0.,0.)
          endif
       end do

       call comp2real(.True.,trace,tracer,ntdata)   
       call FOUR1(tracer,ntdata,+1)
       call comp2real(.False.,trace,tracer,ntdata)

       do it=1,nfft
!          cdum = trace(it)*csqrt(cmplx(0.,-float(it-1)*dw))
          cdum = trace(it)*cexp(cmplx(0,phase))
          trace(it) = cdum
          ifft=ntdata-it+2
          if (ifft.gt.nfft .and. ifft.le.ntdata) &
          &   trace(ifft) = conjg(cdum)
       end do
       call comp2real(.True.,trace,tracer,ntdata) 
       call FOUR1(tracer,ntdata,-1)
       call comp2real(.False.,trace,tracer,ntdata) 

       do it=1,nt
          wave(it,j) = real(trace(it))/float(ntdata)
       end do
    end do
    deallocate(trace,tracer)
  end subroutine SourceRandom
  !
  subroutine SourceFct_inverse(wave,nt,dt,nhx)
    !
    integer :: nt, ntdata, nfft, ifft, it, j, nhx,k,i
    real :: wave(nt,nhx), dt
    real :: tpi, dw
    complex :: cdum
    complex, allocatable :: trace(:)
    real,    allocatable :: tracer(:)
    !
    tpi=6.283185307
    ntdata = 2**(int(log(float(nt+100))/log(2.0))+1)
    nfft = ntdata/2+1
    allocate(trace(ntdata))
    allocate(tracer(2*ntdata))

    dw = tpi/float(ntdata)/dt
    do j=1,nhx
       do it=1,ntdata
          if (it.le.nt) then
             trace(it) = cmplx(wave(it,j),0.)
          else
             trace(it) = cmplx(0.,0.)
          endif
       end do

       call comp2real(.True.,trace,tracer,ntdata)   
       call FOUR1(tracer,ntdata,+1)
       call comp2real(.False.,trace,tracer,ntdata)

       ! Do not start at the zero frequency
       do it=2,nfft
          cdum = trace(it)/csqrt(cmplx(0.,float(it-1)*dw))
          trace(it) = cdum
          ifft=ntdata-it+2
          if (ifft.gt.nfft .and. ifft.le.ntdata) &
          &   trace(ifft) = conjg(cdum)
       end do
       call comp2real(.True.,trace,tracer,ntdata) 
       call FOUR1(tracer,ntdata,-1)
       call comp2real(.False.,trace,tracer,ntdata) 

       do it=1,nt
          wave(it,j) = real(trace(it))/float(ntdata)
       end do
    end do
    deallocate(trace,tracer)
  end subroutine SourceFct_inverse
  !
  subroutine FOUR1(x,n,isign)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    !      FOUR1 calculates either the forward or inverse discrete Fourier
    !  transform of a complex series.
    !               N
    !       F(K) = SUM X(J)*EXP(I*ISIGN*2*PI*(K-1)*(J-1)/N)    K=1,N
    !              J=1
    !
    !  Input  --
    !      X     Complex vector of dimension .GE. N.
    !      N     Length of the vector X which must be an exact power of 2.
    !      ISIGN .EQ. -1 for computation of forward transform
    !            .EQ. +1 for computation of inverse transform
    !
    !  Output  --
    !      X     The transform is stored in the position of the original
    !            series with point 1 the zero frequency and point N/2+1
    !            the Nyquist frequency.  The real part contains the cosine
    !            series and is symmetric about the Nyquist.  The imaginary
    !            part contains the sine series and is asymmetric about the
    !            Nyquist. Positive frequencies are in the range 1 to N/2+1
    !            and negative frequencies are in the range N/2+2 to N.
    !
    !  Notes  --
    !      The scale factor of 1/N is not  applied.
    !
    !      The vector X is overwritten by the transform so it is destroyed.
    !
    !      The vector X is single precision complex (COMPLEX*8), but
    !      the critical internal calculations are performed in double
    !      precision arithmetic.
    !
    !      This is a slightly modified version of the program in
    !     "Numerical Recipes" by Press et al.
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

    integer n,i,j,nn,m,mmax,istep,isign
    real*4 x(2*n),tempr,tempi
    real*8 wr,wi,wpr,wpi,wtemp,theta
    nn=2*n
    j=1

    do i=1,nn,2
       if(j.gt.i) then
          tempr=x(j)
          tempi=x(j+1)
          x(j)=x(i)
          x(j+1)=x(i+1)
          x(i)=tempr
          x(i+1)=tempi
       endif
       m=nn/2
       do while ((m.ge.2).and.(j.gt.m))
          j=j-m
          m=m/2
       end do
       j=j+m
    end do
    mmax=2
    do while (nn.gt.mmax)
       istep=2*mmax
       theta=6.28318530717959d0/dfloat(isign*mmax)
       wpr=-2.0*dsin(0.5d0*theta)**2
       wpi=dsin(theta)
       wr=1.0d0
       wi=0.0d0
       do m=1,mmax,2
          do i=m,nn,istep
             j=i+mmax
             tempr=sngl(wr)*x(j)-sngl(wi)*x(j+1)
             tempi=sngl(wr)*x(j+1)+sngl(wi)*x(j)
             x(j)=x(i)-tempr
             x(j+1)=x(i+1)-tempi
             x(i)=x(i)+tempr
             x(i+1)=x(i+1)+tempi
          end do
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
       end do
       mmax=istep
    end do

  end subroutine FOUR1

!  subroutine invhalfdifa( adj, add, n, xx,    yy  )
!    integer n2, i,     adj, add, n
!    real omega,                     xx(n), yy(n)
!    complex cz, cv(n*2)
!    n2=1
!
!    do while (n2<n)
!       n2=2*n2
!    end do
!    do i= 1, n2  
!       cv(i) = 0.
!    end do
!    do i= 1, n
!       if ( adj .eq. 0) then
!          cv(i) = xx(i)
!       else
!          cv(i) = yy(i)
!       end if
!    end do
!    call adjnull(adj, add,xx,n,yy,n)
!    call ftu( +1., n2, cv)
!    do i= 1, n2  
!       omega = (i-1.) * 2.*3.14159265 / n2
!       cz = csqrt( 1.0000000001 - cexp( cmplx( 0., omega)))
!       if ( adj .ne. 0) then
!          cz = conjg( cz)
!       end if
!       cv(i) = cv(i) / cz
!    end do
!    call ftu( -1., n2, cv)
!    do i= 1, n
!       if ( adj .eq. 0) then
!          yy(i) = yy(i) + cv(i)
!       else
!          xx(i) = xx(i) + cv(i)
!       end if
!    end do
!    return
!  end subroutine invhalfdifa
!
!  subroutine ftu( signi, nx, cx )
!    integer        nx, i, j, k, m, istep,a
!    real        signi, scale, arg
!    complex        cx(nx), cmplx, cw, cdel, ct
!    if ( nx .ne. pad2(nx) ) then
!       call erexit('ftu: nx not a power of 2')
!    end if
!    scale = 1. / sqrt( 1.*nx)
!    do i= 1, nx
!       cx(i) = cx(i) * scale
!    end do
!    j = 1
!    k = 1
!    do i= 1, nx  
!       if (i<=j) then
!          ct = cx(j)
!          cx(j) = cx(i)
!          cx(i) = ct
!       end if
!       m = nx/2
!       do while  (j>m .and. m>1)
!          j = j-m
!          m = m/2
!       end do
!       ! "&" means .AND.
!       j = j+m
!    end do
!    a=0
!    do while (a.eq.0)
!       istep = 2*k
!       cw = 1.
!       arg = signi*3.14159265/k
!       cdel = cmplx( cos(arg), sin(arg))
!       do m= 1, k  
!          do i= m, nx, istep
!             ct=cw*cx(i+k)
!             cx(i+k)=cx(i)-ct
!             cx(i)=cx(i)+ct
!          end do
!          cw = cw * cdel
!       end do
!       k = istep
!       if (k>=nx) then
!          exit
!       end if
!    end do
!    return
!  end subroutine ftu

  integer function pad2(nx)
    integer nx,n
    n=1
    do while (n < nx)
       n=n*2
    end do
    pad2=n
    return
  end function pad2

!  subroutine adjnull( adj, add, x, nx,  y, ny )
!    integer ix, iy,     adj, add,    nx,     ny
!    real                          x( nx), y( ny )
!    if ( add .eq. 0 ) then
!       if ( adj .eq. 0 ) then
!          do iy= 1, ny
!             y(iy) = 0.
!          end do
!       else
!          do ix= 1, nx
!             x(ix) = 0. 
!          end do
!       end if
!    end if
!    return
!  end subroutine adjnull

end module Source_mod

