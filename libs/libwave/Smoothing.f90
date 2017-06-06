! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Smoothing_mod

  use sep
  use Bandpass_mod
  use DataSpace_types
  use ModelSpace_types

  implicit none

  type SmoothingParam
     real :: percx
     real :: percz 
     real :: percy
     integer :: rect1
     integer :: rect2
     integer :: rect3
     integer :: n1
     integer :: n2
     integer :: n3
     
     real :: velref
     real :: freqref

     logical :: dox
     logical :: doy
     logical :: doz
  end type SmoothingParam

contains

  subroutine Init_SmoothingParam(smoothparam,mod,bpparam)
    type(SMoothingParam)::       smoothparam
    type(ModelSpace)    ::                   mod
    type(BandPassParam) ::                       bpparam
    real :: wavex,wavez,wavey

    call from_param('percx',smoothparam%percx,0.5)
    call from_param('percz',smoothparam%percz,0.5)
    call from_param('percy',smoothparam%percy,0.)

    smoothparam%velref=(maxval(mod%vel)-minval(mod%vel))/2+minval(mod%vel)
    smoothparam%freqref=bpparam%flo+(bpparam%fhi-bpparam%flo)/2

    wavex=smoothparam%percx*smoothparam%velref/smoothparam%freqref
    wavey=smoothparam%percy*smoothparam%velref/smoothparam%freqref
    wavez=smoothparam%percz*smoothparam%velref/smoothparam%freqref
   
    smoothparam%rect1=(wavez/mod%dz+1.5)/2
    smoothparam%rect2=(wavex/mod%dx+1.5)/2
    smoothparam%rect3=(wavey/mod%dy+1.5)/2

    smoothparam%dox=.true.;smoothparam%doy=.true.;smoothparam%doz=.true.
    
    if(wavez.le.mod%dz.or.wavez.eq.0.) smoothparam%doz=.false.
    if(wavex.le.mod%dx.or.wavex.eq.0.) smoothparam%dox=.false.
    if(wavey.le.mod%dy.or.wavey.eq.0.) smoothparam%doy=.false.
    
    smoothparam%n1=mod%nz
    smoothparam%n2=mod%nx
    smoothparam%n3=mod%ny
    write(0,*) 'INFO:------------------------------'
    write(0,*) 'INFO: Gradient smoothing parameters '
    write(0,*) 'INFO:------------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  percx   =  ',smoothparam%percx
    write(0,*) 'INFO:  percz   =  ',smoothparam%percz
    write(0,*) 'INFO:  percy   =  ',smoothparam%percy
    write(0,*) 'INFO:  dox     =  ',smoothparam%dox
    write(0,*) 'INFO:  doz     =  ',smoothparam%doz
    write(0,*) 'INFO:  doy     =  ',smoothparam%doy
    write(0,*) 'INFO:  rectz   =  ',smoothparam%rect1
    write(0,*) 'INFO:  rectx   =  ',smoothparam%rect2
    write(0,*) 'INFO:  recty   =  ',smoothparam%rect3
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  vref    =  ',smoothparam%velref
    write(0,*) 'INFO:  wref    =  ',smoothparam%freqref
    write(0,*) 'INFO:'
    write(0,*) 'INFO:------------------------------'
    
  end subroutine Init_SmoothingParam
  
  subroutine triangle2(smoothparam,uu)
    type(SMoothingParam):: smoothparam
    integer  :: i1,i2,i3,j,k
    real, dimension(smoothparam%n1,smoothparam%n2,smoothparam%n3) ::  uu
    real, dimension(:,:,:), allocatable ::  ss

    integer :: n1,n2,n3
    integer :: rect1,rect2,rect3
    logical  :: dox,doy,doz

    n1=smoothparam%n1
    n2=smoothparam%n2
    n3=smoothparam%n3
    rect1=smoothparam%rect1
    rect2=smoothparam%rect2
    rect3=smoothparam%rect3
    dox=smoothparam%dox
    doy=smoothparam%doy
    doz=smoothparam%doz
    
    allocate(ss(n1,n2,n3))

    if (dox) then
       do i3= 1, n3
          do i1= 1, n1
             call triangle( rect2, n2, uu(i1,:,i3), ss(i1,:,i3))
          end do
          uu=ss
       end do
    end if

    if (doz) then
       do i3= 1, n3
          do i2= 1, n2
             call triangle( rect1, n1, uu(:,i2,i3), ss(:,i2,i3))
          end do
       end do
       uu=ss
    end if

    if (doy) then
       do i2= 1, n2
          do i1= 1, n1
             call triangle( rect1, n1, uu(i1,i2,:), ss(i1,i2,:))
          end do
       end do
       uu=ss
    end if

    deallocate(ss)
    return
  end subroutine triangle2

  subroutine triangle( nr, n12, uu, vv)
    ! input:        nr      rectangle width (points) (Triangle base twice as wide.)
    ! input:        uu(m1,i2),i2=1,n12      is a vector of data.
    ! output:       vv(m1,i2),i2=1,n12      may be on top of uu
    integer nr,m1,n12, i,np,nq
    real, dimension(:)     :: uu,  vv
!    real, dimension(n12+nr-1)   :: pp
!    real, dimension(n12+2*nr-1) :: qq
!    real, dimension(n12)        :: tt
    real, dimension(:), allocatable   :: pp
    real, dimension(:), allocatable   :: qq
    real, dimension(:), allocatable   :: tt

    allocate(pp(n12+nr-1),qq(n12+2*nr-1),tt(n12))

    do i=1,n12
       qq(i) = uu(i)
    end do
    if ( n12 .eq. 1 ) then
       call copy( n12, qq, tt)
    else
       call boxconv( nr, n12, qq, pp)
       np = nr+n12-1
       call boxconv( nr, np , pp, qq)
       nq = nr+np-1
       do i= 1, n12
          tt(i) = qq(i+nr-1)
       end do
       do i= 1, nr-1
          ! fold back near end
          tt(i) = tt(i) + qq(nr-i)
       end do
       do i= 1, nr-1                                   ! fold back far end
          tt(n12-i+1) = tt(n12-i+1) + qq(n12+(nr-1)+i)
       end do
    end if
    do i=1,n12
       vv(i) = tt(i)
    end do
    deallocate(pp,qq,tt)
    return
  end subroutine triangle

  subroutine copy( n, xx, yy)
    integer i, n
    real xx(n), yy(n)
    do i= 1, n
       yy(i) = xx(i)
    end do
    return
  end subroutine copy

  subroutine boxconv( nb, nx, xx, yy)
    ! inputs:       nx,  xx(i), i=1,nx      the data
    !               nb                      the box length
    ! output:       yy(i),i=1,nx+nb-1       smoothed data
    integer nx, ny, nb, i
    real xx(nx), yy(1)
!    real bb(nx+nb)
    real, dimension(:), allocatable :: bb

    allocate(bb(nx+nb))
    if ( nb < 1 .or. nb > nx) then
       call erexit('boxconv')  ! "||" means .OR.
    end if
    ny = nx+nb-1
    do i= 1, ny
       bb(i) = 0.
    end do
    bb(1) = xx(1)
    do i= 2, nx
       bb(i) = bb(i-1) + xx(i)         ! make B(Z) = X(Z)/(1-Z)
    end do
    do i= nx+1, ny
       bb(i) = bb(i-1)
    end do
    do i= 1, nb
       yy(i) = bb(i)
    end do
    do i= nb+1, ny
       yy(i) = bb(i) - bb(i-nb)        ! make Y(Z) = B(Z)*(1-Z**nb)
    end do
    do i= 1, ny
       yy(i) = yy(i) / nb
    end do
    deallocate(bb)
    return
  end subroutine boxconv

end module Smoothing_mod
