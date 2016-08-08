module Boundary_mod

  use GeneralParam_types
  use ModelSpace_types
  use Interpolate_mod
  use FD_types
  
  implicit none
  
contains
  
  subroutine Higdon(dt,mod,bounds,hig)
    type(HigdonParam) :: hig
    type(FDbounds)    :: bounds
    type(ModelSpace)  :: mod
    real              :: dt
    !
    integer :: i, j, k
    real :: beta1, beta2, weight, dxi, dyi, dzi
    real :: nu
    real :: qx, qt, qxt, rx, rt, rxt, fdum
    
    beta1 = 1.0
    beta2 = 0.7
    weight = 0.5
    
    dxi=1./mod%dx
    dyi=1./mod%dy
    dzi=1./mod%dz
    !
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          nu = mod%vel(bounds%nmin1,j,k)*dt*dzi
          call Higdon_coeff(nu,beta1,beta2,weight,hig%gz(1,j,k))
          nu = mod%vel(bounds%nmax1,j,k)*dt*dzi
          call Higdon_coeff(nu,beta1,beta2,weight,hig%gz(9,j,k))
       end do
    end do
    do k=bounds%nmin3,bounds%nmax3
       do i=bounds%nmin1,bounds%nmax1
          nu = mod%vel(i,bounds%nmin2,k)*dt*dxi
          call Higdon_coeff(nu,beta1,beta2,weight,hig%gx(1,i,k))
          nu = mod%vel(i,bounds%nmax2,k)*dt*dxi
          call Higdon_coeff(nu,beta1,beta2,weight,hig%gx(9,i,k))
       end do
    end do
    if (bounds%nmax3.gt.1) then
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             nu = mod%vel(i,j,bounds%nmin3)*dt*dyi
             call Higdon_coeff(nu,beta1,beta2,weight,hig%gy(1,i,j))
             nu = mod%vel(i,j,bounds%nmax3)*dt*dyi
             call Higdon_coeff(nu,beta1,beta2,weight,hig%gy(9,i,j))
          end do
       end do
    endif
    !
  end subroutine Higdon
  !
  !--------------------------------------------------------------------
  subroutine Higdon_coeff(nu,beta1,beta2,weight,coeff)
    real :: nu, beta1, beta2, weight, coeff(8)
    real :: qx, qt, qxt, rx, rt, rxt, fdum
    !
    fdum = (beta1+nu)*(1-weight)
    qx = (weight*(beta1+nu)-nu)/fdum
    qt = (weight*(beta1+nu)-beta1)/fdum
    qxt = weight/(weight-1)
    fdum = (beta2+nu)*(1-weight)
    rx = (weight*(beta2+nu)-nu)/fdum
    rt = (weight*(beta2+nu)-beta2)/fdum
    rxt = weight/(weight-1)
    coeff(1) = qx+rx
    coeff(2) = qx*rx
    coeff(3) = qt+rt
    coeff(4) = qx*rt+qt*rx+qxt+rxt
    coeff(5) = qx*rxt+rx*qxt
    coeff(6) = qt*rt
    coeff(7) = qt*rxt+rt*qxt
    coeff(8) = qxt*rxt
    !
  end subroutine Higdon_coeff

  subroutine Boundary0(genpar,bounds,u,mod,hig)
    !
    type(GeneralParam) :: genpar
    type(ModelSpace)   :: mod
    type(FDbounds)     :: bounds
    type(HigdonParam)  :: hig

    real :: u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &         bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)

    integer :: i, j, k, nz0
    integer :: fscale0
    real    :: fscale, fdum

    nz0 = 1
    fscale0 = 3.
    !
    ! Update the taper zone and update the regions used for operator padding
    ! Since everything is outside the main grid, constant density applies
    !
    ! Axis 1 ----------------------
    !
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          ! Update the top
          ! Taper
          if (genpar%surf_type.le.0) then
             do i=nz0-1,bounds%nmin1,-1
                fscale = float(nz0-i)/float(nz0-1-bounds%nmin1)/fscale0
                fdum = hig%gz(1,j,k)*u(i+1,j,k,3)+ &
                &   hig%gz(2,j,k)*u(i+2,j,k,3)+ &
                &   hig%gz(3,j,k)*u(i  ,j,k,2)+ &
                &   hig%gz(4,j,k)*u(i+1,j,k,2)+ &
                &   hig%gz(5,j,k)*u(i+2,j,k,2)+ &
                &   hig%gz(6,j,k)*u(i  ,j,k,1)+ &
                &   hig%gz(7,j,k)*u(i+1,j,k,1)+ &
                &   hig%gz(8,j,k)*u(i+2,j,k,1)
                u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
             end do
             ! Operator padding
             do i=bounds%nmin1-1,bounds%nmin1-4,-1
                fdum = hig%gz(1,j,k)*u(i+1,j,k,3)+ &
                &   hig%gz(2,j,k)*u(i+2,j,k,3)+ &
                &   hig%gz(3,j,k)*u(i  ,j,k,2)+ &
                &   hig%gz(4,j,k)*u(i+1,j,k,2)+ &
                &   hig%gz(5,j,k)*u(i+2,j,k,2)+ &
                &   hig%gz(6,j,k)*u(i  ,j,k,1)+ &
                &   hig%gz(7,j,k)*u(i+1,j,k,1)+ &
                &   hig%gz(8,j,k)*u(i+2,j,k,1)
                u(i,j,k,3) = -fdum
             end do
          end if
          !
          ! Update the bottom operator pad
          ! Taper
          do i=mod%nz+1,bounds%nmax1
             fscale = float(i-mod%nz)/float(bounds%nmax1-mod%nz-1)/fscale0
             fdum = hig%gz( 9,j,k)*u(i-1,j,k,3)+ &
             &   hig%gz(10,j,k)*u(i-2,j,k,3)+ &
             &   hig%gz(11,j,k)*u(i  ,j,k,2)+ &
             &   hig%gz(12,j,k)*u(i-1,j,k,2)+ &
             &   hig%gz(13,j,k)*u(i-2,j,k,2)+ &
             &   hig%gz(14,j,k)*u(i  ,j,k,1)+ &
             &   hig%gz(15,j,k)*u(i-1,j,k,1)+ &
             &   hig%gz(16,j,k)*u(i-2,j,k,1)
             u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
          end do
          ! Operator padding
          do i=bounds%nmax1+1,bounds%nmax1+4
             fdum = hig%gz( 9,j,k)*u(i-1,j,k,3)+ &
             &   hig%gz(10,j,k)*u(i-2,j,k,3)+ &
             &   hig%gz(11,j,k)*u(i  ,j,k,2)+ &
             &   hig%gz(12,j,k)*u(i-1,j,k,2)+ &
             &   hig%gz(13,j,k)*u(i-2,j,k,2)+ &
             &   hig%gz(14,j,k)*u(i  ,j,k,1)+ &
             &   hig%gz(15,j,k)*u(i-1,j,k,1)+ &
             &   hig%gz(16,j,k)*u(i-2,j,k,1)
             u(i,j,k,3) = -fdum
          end do
       end do
    end do
    !
    ! Axis 2 ----------------------
    !
    do k=bounds%nmin3,bounds%nmax3
       !
       ! Update the left side operator pad
       ! Taper
       do j=0,bounds%nmin2,-1
          fscale = float(1-j)/float(-bounds%nmin2)/fscale0
          do i=bounds%nmin1,bounds%nmax1
             fdum = hig%gx(1,i,k)*u(i,j+1,k,3)+ &
             &   hig%gx(2,i,k)*u(i,j+2,k,3)+ &
             &   hig%gx(3,i,k)*u(i,j  ,k,2)+ &
             &   hig%gx(4,i,k)*u(i,j+1,k,2)+ &
             &   hig%gx(5,i,k)*u(i,j+2,k,2)+ &
             &   hig%gx(6,i,k)*u(i,j  ,k,1)+ &
             &   hig%gx(7,i,k)*u(i,j+1,k,1)+ &
             &   hig%gx(8,i,k)*u(i,j+2,k,1)
             u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
          end do
       end do
       ! Operator padding
       do j=bounds%nmin2-1,bounds%nmin2-4,-1
          do i=bounds%nmin1,bounds%nmax1
             fdum = hig%gx(1,i,k)*u(i,j+1,k,3)+ &
             &   hig%gx(2,i,k)*u(i,j+2,k,3)+ &
             &   hig%gx(3,i,k)*u(i,j  ,k,2)+ &
             &   hig%gx(4,i,k)*u(i,j+1,k,2)+ &
             &   hig%gx(5,i,k)*u(i,j+2,k,2)+ &
             &   hig%gx(6,i,k)*u(i,j  ,k,1)+ &
             &   hig%gx(7,i,k)*u(i,j+1,k,1)+ &
             &   hig%gx(8,i,k)*u(i,j+2,k,1)
             u(i,j,k,3) = -fdum
          end do
       end do
       !
       ! Update the right side operator pad
       ! Taper
       do j=mod%nx+1,bounds%nmax2
          fscale = float(j-mod%nx)/float(bounds%nmax2-mod%nx-1)/fscale0
          do i=bounds%nmin1,bounds%nmax1
             fdum = hig%gx( 9,i,k)*u(i,j-1,k,3)+ &
             &   hig%gx(10,i,k)*u(i,j-2,k,3)+ &
             &   hig%gx(11,i,k)*u(i,j  ,k,2)+ &
             &   hig%gx(12,i,k)*u(i,j-1,k,2)+ &
             &   hig%gx(13,i,k)*u(i,j-2,k,2)+ &
             &   hig%gx(14,i,k)*u(i,j  ,k,1)+ &
             &   hig%gx(15,i,k)*u(i,j-1,k,1)+ &
             &   hig%gx(16,i,k)*u(i,j-2,k,1)
             u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
          end do
       end do
       do j=bounds%nmax2+1,bounds%nmax2+4
          do i=bounds%nmin1,bounds%nmax1
             fdum = hig%gx( 9,i,k)*u(i,j-1,k,3)+ &
             &   hig%gx(10,i,k)*u(i,j-2,k,3)+ &
             &   hig%gx(11,i,k)*u(i,j  ,k,2)+ &
             &   hig%gx(12,i,k)*u(i,j-1,k,2)+ &
             &   hig%gx(13,i,k)*u(i,j-2,k,2)+ &
             &   hig%gx(14,i,k)*u(i,j  ,k,1)+ &
             &   hig%gx(15,i,k)*u(i,j-1,k,1)+ &
             &   hig%gx(16,i,k)*u(i,j-2,k,1)
             u(i,j,k,3) = -fdum
          end do
       end do
    end do
    !
    ! Axis 3 ----------------------
    !
    if (genpar%nbound.gt.0) then
       !
       ! Update the front side operator pad
       ! Taper
       do k=0,bounds%nmin3,-1
          fscale = float(1-k)/float(-bounds%nmin3)/fscale0
          do j=bounds%nmin2,bounds%nmax2
             do i=bounds%nmin1,bounds%nmax1
                fdum = hig%gy(1,i,j)*u(i,j,k+1,3)+ &
                &   hig%gy(2,i,j)*u(i,j,k+2,3)+ &
                &   hig%gy(3,i,j)*u(i,j,k  ,2)+ &
                &   hig%gy(4,i,j)*u(i,j,k+1,2)+ &
                &   hig%gy(5,i,j)*u(i,j,k+2,2)+ &
                &   hig%gy(6,i,j)*u(i,j,k  ,1)+ &
                &   hig%gy(7,i,j)*u(i,j,k+1,1)+ &
                &   hig%gy(8,i,j)*u(i,j,k+2,1)
                u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
             end do
          end do
       end do
       ! Operator padding
       do k=bounds%nmin3-1,bounds%nmin3-4,-1
          do j=bounds%nmin2,bounds%nmax2
             do i=bounds%nmin1,bounds%nmax1
                fdum = hig%gy(1,i,j)*u(i,j,k+1,3)+ &
                &   hig%gy(2,i,j)*u(i,j,k+2,3)+ &
                &   hig%gy(3,i,j)*u(i,j,k  ,2)+ &
                &   hig%gy(4,i,j)*u(i,j,k+1,2)+ &
                &   hig%gy(5,i,j)*u(i,j,k+2,2)+ &
                &   hig%gy(6,i,j)*u(i,j,k  ,1)+ &
                &   hig%gy(7,i,j)*u(i,j,k+1,1)+ &
                &   hig%gy(8,i,j)*u(i,j,k+2,1)
                u(i,j,k,3) = -fdum
             end do
          end do
       end do
       !
       ! Update the back side operator pad
       ! Taper
       do k=mod%ny+1,bounds%nmax3
          fscale = float(k-mod%ny)/float(bounds%nmax3-mod%ny-1)/fscale0
          do j=bounds%nmin2,bounds%nmax2
             do i=bounds%nmin1,bounds%nmax1
                fdum = hig%gy( 9,i,j)*u(i,j,k-1,3)+ &
                &   hig%gy(10,i,j)*u(i,j,k-2,3)+ &
                &   hig%gy(11,i,j)*u(i,j,k  ,2)+ &
                &   hig%gy(12,i,j)*u(i,j,k-1,2)+ &
                &   hig%gy(13,i,j)*u(i,j,k-2,2)+ &
                &   hig%gy(14,i,j)*u(i,j,k  ,1)+ &
                &   hig%gy(15,i,j)*u(i,j,k-1,1)+ &
                &   hig%gy(16,i,j)*u(i,j,k-2,1)
                u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
             end do
          end do
       end do
       ! Operator padding
       do k=bounds%nmax3+1,bounds%nmax3+4
          do j=bounds%nmin2,bounds%nmax2
             do i=bounds%nmin1,bounds%nmax1
                fdum = hig%gy( 9,i,j)*u(i,j,k-1,3)+ &
                &   hig%gy(10,i,j)*u(i,j,k-2,3)+ &
                &   hig%gy(11,i,j)*u(i,j,k  ,2)+ &
                &   hig%gy(12,i,j)*u(i,j,k-1,2)+ &
                &   hig%gy(13,i,j)*u(i,j,k-2,2)+ &
                &   hig%gy(14,i,j)*u(i,j,k  ,1)+ &
                &   hig%gy(15,i,j)*u(i,j,k-1,1)+ &
                &   hig%gy(16,i,j)*u(i,j,k-2,1)
                u(i,j,k,3) = -fdum
             end do
          end do
       end do
    endif
    !
  end subroutine Boundary0
!
!--------------------------------------------------------------------
  subroutine Boundary1(genpar,bounds,u,mod,hig)
    !
    type(GeneralParam) :: genpar
    type(ModelSpace)   :: mod
    type(FDbounds)     :: bounds
    type(HigdonParam)  :: hig

    real :: u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &         bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)

    integer :: i, j, k, nz0,  ntaper
    real :: fscale0, fscale, fdum
    !

    ntaper = bounds%nmax2 - mod%nx
    if (ntaper.lt.genpar%lsinc/2+2) then
       nz0 = -(genpar%lsinc/2+2)+1
    else
       nz0 = 1
    endif
    fscale0 = 2.
    !
    !
    ! Update the taper zone and update the regions used for operator padding
    ! Since everything is outside the main grid, constant density applies
    !
    ! Axis 1 ----------------------
    !
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          ! Update the top
          if (genpar%surf_type.le.0) then
             do i=nz0-1,bounds%nmin1-4,-1
                fscale = float(nz0-i)/float(nz0-1-bounds%nmin1+4)/fscale0
                fdum = hig%gz(1,j,k)*u(i+1,j,k,3)+ &
                &   hig%gz(2,j,k)*u(i+2,j,k,3)+ &
                &   hig%gz(3,j,k)*u(i  ,j,k,2)+ &
                &   hig%gz(4,j,k)*u(i+1,j,k,2)+ &
                &   hig%gz(5,j,k)*u(i+2,j,k,2)+ &
                &   hig%gz(6,j,k)*u(i  ,j,k,1)+ &
                &   hig%gz(7,j,k)*u(i+1,j,k,1)+ &
                &   hig%gz(8,j,k)*u(i+2,j,k,1)
                u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
             end do
          endif
          ! Update the bottom operator pad
          do i=mod%nz+1,bounds%nmax1+4
             fscale = float(i-mod%nz)/float(bounds%nmax1+4-mod%nz-1)/fscale0
             fdum = hig%gz( 9,j,k)*u(i-1,j,k,3)+ &
             &   hig%gz(10,j,k)*u(i-2,j,k,3)+ &
             &   hig%gz(11,j,k)*u(i  ,j,k,2)+ &
             &   hig%gz(12,j,k)*u(i-1,j,k,2)+ &
             &   hig%gz(13,j,k)*u(i-2,j,k,2)+ &
             &   hig%gz(14,j,k)*u(i  ,j,k,1)+ &
             &   hig%gz(15,j,k)*u(i-1,j,k,1)+ &
             &   hig%gz(16,j,k)*u(i-2,j,k,1)
             u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
          end do
       end do
    end do
    !
    ! Axis 2 ----------------------
    !
    do k=bounds%nmin3,bounds%nmax3
       ! Update the left side operator pad
       do j=0,bounds%nmin2-4,-1
          fscale = float(1-j)/float(-(bounds%nmin2-4))/fscale0
          do i=bounds%nmin1,bounds%nmax1
             fdum = hig%gx(1,i,k)*u(i,j+1,k,3)+ &
             &   hig%gx(2,i,k)*u(i,j+2,k,3)+ &
             &   hig%gx(3,i,k)*u(i,j  ,k,2)+ &
             &   hig%gx(4,i,k)*u(i,j+1,k,2)+ &
             &   hig%gx(5,i,k)*u(i,j+2,k,2)+ &
             &   hig%gx(6,i,k)*u(i,j  ,k,1)+ &
             &   hig%gx(7,i,k)*u(i,j+1,k,1)+ &
             &   hig%gx(8,i,k)*u(i,j+2,k,1)
             u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
          end do
       end do
       ! Update the right side operator pad
       do j=mod%nx+1,bounds%nmax2+4
          fscale = float(j-mod%nx)/float(bounds%nmax2+4-mod%nx-1)/fscale0
          do i=bounds%nmin1,bounds%nmax1
             fdum = hig%gx( 9,i,k)*u(i,j-1,k,3)+ &
             &   hig%gx(10,i,k)*u(i,j-2,k,3)+ &
             &   hig%gx(11,i,k)*u(i,j  ,k,2)+ &
             &   hig%gx(12,i,k)*u(i,j-1,k,2)+ &
             &   hig%gx(13,i,k)*u(i,j-2,k,2)+ &
             &   hig%gx(14,i,k)*u(i,j  ,k,1)+ &
             &   hig%gx(15,i,k)*u(i,j-1,k,1)+ &
             &   hig%gx(16,i,k)*u(i,j-2,k,1)
             u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
          end do
       end do
    end do
    !
    ! Axis 3 ----------------------
    !
    if (genpar%nbound.gt.0) then
       ! Update the front side operator pad
       do k=0,bounds%nmin3-4,-1
          fscale = float(1-k)/float(-(bounds%nmin3-4))/fscale0
          do j=bounds%nmin2,bounds%nmax2
             do i=bounds%nmin1,bounds%nmax1
                fdum = hig%gy(1,i,j)*u(i,j,k+1,3)+ &
                &   hig%gy(2,i,j)*u(i,j,k+2,3)+ &
                &   hig%gy(3,i,j)*u(i,j,k  ,2)+ &
                &   hig%gy(4,i,j)*u(i,j,k+1,2)+ &
                &   hig%gy(5,i,j)*u(i,j,k+2,2)+ &
                &   hig%gy(6,i,j)*u(i,j,k  ,1)+ &
                &   hig%gy(7,i,j)*u(i,j,k+1,1)+ &
                &   hig%gy(8,i,j)*u(i,j,k+2,1)
                u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
             end do
          end do
       end do
       ! Update the back side operator pad
       do k=mod%ny+1,bounds%nmax3+4
          fscale = float(k-mod%ny)/float(bounds%nmax3+4-mod%ny-1)/fscale0
          do j=bounds%nmin2,bounds%nmax2
             do i=bounds%nmin1,bounds%nmax1
                fdum = hig%gy( 9,i,j)*u(i,j,k-1,3)+ &
                &   hig%gy(10,i,j)*u(i,j,k-2,3)+ &
                &   hig%gy(11,i,j)*u(i,j,k  ,2)+ &
                &   hig%gy(12,i,j)*u(i,j,k-1,2)+ &
                &   hig%gy(13,i,j)*u(i,j,k-2,2)+ &
                &   hig%gy(14,i,j)*u(i,j,k  ,1)+ &
                &   hig%gy(15,i,j)*u(i,j,k-1,1)+ &
                &   hig%gy(16,i,j)*u(i,j,k-2,1)
                u(i,j,k,3) = u(i,j,k,3)*(1.-fscale) - fdum*fscale
             end do
          end do
       end do
    endif
    !
  end subroutine Boundary1

  subroutine Boundary_set_free_surface(bounds,model,elev,u,genpar)
    type(FDbounds)    ::               bounds
    type(ModelSpace)  ::                      model
    type(ModelSpace_elevation) ::                   elev
    type(GeneralParam)::                                   genpar
    real              ::                                 u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                                      bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)

    integer:: i,j,k,l

    real, allocatable :: sinc(:),uinterp(:)
    real :: fdum
    
    allocate(sinc(genpar%lsinc),uinterp(genpar%lsinc))

    if (genpar%surf_type.ne.0) then
       do k=bounds%nmin3,bounds%nmax3
          do j=bounds%nmin2,bounds%nmax2
             ! Set the pressure to zero above the free surface
             do i=bounds%nmin1,elev%ielev_z(j,k)-1
                u(i,j,k,2) = 0.
             end do
             if (elev%delev_z(j,k).eq.0.) then
                do i=0,2
                   u(elev%ielev_z(j,k)-2-i,j,k,2) = -u(elev%ielev_z(j,k)+i,j,k,2)
                end do
             else
                ! Interpolate between grid
                call mksinc(sinc,genpar%lsinc,elev%delev_z(j,k)/model%dz)
                do i=-genpar%lsinc/2,genpar%lsinc/2
                   fdum = 0.
                   do l=-genpar%lsinc/2,genpar%lsinc/2
                      fdum = fdum + u(elev%ielev_z(j,k)+l+i,j,k,2)*sinc(genpar%lsinc/2+l+1)
                   end do
                   uinterp(genpar%lsinc/2+i+1) = fdum
                end do
                ! Add inverse mirror image
                do i=2,genpar%lsinc/2+1
                   do l=-genpar%lsinc/2,i-1
                      u(elev%ielev_z(j,k)-i+l,j,k,2) = u(elev%ielev_z(j,k)-i+l,j,k,2)- &
                      &  uinterp(genpar%lsinc/2+i-1)*sinc(genpar%lsinc/2+l+1)
                   end do
                end do
             endif
          end do
       end do
    endif

    deallocate(sinc,uinterp)

  end subroutine Boundary_set_free_surface

end module Boundary_mod
