! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module FD_derivatives

  use FD_types
  use ModelSpace_types
  use GeneralParam_types
  use omp_lib

  implicit none

  type(ScaledFDcoefs), pointer, private  :: coefs
  real, dimension(:,:,:), allocatable, private :: s3dxx,s3dzz,s3dyy,delp3d
  real, dimension(:,:,:), allocatable, private :: s2dxx,s2dzz,delp2d
  
  !$omp threadprivate(s2dxx,s2dzz,delp2d)

contains

  subroutine FD_derivatives_coef_init(coef_init)
    type(ScaledFDcoefs), target :: coef_init
    coefs => coef_init
  end subroutine FD_derivatives_coef_init

  subroutine FD_2nd_3D_derivatives_scalar_forward(genpar,bounds,u,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j,k

    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz= coefs%c0z*(u(i  ,j,k,2)) +&
             &      coefs%c1z*(u(i+1,j,k,2)+u(i-1,j,k,2))+ &
             &      coefs%c2z*(u(i+2,j,k,2)+u(i-2,j,k,2))+ &
             &      coefs%c3z*(u(i+3,j,k,2)+u(i-3,j,k,2))+ &
             &      coefs%c4z*(u(i+4,j,k,2)+u(i-4,j,k,2))
             tmpxx= coefs%c0x*(u(i,j  ,k,2)) +&
             &      coefs%c1x*(u(i,j+1,k,2)+u(i,j-1,k,2))+ &
             &      coefs%c2x*(u(i,j+2,k,2)+u(i,j-2,k,2))+ &
             &      coefs%c3x*(u(i,j+3,k,2)+u(i,j-3,k,2))+ &
             &      coefs%c4x*(u(i,j+4,k,2)+u(i,j-4,k,2))
             tmpyy= coefs%c0y*(u(i,j,k  ,2)) +&
             &      coefs%c1y*(u(i,j,k+1,2)+u(i,j,k-1,2))+ &
             &      coefs%c2y*(u(i,j,k+2,2)+u(i,j,k-2,2))+ &
             &      coefs%c3y*(u(i,j,k+3,2)+u(i,j,k-3,2))+ &
             &      coefs%c4y*(u(i,j,k+4,2)+u(i,j,k-4,2))
             u(i,j,k,3)=mod%vel(i,j,k)**2*(tmpxx+tmpzz+tmpyy)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine FD_2nd_3D_derivatives_scalar_forward
  
  subroutine FD_2nd_3D_derivatives_scalar_forward_grid_noomp(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                        genpar
    type(ModelSpace)     ::                                            mod
    type(FDbounds)       ::                               bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j,k

    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz= coefs%c0z*(u2(i  ,j,k)) +&
             &      coefs%c1z*(u2(i+1,j,k)+u2(i-1,j,k))+ &
             &      coefs%c2z*(u2(i+2,j,k)+u2(i-2,j,k))+ &
             &      coefs%c3z*(u2(i+3,j,k)+u2(i-3,j,k))+ &
             &      coefs%c4z*(u2(i+4,j,k)+u2(i-4,j,k))
             tmpxx= coefs%c0x*(u2(i,j  ,k)) +&
             &      coefs%c1x*(u2(i,j+1,k)+u2(i,j-1,k))+ &
             &      coefs%c2x*(u2(i,j+2,k)+u2(i,j-2,k))+ &
             &      coefs%c3x*(u2(i,j+3,k)+u2(i,j-3,k))+ &
             &      coefs%c4x*(u2(i,j+4,k)+u2(i,j-4,k))
             tmpyy= coefs%c0y*(u2(i,j,k  )) +&
             &      coefs%c1y*(u2(i,j,k+1)+u2(i,j,k-1))+ &
             &      coefs%c2y*(u2(i,j,k+2)+u2(i,j,k-2))+ &
             &      coefs%c3y*(u2(i,j,k+3)+u2(i,j,k-3))+ &
             &      coefs%c4y*(u2(i,j,k+4)+u2(i,j,k-4))
             u3(i,j,k)=mod%vel(i,j,k)**2*(tmpxx+tmpzz+tmpyy)
          end do
       end do
    end do
    
  end subroutine FD_2nd_3D_derivatives_scalar_forward_grid_noomp
  
  subroutine FD_2nd_3D_derivatives_scalar_forward_grid(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                        genpar
    type(ModelSpace)     ::                                            mod
    type(FDbounds)       ::                               bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j,k

    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz= coefs%c0z*(u2(i  ,j,k)) +&
             &      coefs%c1z*(u2(i+1,j,k)+u2(i-1,j,k))+ &
             &      coefs%c2z*(u2(i+2,j,k)+u2(i-2,j,k))+ &
             &      coefs%c3z*(u2(i+3,j,k)+u2(i-3,j,k))+ &
             &      coefs%c4z*(u2(i+4,j,k)+u2(i-4,j,k))
             tmpxx= coefs%c0x*(u2(i,j  ,k)) +&
             &      coefs%c1x*(u2(i,j+1,k)+u2(i,j-1,k))+ &
             &      coefs%c2x*(u2(i,j+2,k)+u2(i,j-2,k))+ &
             &      coefs%c3x*(u2(i,j+3,k)+u2(i,j-3,k))+ &
             &      coefs%c4x*(u2(i,j+4,k)+u2(i,j-4,k))
             tmpyy= coefs%c0y*(u2(i,j,k  )) +&
             &      coefs%c1y*(u2(i,j,k+1)+u2(i,j,k-1))+ &
             &      coefs%c2y*(u2(i,j,k+2)+u2(i,j,k-2))+ &
             &      coefs%c3y*(u2(i,j,k+3)+u2(i,j,k-3))+ &
             &      coefs%c4y*(u2(i,j,k+4)+u2(i,j,k-4))
             u3(i,j,k)=mod%vel(i,j,k)**2*(tmpxx+tmpzz+tmpyy)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine FD_2nd_3D_derivatives_scalar_forward_grid
  
  subroutine FD_3D_derivatives_acoustic_forward(genpar,bounds,u,mod)
    type(GeneralParam)   ::                     genpar
    type(ModelSpace)     ::                                     mod
    type(FDbounds)       ::                            bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx, tmpyy
    real, dimension(:,:,:), allocatable :: sxx,szz,syy,delp
    real                 :: dxi,dyi,dzi
    integer              :: i,j,k

    allocate(sxx(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))
    allocate(syy(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))
    allocate(szz(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))
    allocate(delp(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))  
    
    delp=1./mod%rho2
    dzi =1./genpar%delta(1)
    dxi =1./genpar%delta(2)
    dyi =1./genpar%delta(3)

    szz=0.
    syy=0.
    sxx=0.

    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1   
             szz(i,j,k)=           (coefs%c1z*(u(i+1,j,k,2)-u(i,j,k,2))+ &
             &                      coefs%c2z*(u(i+2,j,k,2)-u(i-1,j,k,2))+ &
             &                      coefs%c3z*(u(i+3,j,k,2)-u(i-2,j,k,2))+ &
             &                      coefs%c4z*(u(i+4,j,k,2)-u(i-3,j,k,2)))*delp(i,j,k)
             sxx(i,j,k)=           (coefs%c1x*(u(i,j+1,k,2)-u(i,j,k,2))+ &
             &                      coefs%c2x*(u(i,j+2,k,2)-u(i,j-1,k,2))+ &
             &                      coefs%c3x*(u(i,j+3,k,2)-u(i,j-2,k,2))+ &
             &                      coefs%c4x*(u(i,j+4,k,2)-u(i,j-3,k,2)))*delp(i,j,k)
             syy(i,j,k)=           (coefs%c1y*(u(i,j,k+1,2)-u(i,j,k,2))+ &
             &                      coefs%c2y*(u(i,j,k+2,2)-u(i,j,k-1,2))+ &
             &                      coefs%c3y*(u(i,j,k+3,2)-u(i,j,k-2,2))+ &
             &                      coefs%c4y*(u(i,j,k+4,2)-u(i,j,k-3,2)))*delp(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
     
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          ! Gradient of the p-wave potential in the top operator padding 
          do i=bounds%nmin1-1,bounds%nmin1-3,-1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(bounds%nmin1,j,k)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(bounds%nmin1,j,k)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(bounds%nmin1,j,k)
          end do
          ! Gradient of the p-wave potential in the bottom operator padding
          do i=bounds%nmax1+1,bounds%nmax1+3
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(bounds%nmax1,j,k)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(bounds%nmax1,j,k)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(bounds%nmax1,j,k)
          end do
       end do
       ! Gradient of the p-wave potential in the left side operator padding 
       do j=bounds%nmin2-1,bounds%nmin2-3,-1
          do i=bounds%nmin1,bounds%nmax1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(i,bounds%nmin2,k)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(i,bounds%nmin2,k)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(i,bounds%nmin2,k)
          end do
       end do
       ! Gradient of the p-wave potential in the right side operator padding 
       do j=bounds%nmax2+1,bounds%nmax2+3
          do i=bounds%nmin1,bounds%nmax1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(i,bounds%nmax2,k)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(i,bounds%nmax2,k)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(i,bounds%nmax2,k)
          end do
       end do
    end do
    
    ! Gradient of the p-wave potential in the front side operator padding 
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3-1,bounds%nmin3-3,-1
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmin3)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmin3)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(i,j,bounds%nmin3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Gradient of the p-wave potential in the back side operator padding 
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmax3+1,bounds%nmax3+3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmax3)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmax3)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(i,j,bounds%nmax3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz=         coefs%c1z*(szz(i  ,j,k)-szz(i-1,j,k))+ &
             &              coefs%c2z*(szz(i+1,j,k)-szz(i-2,j,k))+ &
             &              coefs%c3z*(szz(i+2,j,k)-szz(i-3,j,k))+ &
             &              coefs%c4z*(szz(i+3,j,k)-szz(i-4,j,k))
             tmpxx=         coefs%c1x*(sxx(i  ,j,k)-sxx(i,j-1,k))+ &
             &              coefs%c2x*(sxx(i,j+1,k)-sxx(i,j-2,k))+ &
             &              coefs%c3x*(sxx(i,j+2,k)-sxx(i,j-3,k))+ &
             &              coefs%c4x*(sxx(i,j+3,k)-sxx(i,j-4,k))
             tmpyy=         coefs%c1y*(syy(i,j,k  )-syy(i,j,k-1))+ &
             &              coefs%c2y*(syy(i,j,k+1)-syy(i,j,k-2))+ &
             &              coefs%c3y*(syy(i,j,k+2)-syy(i,j,k-3))+ &
             &              coefs%c4y*(syy(i,j,k+3)-syy(i,j,k-4))
             u(i,j,k,3)=mod%vel(i,j,k)**2* &
             &              mod%rho(i,j,k)*(tmpxx+tmpyy+tmpzz)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    deallocate(sxx,syy,szz,delp)

  end subroutine FD_3D_derivatives_acoustic_forward
  
  subroutine FD_3D_derivatives_acoustic_forward_grid(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                     genpar
    type(ModelSpace)     ::                                     mod
    type(FDbounds)       ::                            bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx, tmpyy
    real                 :: dxi,dyi,dzi
    integer              :: i,j,k
  
    dzi =1./genpar%delta(1)
    dxi =1./genpar%delta(2)
    dyi =1./genpar%delta(3)

    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1   
             s3dzz(i,j,k)=           (coefs%c1z*(u2(i+1,j,k)-u2(i,j,k))+ &
             &                      coefs%c2z*(u2(i+2,j,k)-u2(i-1,j,k))+ &
             &                      coefs%c3z*(u2(i+3,j,k)-u2(i-2,j,k))+ &
             &                      coefs%c4z*(u2(i+4,j,k)-u2(i-3,j,k)))*delp3d(i,j,k)
             s3dxx(i,j,k)=           (coefs%c1x*(u2(i,j+1,k)-u2(i,j,k))+ &
             &                      coefs%c2x*(u2(i,j+2,k)-u2(i,j-1,k))+ &
             &                      coefs%c3x*(u2(i,j+3,k)-u2(i,j-2,k))+ &
             &                      coefs%c4x*(u2(i,j+4,k)-u2(i,j-3,k)))*delp3d(i,j,k)
             s3dyy(i,j,k)=           (coefs%c1y*(u2(i,j,k+1)-u2(i,j,k))+ &
             &                      coefs%c2y*(u2(i,j,k+2)-u2(i,j,k-1))+ &
             &                      coefs%c3y*(u2(i,j,k+3)-u2(i,j,k-2))+ &
             &                      coefs%c4y*(u2(i,j,k+4)-u2(i,j,k-3)))*delp3d(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
     
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          ! Gradient of the p-wave potential in the top operator padding 
          do i=bounds%nmin1-1,bounds%nmin1-3,-1
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*delp3d(bounds%nmin1,j,k)
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*delp3d(bounds%nmin1,j,k)
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*delp3d(bounds%nmin1,j,k)
          end do
          ! Gradient of the p-wave potential in the bottom operator padding
          do i=bounds%nmax1+1,bounds%nmax1+3
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*delp3d(bounds%nmax1,j,k)
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*delp3d(bounds%nmax1,j,k)
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*delp3d(bounds%nmax1,j,k)
          end do
       end do
       ! Gradient of the p-wave potential in the left side operator padding 
       do j=bounds%nmin2-1,bounds%nmin2-3,-1
          do i=bounds%nmin1,bounds%nmax1
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*delp3d(i,bounds%nmin2,k)
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*delp3d(i,bounds%nmin2,k)
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*delp3d(i,bounds%nmin2,k)
          end do
       end do
       ! Gradient of the p-wave potential in the right side operator padding 
       do j=bounds%nmax2+1,bounds%nmax2+3
          do i=bounds%nmin1,bounds%nmax1
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*delp3d(i,bounds%nmax2,k)
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*delp3d(i,bounds%nmax2,k)
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*delp3d(i,bounds%nmax2,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Gradient of the p-wave potential in the front side operator padding 
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3-1,bounds%nmin3-3,-1
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*delp3d(i,j,bounds%nmin3)
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*delp3d(i,j,bounds%nmin3)
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*delp3d(i,j,bounds%nmin3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    ! Gradient of the p-wave potential in the back side operator padding 
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmax3+1,bounds%nmax3+3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*delp3d(i,j,bounds%nmax3)
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*delp3d(i,j,bounds%nmax3)
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*delp3d(i,j,bounds%nmax3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz=         coefs%c1z*(s3dzz(i  ,j,k)-s3dzz(i-1,j,k))+ &
             &              coefs%c2z*(s3dzz(i+1,j,k)-s3dzz(i-2,j,k))+ &
             &              coefs%c3z*(s3dzz(i+2,j,k)-s3dzz(i-3,j,k))+ &
             &              coefs%c4z*(s3dzz(i+3,j,k)-s3dzz(i-4,j,k))
             tmpxx=         coefs%c1x*(s3dxx(i  ,j,k)-s3dxx(i,j-1,k))+ &
             &              coefs%c2x*(s3dxx(i,j+1,k)-s3dxx(i,j-2,k))+ &
             &              coefs%c3x*(s3dxx(i,j+2,k)-s3dxx(i,j-3,k))+ &
             &              coefs%c4x*(s3dxx(i,j+3,k)-s3dxx(i,j-4,k))
             tmpyy=         coefs%c1y*(s3dyy(i,j,k  )-s3dyy(i,j,k-1))+ &
             &              coefs%c2y*(s3dyy(i,j,k+1)-s3dyy(i,j,k-2))+ &
             &              coefs%c3y*(s3dyy(i,j,k+2)-s3dyy(i,j,k-3))+ &
             &              coefs%c4y*(s3dyy(i,j,k+3)-s3dyy(i,j,k-4))
             u3(i,j,k)=mod%vel(i,j,k)**2* &
             &              mod%rho(i,j,k)*(tmpxx+tmpyy+tmpzz)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine FD_3D_derivatives_acoustic_forward_grid
  

  subroutine allocate_sxx3d_syy3d_szz3d_delp3d(bounds,mod)    
    type(FDbounds)       ::                    bounds
    type(ModelSpace)     ::                           mod

    allocate(s3dxx(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))
    allocate(s3dyy(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))
    allocate(s3dzz(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4)) 
    allocate(delp3d(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-4:bounds%nmax3+4))     
    delp3d=1./mod%rho2

  end subroutine allocate_sxx3d_syy3d_szz3d_delp3d

  subroutine deallocate_sxx3d_syy3d_szz3d_delp3d()
    deallocate(s3dxx,s3dyy,s3dzz,delp3d)
  end subroutine deallocate_sxx3d_syy3d_szz3d_delp3d

  subroutine FD_3D_derivatives_acoustic_forward_grid_opt(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                     genpar
    type(ModelSpace)     ::                                     mod
    type(FDbounds)       ::                            bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx, tmpyy
    real                 :: dxi,dyi,dzi,alpha
    integer              :: i,j,k,iz,iy,ix

    dzi =1./genpar%delta(1)
    dxi =1./genpar%delta(2)
    dyi =1./genpar%delta(3)

    ! First compute derivatives everywhere with first order stencil
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3-4,bounds%nmax3+4
       do j=bounds%nmin2-4,bounds%nmax2+4
          do i=bounds%nmin1-4,bounds%nmax1+4  
             s3dzz(i,j,k)=0.
             s3dyy(i,j,k)=0.
             s3dxx(i,j,k)=0.
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! First compute derivatives everywhere with first order stencil
    !$OMP PARALLEL DO PRIVATE(k,j,i,ix,iy,iz,alpha)
    do k=bounds%nmin3-3,bounds%nmax3+3
       do j=bounds%nmin2-3,bounds%nmax2+3
          do i=bounds%nmin1-3,bounds%nmax1+3   
             iz= min(max(i,bounds%nmin1),bounds%nmax1)
             ix= min(max(j,bounds%nmin2),bounds%nmax2)
             iy= min(max(k,bounds%nmin3),bounds%nmax3)
             
             alpha=delp3d(iz,ix,iy)
             s3dzz(i,j,k)=dzi*(u2(i+1,j,k)-u2(i,j,k))*alpha
             s3dxx(i,j,k)=dxi*(u2(i,j+1,k)-u2(i,j,k))*alpha
             s3dyy(i,j,k)=dyi*(u2(i,j,k+1)-u2(i,j,k))*alpha
          end do
       end do
    end do
    !$OMP END PARALLEL DO
 
    ! First compute derivatives everywhere with first order stencil
    !$OMP PARALLEL DO PRIVATE(ix,iy,iz,alpha)
    do iy=bounds%nmin3,bounds%nmax3
       do ix=bounds%nmin2,bounds%nmax2
          do iz=bounds%nmin1,bounds%nmax1
             alpha=delp3d(iz,ix,iy)
             s3dzz(iz,ix,iy)=      (coefs%c1z*(u2(iz+1,ix  ,iy  )-u2(iz  ,ix  ,iy  ))+ &
             &                      coefs%c2z*(u2(iz+2,ix  ,iy  )-u2(iz-1,ix  ,iy  ))+ &
             &                      coefs%c3z*(u2(iz+3,ix  ,iy  )-u2(iz-2,ix  ,iy  ))+ &
             &                      coefs%c4z*(u2(iz+4,ix  ,iy  )-u2(iz-3,ix  ,iy  )))*alpha
             s3dxx(iz,ix,iy)=      (coefs%c1x*(u2(iz  ,ix+1,iy  )-u2(iz  ,ix  ,iy  ))+ &
             &                      coefs%c2x*(u2(iz  ,ix+2,iy  )-u2(iz  ,ix-1,iy  ))+ &
             &                      coefs%c3x*(u2(iz  ,ix+3,iy  )-u2(iz  ,ix-2,iy  ))+ &
             &                      coefs%c4x*(u2(iz  ,ix+4,iy  )-u2(iz  ,ix-3,iy  )))*alpha
             s3dyy(iz,ix,iy)=      (coefs%c1y*(u2(iz  ,ix  ,iy+1)-u2(iz  ,ix  ,iy  ))+ &
             &                      coefs%c2y*(u2(iz  ,ix  ,iy+2)-u2(iz  ,ix  ,iy-1))+ &
             &                      coefs%c3y*(u2(iz  ,ix  ,iy+3)-u2(iz  ,ix  ,iy-2))+ &
             &                      coefs%c4y*(u2(iz  ,ix  ,iy+4)-u2(iz  ,ix  ,iy-3)))*alpha
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz=         coefs%c1z*(s3dzz(i  ,j,k)-s3dzz(i-1,j,k))+ &
             &              coefs%c2z*(s3dzz(i+1,j,k)-s3dzz(i-2,j,k))+ &
             &              coefs%c3z*(s3dzz(i+2,j,k)-s3dzz(i-3,j,k))+ &
             &              coefs%c4z*(s3dzz(i+3,j,k)-s3dzz(i-4,j,k))
             tmpxx=         coefs%c1x*(s3dxx(i  ,j,k)-s3dxx(i,j-1,k))+ &
             &              coefs%c2x*(s3dxx(i,j+1,k)-s3dxx(i,j-2,k))+ &
             &              coefs%c3x*(s3dxx(i,j+2,k)-s3dxx(i,j-3,k))+ &
             &              coefs%c4x*(s3dxx(i,j+3,k)-s3dxx(i,j-4,k))
             tmpyy=         coefs%c1y*(s3dyy(i,j,k  )-s3dyy(i,j,k-1))+ &
             &              coefs%c2y*(s3dyy(i,j,k+1)-s3dyy(i,j,k-2))+ &
             &              coefs%c3y*(s3dyy(i,j,k+2)-s3dyy(i,j,k-3))+ &
             &              coefs%c4y*(s3dyy(i,j,k+3)-s3dyy(i,j,k-4))
             u3(i,j,k)=mod%vel(i,j,k)**2* &
             &              mod%rho(i,j,k)*(tmpxx+tmpyy+tmpzz)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine FD_3D_derivatives_acoustic_forward_grid_opt
  
  subroutine FD_2D_derivatives_acoustic_forward(genpar,bounds,u,mod)
    type(GeneralParam)   ::                     genpar
    type(ModelSpace)     ::                                     mod
    type(FDbounds)       ::                            bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx
    real, dimension(:,:,:), allocatable :: sxx,szz,delp
    real                 :: dxi,dzi
    integer              :: i,j

    allocate(sxx(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3:bounds%nmax3))
    allocate(szz(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3:bounds%nmax3))  
    allocate(delp(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3:bounds%nmax3))
    
    delp=1./mod%rho2
    dzi =1./genpar%delta(1)
    dxi =1./genpar%delta(2)

    sxx=0.
    szz=0. 

    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1   
          szz(i,j,1)=           (coefs%c1z*(u(i+1,j,1,2)-u(i  ,j,1,2))+ &
          &                      coefs%c2z*(u(i+2,j,1,2)-u(i-1,j,1,2))+ &
          &                      coefs%c3z*(u(i+3,j,1,2)-u(i-2,j,1,2))+ &
          &                      coefs%c4z*(u(i+4,j,1,2)-u(i-3,j,1,2)))*delp(i,j,1)
          sxx(i,j,1)=           (coefs%c1x*(u(i,j+1,1,2)-u(i  ,j,1,2))+ &
          &                      coefs%c2x*(u(i,j+2,1,2)-u(i,j-1,1,2))+ &
          &                      coefs%c3x*(u(i,j+3,1,2)-u(i,j-2,1,2))+ &
          &                      coefs%c4x*(u(i,j+4,1,2)-u(i,j-3,1,2)))*delp(i,j,1)
       end do
    end do
     
    do j=bounds%nmin2,bounds%nmax2
       ! Gradient of the p-wave potential in the top operator padding 
       do i=bounds%nmin1-1,bounds%nmin1-3,-1
          szz(i,j,1)=dzi*(u(i+1,j,1,2)-u(i,j,1,2))*delp(bounds%nmin1,j,1)
          sxx(i,j,1)=dxi*(u(i,j+1,1,2)-u(i,j,1,2))*delp(bounds%nmin1,j,1)
       end do
       ! Gradient of the p-wave potential in the bottom operator padding
       do i=bounds%nmax1+1,bounds%nmax1+3
          szz(i,j,1)=dzi*(u(i+1,j,1,2)-u(i,j,1,2))*delp(bounds%nmax1,j,1)
          sxx(i,j,1)=dxi*(u(i,j+1,1,2)-u(i,j,1,2))*delp(bounds%nmax1,j,1)
       end do
    end do
    ! Gradient of the p-wave potential in the left side operator padding 
    do j=bounds%nmin2-1,bounds%nmin2-3,-1
       do i=bounds%nmin1,bounds%nmax1
          szz(i,j,1)=dzi*(u(i+1,j,1,2)-u(i,j,1,2))*delp(i,bounds%nmin2,1)
          sxx(i,j,1)=dxi*(u(i,j+1,1,2)-u(i,j,1,2))*delp(i,bounds%nmin2,1)
       end do
    end do
    
    ! Gradient of the p-wave potential in the right side operator padding 
    do j=bounds%nmax2+1,bounds%nmax2+3
       do i=bounds%nmin1,bounds%nmax1
          szz(i,j,1)=dzi*(u(i+1,j,1,2)-u(i,j,1,2))*delp(i,bounds%nmax2,1)
          sxx(i,j,1)=dxi*(u(i,j+1,1,2)-u(i,j,1,2))*delp(i,bounds%nmax2,1)
       end do
    end do 
    
    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz=         coefs%c1z*(szz(i,j,1)-szz(i-1,j,1))+ &
          &              coefs%c2z*(szz(i+1,j,1)-szz(i-2,j,1))+ &
          &              coefs%c3z*(szz(i+2,j,1)-szz(i-3,j,1))+ &
          &              coefs%c4z*(szz(i+3,j,1)-szz(i-4,j,1))
          tmpxx=         coefs%c1x*(sxx(i,j,1)-sxx(i,j-1,1))+ &
          &              coefs%c2x*(sxx(i,j+1,1)-sxx(i,j-2,1))+ &
          &              coefs%c3x*(sxx(i,j+2,1)-sxx(i,j-3,1))+ &
          &              coefs%c4x*(sxx(i,j+3,1)-sxx(i,j-4,1))
          u(i,j,1,3)=mod%vel(i,j,1)*mod%vel(i,j,1)* &
          &              mod%rho(i,j,1)*(tmpxx+tmpzz)
       end do
    end do
    deallocate(sxx,szz,delp)

  end subroutine FD_2D_derivatives_acoustic_forward
  
  subroutine allocate_sxx2d_szz2d_delp2d(bounds,mod)    
    type(FDbounds)       ::              bounds
    type(ModelSpace)     ::                     mod

    allocate(s2dxx(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3:bounds%nmax3))
    allocate(s2dzz(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3:bounds%nmax3)) 
    allocate(delp2d(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3:bounds%nmax3))     
    delp2d=1./mod%rho2
  end subroutine allocate_sxx2d_szz2d_delp2d

  subroutine deallocate_sxx2d_szz2d_delp2d()
    deallocate(s2dxx,s2dzz,delp2d)
  end subroutine deallocate_sxx2d_szz2d_delp2d

  subroutine FD_2D_derivatives_acoustic_forward_grid(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                     genpar
    type(ModelSpace)     ::                                     mod
    type(FDbounds)       ::                            bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx
    real                 :: dxi,dzi
    integer              :: i,j
    
    dzi =1./genpar%delta(1)
    dxi =1./genpar%delta(2)

    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1   
          s2dzz(i,j,1)=           (coefs%c1z*(u2(i+1,j,1)-u2(i  ,j,1))+ &
          &                      coefs%c2z*(u2(i+2,j,1)-u2(i-1,j,1))+ &
          &                      coefs%c3z*(u2(i+3,j,1)-u2(i-2,j,1))+ &
          &                      coefs%c4z*(u2(i+4,j,1)-u2(i-3,j,1)))*delp2d(i,j,1)
          s2dxx(i,j,1)=           (coefs%c1x*(u2(i,j+1,1)-u2(i  ,j,1))+ &
          &                      coefs%c2x*(u2(i,j+2,1)-u2(i,j-1,1))+ &
          &                      coefs%c3x*(u2(i,j+3,1)-u2(i,j-2,1))+ &
          &                      coefs%c4x*(u2(i,j+4,1)-u2(i,j-3,1)))*delp2d(i,j,1)
       end do
    end do
     
    do j=bounds%nmin2,bounds%nmax2
       ! Gradient of the p-wave potential in the top operator padding 
       do i=bounds%nmin1-1,bounds%nmin1-3,-1
          s2dzz(i,j,1)=dzi*(u2(i+1,j,1)-u2(i,j,1))*delp2d(bounds%nmin1,j,1)
          s2dxx(i,j,1)=dxi*(u2(i,j+1,1)-u2(i,j,1))*delp2d(bounds%nmin1,j,1)
       end do
       ! Gradient of the p-wave potential in the bottom operator padding
       do i=bounds%nmax1+1,bounds%nmax1+3
          s2dzz(i,j,1)=dzi*(u2(i+1,j,1)-u2(i,j,1))*delp2d(bounds%nmax1,j,1)
          s2dxx(i,j,1)=dxi*(u2(i,j+1,1)-u2(i,j,1))*delp2d(bounds%nmax1,j,1)
       end do
    end do
    ! Gradient of the p-wave potential in the left side operator padding 
    do j=bounds%nmin2-1,bounds%nmin2-3,-1
       do i=bounds%nmin1,bounds%nmax1
          s2dzz(i,j,1)=dzi*(u2(i+1,j,1)-u2(i,j,1))*delp2d(i,bounds%nmin2,1)
          s2dxx(i,j,1)=dxi*(u2(i,j+1,1)-u2(i,j,1))*delp2d(i,bounds%nmin2,1)
       end do
    end do
    
    ! Gradient of the p-wave potential in the right side operator padding 
    do j=bounds%nmax2+1,bounds%nmax2+3
       do i=bounds%nmin1,bounds%nmax1
          s2dzz(i,j,1)=dzi*(u2(i+1,j,1)-u2(i,j,1))*delp2d(i,bounds%nmax2,1)
          s2dxx(i,j,1)=dxi*(u2(i,j+1,1)-u2(i,j,1))*delp2d(i,bounds%nmax2,1)
       end do
    end do 
    
    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz=         coefs%c1z*(s2dzz(i,j,1)-s2dzz(i-1,j,1))+ &
          &              coefs%c2z*(s2dzz(i+1,j,1)-s2dzz(i-2,j,1))+ &
          &              coefs%c3z*(s2dzz(i+2,j,1)-s2dzz(i-3,j,1))+ &
          &              coefs%c4z*(s2dzz(i+3,j,1)-s2dzz(i-4,j,1))
          tmpxx=         coefs%c1x*(s2dxx(i,j,1)-s2dxx(i,j-1,1))+ &
          &              coefs%c2x*(s2dxx(i,j+1,1)-s2dxx(i,j-2,1))+ &
          &              coefs%c3x*(s2dxx(i,j+2,1)-s2dxx(i,j-3,1))+ &
          &              coefs%c4x*(s2dxx(i,j+3,1)-s2dxx(i,j-4,1))
          u3(i,j,1)=mod%vel(i,j,1)*mod%vel(i,j,1)* &
          &              mod%rho(i,j,1)*(tmpxx+tmpzz)
       end do
    end do

  end subroutine FD_2D_derivatives_acoustic_forward_grid
  
  subroutine FD_2nd_2D_derivatives_scalar_forward(genpar,bounds,u,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j

    !!!$OMP PARALLEL DO PRIVATE(j,i,tmpxx,tmpzz)
    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz= coefs%c0z*(u(i  ,j,1,2)) +&
          &      coefs%c1z*(u(i+1,j,1,2)+u(i-1,j,1,2))+ &
          &      coefs%c2z*(u(i+2,j,1,2)+u(i-2,j,1,2))+ &
          &      coefs%c3z*(u(i+3,j,1,2)+u(i-3,j,1,2))+ &
          &      coefs%c4z*(u(i+4,j,1,2)+u(i-4,j,1,2))
          tmpxx= coefs%c0x*(u(i,j  ,1,2)) +&
          &      coefs%c1x*(u(i,j+1,1,2)+u(i,j-1,1,2))+ &
          &      coefs%c2x*(u(i,j+2,1,2)+u(i,j-2,1,2))+ &
          &      coefs%c3x*(u(i,j+3,1,2)+u(i,j-3,1,2))+ &
          &      coefs%c4x*(u(i,j+4,1,2)+u(i,j-4,1,2))
          u(i,j,1,3)=mod%vel(i,j,1)**2*(tmpxx+tmpzz)
       end do
    end do
    !!!$OMP END PARALLEL DO
    
  end subroutine FD_2nd_2D_derivatives_scalar_forward
  
  subroutine FD_2nd_2D_derivatives_scalar_forward_grid_noomp(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    type(USpace)         :: grid
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j

    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz= coefs%c0z*(u2(i  ,j,1)) +&
          &      coefs%c1z*(u2(i+1,j,1)+u2(i-1,j,1))+ &
          &      coefs%c2z*(u2(i+2,j,1)+u2(i-2,j,1))+ &
          &      coefs%c3z*(u2(i+3,j,1)+u2(i-3,j,1))+ &
          &      coefs%c4z*(u2(i+4,j,1)+u2(i-4,j,1))
          tmpxx= coefs%c0x*(u2(i,j  ,1)) +&
          &      coefs%c1x*(u2(i,j+1,1)+u2(i,j-1,1))+ &
          &      coefs%c2x*(u2(i,j+2,1)+u2(i,j-2,1))+ &
          &      coefs%c3x*(u2(i,j+3,1)+u2(i,j-3,1))+ &
          &      coefs%c4x*(u2(i,j+4,1)+u2(i,j-4,1))
          u3(i,j,1)=mod%vel(i,j,1)**2*(tmpxx+tmpzz)
       end do
    end do
    
  end subroutine FD_2nd_2D_derivatives_scalar_forward_grid_noomp
  
  subroutine FD_2nd_2D_derivatives_scalar_forward_grid(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    type(USpace)         :: grid
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j

    !$OMP PARALLEL DO PRIVATE(j,i,tmpxx,tmpzz)
    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz= coefs%c0z*(u2(i  ,j,1)) +&
          &      coefs%c1z*(u2(i+1,j,1)+u2(i-1,j,1))+ &
          &      coefs%c2z*(u2(i+2,j,1)+u2(i-2,j,1))+ &
          &      coefs%c3z*(u2(i+3,j,1)+u2(i-3,j,1))+ &
          &      coefs%c4z*(u2(i+4,j,1)+u2(i-4,j,1))
          tmpxx= coefs%c0x*(u2(i,j  ,1)) +&
          &      coefs%c1x*(u2(i,j+1,1)+u2(i,j-1,1))+ &
          &      coefs%c2x*(u2(i,j+2,1)+u2(i,j-2,1))+ &
          &      coefs%c3x*(u2(i,j+3,1)+u2(i,j-3,1))+ &
          &      coefs%c4x*(u2(i,j+4,1)+u2(i,j-4,1))
          u3(i,j,1)=mod%vel(i,j,1)**2*(tmpxx+tmpzz)
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine FD_2nd_2D_derivatives_scalar_forward_grid
  
  subroutine FD_2nd_3D_derivatives_scalar_adjoint(genpar,bounds,u,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j,k

    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz= coefs%c0z*(mod%vel(i  ,j,k)**2*u(i  ,j,k,2)) +&
             &      coefs%c1z*(mod%vel(i+1,j,k)**2*u(i+1,j,k,2)+mod%vel(i-1,j,k)**2*u(i-1,j,k,2))+ &
             &      coefs%c2z*(mod%vel(i+2,j,k)**2*u(i+2,j,k,2)+mod%vel(i-2,j,k)**2*u(i-2,j,k,2))+ &
             &      coefs%c3z*(mod%vel(i+3,j,k)**2*u(i+3,j,k,2)+mod%vel(i-3,j,k)**2*u(i-3,j,k,2))+ &
             &      coefs%c4z*(mod%vel(i+4,j,k)**2*u(i+4,j,k,2)+mod%vel(i-4,j,k)**2*u(i-4,j,k,2))
             tmpxx= coefs%c0x*(mod%vel(i  ,j,k)**2*u(i,j  ,k,2)) +&
             &      coefs%c1x*(mod%vel(i,j+1,k)**2*u(i,j+1,k,2)+mod%vel(i,j-1,k)**2*u(i,j-1,k,2))+ &
             &      coefs%c2x*(mod%vel(i,j+2,k)**2*u(i,j+2,k,2)+mod%vel(i,j-2,k)**2*u(i,j-2,k,2))+ &
             &      coefs%c3x*(mod%vel(i,j+3,k)**2*u(i,j+3,k,2)+mod%vel(i,j-3,k)**2*u(i,j-3,k,2))+ &
             &      coefs%c4x*(mod%vel(i,j+4,k)**2*u(i,j+4,k,2)+mod%vel(i,j-4,k)**2*u(i,j-4,k,2))
             tmpyy= coefs%c0y*(mod%vel(i  ,j,k)**2*u(i,j,k  ,2)) +&
             &      coefs%c1y*(mod%vel(i,j,k+1)**2*u(i,j,k+1,2)+mod%vel(i,j,k-1)**2*u(i,j,k-1,2))+ &
             &      coefs%c2y*(mod%vel(i,j,k+2)**2*u(i,j,k+2,2)+mod%vel(i,j,k-2)**2*u(i,j,k-2,2))+ &
             &      coefs%c3y*(mod%vel(i,j,k+3)**2*u(i,j,k+3,2)+mod%vel(i,j,k-3)**2*u(i,j,k-3,2))+ &
             &      coefs%c4y*(mod%vel(i,j,k+4)**2*u(i,j,k+4,2)+mod%vel(i,j,k-4)**2*u(i,j,k-4,2))
             u(i,j,k,3)=tmpxx+tmpzz+tmpyy
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine FD_2nd_3D_derivatives_scalar_adjoint
  
  subroutine FD_2nd_3D_derivatives_scalar_adjoint_grid(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j,k
    real, allocatable, dimension(:,:,:) :: tmp

    allocate(tmp(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
    tmp=mod%vel**2*u2

    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz= coefs%c0z*tmp(i,j,k) +&
             &      coefs%c1z*(tmp(i+1,j,k)+tmp(i-1,j,k))+ &
             &      coefs%c2z*(tmp(i+2,j,k)+tmp(i-2,j,k))+ &
             &      coefs%c3z*(tmp(i+3,j,k)+tmp(i-3,j,k))+ &
             &      coefs%c4z*(tmp(i+4,j,k)+tmp(i-4,j,k))
             tmpxx= coefs%c0x*tmp(i,j,k) +&
             &      coefs%c1x*(tmp(i,j+1,k)+tmp(i,j-1,k))+ &
             &      coefs%c2x*(tmp(i,j+2,k)+tmp(i,j-2,k))+ &
             &      coefs%c3x*(tmp(i,j+3,k)+tmp(i,j-3,k))+ &
             &      coefs%c4x*(tmp(i,j+4,k)+tmp(i,j-4,k))
             tmpyy= coefs%c0y*tmp(i,j,k) +&
             &      coefs%c1y*(tmp(i,j,k+1)+tmp(i,j,k-1))+ &
             &      coefs%c2y*(tmp(i,j,k+2)+tmp(i,j,k-2))+ &
             &      coefs%c3y*(tmp(i,j,k+3)+tmp(i,j,k-3))+ &
             &      coefs%c4y*(tmp(i,j,k+4)+tmp(i,j,k-4))
             u3(i,j,k)=tmpxx+tmpzz+tmpyy
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    deallocate(tmp)

  end subroutine FD_2nd_3D_derivatives_scalar_adjoint_grid
  
  subroutine FD_2nd_2D_derivatives_scalar_adjoint(genpar,bounds,u,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx
    integer              :: i,j
    
    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz= coefs%c0z*(mod%vel(i  ,j,1)**2*u(i  ,j,1,2)) +&
          &      coefs%c1z*(mod%vel(i+1,j,1)**2*u(i+1,j,1,2)+mod%vel(i-1,j,1)**2*u(i-1,j,1,2))+ &
          &      coefs%c2z*(mod%vel(i+2,j,1)**2*u(i+2,j,1,2)+mod%vel(i-2,j,1)**2*u(i-2,j,1,2))+ &
          &      coefs%c3z*(mod%vel(i+3,j,1)**2*u(i+3,j,1,2)+mod%vel(i-3,j,1)**2*u(i-3,j,1,2))+ &
          &      coefs%c4z*(mod%vel(i+4,j,1)**2*u(i+4,j,1,2)+mod%vel(i-4,j,1)**2*u(i-4,j,1,2))
          tmpxx= coefs%c0x*(mod%vel(i  ,j,1)**2*u(i,j  ,1,2)) +&
          &      coefs%c1x*(mod%vel(i,j+1,1)**2*u(i,j+1,1,2)+mod%vel(i,j-1,1)**2*u(i,j-1,1,2))+ &
          &      coefs%c2x*(mod%vel(i,j+2,1)**2*u(i,j+2,1,2)+mod%vel(i,j-2,1)**2*u(i,j-2,1,2))+ &
          &      coefs%c3x*(mod%vel(i,j+3,1)**2*u(i,j+3,1,2)+mod%vel(i,j-3,1)**2*u(i,j-3,1,2))+ &
          &      coefs%c4x*(mod%vel(i,j+4,1)**2*u(i,j+4,1,2)+mod%vel(i,j-4,1)**2*u(i,j-4,1,2))
          u(i,j,1,3)=real(tmpxx+tmpzz)
       end do
    end do

  end subroutine FD_2nd_2D_derivatives_scalar_adjoint
  
  subroutine FD_2nd_2D_derivatives_scalar_adjoint_grid(genpar,bounds,u2,u3,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    real                 :: tmpzz, tmpxx
    integer              :: i,j
    integer              :: ip1,ip2,ip3,ip4,jp1,jp2,jp3,jp4
    integer              :: im1,im2,im3,im4,jm1,jm2,jm3,jm4

    
    do j=bounds%nmin2,bounds%nmax2
       jp1=min(bounds%nmax2,j+1)
       jp2=min(bounds%nmax2,j+2)
       jp3=min(bounds%nmax2,j+3)
       jp4=min(bounds%nmax2,j+4)
       jm1=max(bounds%nmin2,j-1)
       jm2=max(bounds%nmin2,j-2)
       jm3=max(bounds%nmin2,j-3)
       jm4=max(bounds%nmin2,j-4)
       do i=bounds%nmin1,bounds%nmax1
          ip1=min(bounds%nmax1,i+1)
          ip2=min(bounds%nmax1,i+2)
          ip3=min(bounds%nmax1,i+3)
          ip4=min(bounds%nmax1,i+4)
          im1=max(bounds%nmin1,i-1)
          im2=max(bounds%nmin1,i-2)
          im3=max(bounds%nmin1,i-3)
          im4=max(bounds%nmin1,i-4)
          tmpzz= coefs%c0z*(mod%vel(i  ,j,1)**2*u2(i  ,j,1)) +&
          &      coefs%c1z*(mod%vel(ip1,j,1)**2*u2(i+1,j,1)+mod%vel(im1,j,1)**2*u2(i-1,j,1))+ &
          &      coefs%c2z*(mod%vel(ip2,j,1)**2*u2(i+2,j,1)+mod%vel(im2,j,1)**2*u2(i-2,j,1))+ &
          &      coefs%c3z*(mod%vel(ip3,j,1)**2*u2(i+3,j,1)+mod%vel(im3,j,1)**2*u2(i-3,j,1))+ &
          &      coefs%c4z*(mod%vel(ip4,j,1)**2*u2(i+4,j,1)+mod%vel(im4,j,1)**2*u2(i-4,j,1))
          tmpxx= coefs%c0x*(mod%vel(i  ,j,1)**2*u2(i,j  ,1)) +&
          &      coefs%c1x*(mod%vel(i,jp1,1)**2*u2(i,j+1,1)+mod%vel(i,jm1,1)**2*u2(i,j-1,1))+ &
          &      coefs%c2x*(mod%vel(i,jp2,1)**2*u2(i,j+2,1)+mod%vel(i,jm2,1)**2*u2(i,j-2,1))+ &
          &      coefs%c3x*(mod%vel(i,jp3,1)**2*u2(i,j+3,1)+mod%vel(i,jm3,1)**2*u2(i,j-3,1))+ &
          &      coefs%c4x*(mod%vel(i,jp4,1)**2*u2(i,j+4,1)+mod%vel(i,jm4,1)**2*u2(i,j-4,1))
          u3(i,j,1)=real(tmpxx+tmpzz)
       end do
    end do

  end subroutine FD_2nd_2D_derivatives_scalar_adjoint_grid
  
  subroutine FD_2nd_time_derivative_omp(genpar,bounds,u)
    type(GeneralParam)   ::             genpar
    type(FDbounds)       ::                    bounds
    real                 ::                           u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    integer              :: i,j,k
    real                 :: div11,dt2

    div11=1./11
    dt2=genpar%dt**2
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             u(i,j,k,3) = ( 20.*u(i,j,k,2) - 6.*u(i,j,k,1) - &
             &               4.*u(i,j,k,0) +    u(i,j,k,-1) + &
             &              12.*dt2*u(i,j,k,3) ) * div11
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine FD_2nd_time_derivative_omp
  
  subroutine FD_2nd_time_derivative(genpar,bounds,u)
    type(GeneralParam)   ::         genpar
    type(FDbounds)       ::                bounds
    real                 ::                       u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    integer              :: i,j,k
    real                 :: div11,dt2

    div11=1./11
    dt2=genpar%dt**2

    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             u(i,j,k,3) = ( 20.*u(i,j,k,2) - 6.*u(i,j,k,1) - &
             &               4.*u(i,j,k,0) +    u(i,j,k,-1) + &
             &              12.*dt2*u(i,j,k,3) ) * div11
          end do
       end do
    end do

  end subroutine FD_2nd_time_derivative
  
  subroutine FD_2nd_time_derivative_grid_noomp(genpar,bounds,grid)
    type(GeneralParam)   ::         genpar
    type(FDbounds)       ::                bounds
    type(USpace)         ::                        grid
    integer              :: i,j,k
    real                 :: div11,dt2

    div11=1./11
    dt2=genpar%dt**2

    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             grid%u3(i,j,k) = ( 20.*grid%u2(i,j,k) - 6.*grid%u1(i,j,k) - &
             &                   4.*grid%u0(i,j,k) +    grid%u_1(i,j,k) + &
             &              12.*dt2*grid%u3(i,j,k) ) * div11
          end do
       end do
    end do

  end subroutine FD_2nd_time_derivative_grid_noomp
  
  subroutine FD_2nd_time_derivative_grid(genpar,bounds,grid)
    type(GeneralParam)   ::         genpar
    type(FDbounds)       ::                bounds
    type(USpace)         ::                        grid
    integer              :: i,j,k
    real                 :: div11,dt2

    div11=1./11
    dt2=genpar%dt**2

    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             grid%u3(i,j,k) = ( 20.*grid%u2(i,j,k) - 6.*grid%u1(i,j,k) - &
             &                   4.*grid%u0(i,j,k) +    grid%u_1(i,j,k) + &
             &              12.*dt2*grid%u3(i,j,k) ) * div11
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine FD_2nd_time_derivative_grid
  
  subroutine FD_2nd_time_derivative_grid2(genpar,bounds,grid)
    type(GeneralParam)   ::         genpar
    type(FDbounds)       ::                bounds
    type(USpace)         ::                        grid
    integer              :: i,j,k,ii,ij,ik,kstep,jstep,istep
    real                 :: div11,dt2

    div11=1./11
    dt2=genpar%dt**2

    kstep=2
    jstep=4
    istep=8

    !$OMP PARALLEL DO PRIVATE(ik,ij,ii)

    do ik=bounds%nmin3,bounds%nmax3,kstep
       do ij=bounds%nmin2,bounds%nmax2,jstep
          do ii=bounds%nmin1,bounds%nmax1,istep

             do k=ik,min(bounds%nmax3,ik+kstep-1)
                do j=ij,min(bounds%nmax2,ij+jstep-1)
                   do i=ii,min(bounds%nmax1,ii+istep-1)

                      grid%u3(i,j,k) = ( 20.*grid%u2(i,j,k) - 6.*grid%u1(i,j,k) - &
                      &                   4.*grid%u0(i,j,k) +    grid%u_1(i,j,k) + &
                      &              12.*dt2*grid%u3(i,j,k) ) * div11

                   end do
                end do
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine FD_2nd_time_derivative_grid2
  
end module FD_derivatives
