module FD_derivatives

  use FD_types
  use ModelSpace_types
  use GeneralParam_types
  
  implicit none

  type(ScaledFDcoefs), pointer, private  :: coefs
  

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
    allocate(delp(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
    
    delp=1./mod%rho2
    dxi =1./mod%dx
    dyi =1./mod%dy
    dzi =1./mod%dz

    !$OMP PARALLEL DO PRIVATE(k,j,i,sxx,szz,syy)
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
    do k=bounds%nmin3-1,bounds%nmin3-3,-1
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmin3)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmin3)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(i,j,bounds%nmin3)
          end do
       end do
    end do
    ! Gradient of the p-wave potential in the back side operator padding 
    do k=bounds%nmax3+1,bounds%nmax3+3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             szz(i,j,k)=dzi*(u(i+1,j,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmax3)
             sxx(i,j,k)=dxi*(u(i,j+1,k,2)-u(i,j,k,2))*delp(i,j,bounds%nmax3)
             syy(i,j,k)=dyi*(u(i,j,k+1,2)-u(i,j,k,2))*delp(i,j,bounds%nmax3)
          end do
       end do
    end do
    
    !$OMP PARALLEL DO PRIVATE(k,j,i,tmpxx,tmpzz,tmpyy)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz=         coefs%c1z*(szz(i,j,k)-szz(i-1,j,k))+ &
             &              coefs%c2z*(szz(i+1,j,k)-szz(i-2,j,k))+ &
             &              coefs%c3z*(szz(i+2,j,k)-szz(i-3,j,k))+ &
             &              coefs%c4z*(szz(i+3,j,k)-szz(i-4,j,k))
             tmpxx=         coefs%c1x*(sxx(i,j,k)-sxx(i,j-1,k))+ &
             &              coefs%c2x*(sxx(i,j+1,k)-sxx(i,j-2,k))+ &
             &              coefs%c3x*(sxx(i,j+2,k)-sxx(i,j-3,k))+ &
             &              coefs%c4x*(sxx(i,j+3,k)-sxx(i,j-4,k))
             tmpyy=         coefs%c1y*(syy(i,j,k)-syy(i,j,k-1))+ &
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
    allocate(delp(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
    
    delp=1./mod%rho2
    dxi =1./mod%dx
    dzi =1./mod%dz

    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1   
          szz(i,j,1)=           (coefs%c1z*(u(i+1,j,1,2)-u(i,j,1,2))+ &
          &                      coefs%c2z*(u(i+2,j,1,2)-u(i-1,j,1,2))+ &
          &                      coefs%c3z*(u(i+3,j,1,2)-u(i-2,j,1,2))+ &
          &                      coefs%c4z*(u(i+4,j,1,2)-u(i-3,j,1,2)))*delp(i,j,1)
          sxx(i,j,1)=           (coefs%c1x*(u(i,j+1,1,2)-u(i,j,1,2))+ &
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
          u(i,j,1,3)=mod%vel(i,j,1)**2* &
          &              mod%rho(i,j,1)*(tmpxx+tmpzz)
       end do
    end do
    deallocate(sxx,szz,delp)

  end subroutine FD_2D_derivatives_acoustic_forward
  
  subroutine FD_2nd_2D_derivatives_scalar_forward(genpar,bounds,u,mod)
    type(GeneralParam)   ::                       genpar
    type(ModelSpace)     ::                                       mod
    type(FDbounds)       ::                              bounds
    real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    real                 :: tmpzz, tmpxx, tmpyy
    integer              :: i,j

    !$OMP PARALLEL DO PRIVATE(j,i,tmpxx,tmpzz)
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
    !$OMP END PARALLEL DO
    
  end subroutine FD_2nd_2D_derivatives_scalar_forward
  
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
          u(i,j,1,3)=tmpxx+tmpzz
       end do
    end do

  end subroutine FD_2nd_2D_derivatives_scalar_adjoint
  
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
  
end module FD_derivatives
