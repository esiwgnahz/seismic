module FD_derivatives

  use FD_types
  
  type(ScaledFDcoefs), pointer, private  ::                      coefs
  
  implicit none

contains

  subroutine FD_derivatives_coef_init(coef_init)
    type(ScaledFDcoefs), target :: coef_init
    coefs => coef_init
  end subroutine FD_derivatives_coef_init

  subroutine FD_2nd_3D_derivatives_scalar_forward(bounds,u,v2)
    type(FDbounds)       ::                       bounds
    real                 :: u(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3,-1,3)

    real                 ::v2(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)

    real                 :: tmpzz, tmpxx, tmpyy

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
             u(i,j,k,3)=v2(i,j,k)*(tmpxx+tmpzz+tmpyy)
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine FD_2nd_3D_derivatives_scalar_forward
  
  subroutine FD_2nd_2D_derivatives_scalar_forward(bounds,u,v2)
    type(FDbounds)       ::                       bounds
    real                 :: u(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3,-1,3)

    real                 ::v2(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)

    real                 :: tmpzz, tmpxx, tmpyy

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
          u(i,j,1,3)=v2(i,j,1)*(tmpxx+tmpzz)
       end do
    end do
    !$OMP END PARALLEL DO
    
  end subroutine FD_2nd_2D_derivatives_scalar_forward
  
  subroutine FD_2nd_3D_derivatives_scalar_adjoint(bounds,u,v2)
    type(FDbounds)       ::                       bounds
    real                 :: u(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3,-1,3)
    
    real                 ::v2(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)
    
    real                 :: tmpzz, tmpxx, tmpyy
    
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             tmpzz= coefs%c0z*(v2(i  ,j,k)*u(i  ,j,k,2)) +&
             &      coefs%c1z*(v2(i+1,j,k)*u(i+1,j,k,2)+v2(i-1,j,k)*u(i-1,j,k,2))+ &
             &      coefs%c2z*(v2(i+2,j,k)*u(i+2,j,k,2)+v2(i-2,j,k)*u(i-2,j,k,2))+ &
             &      coefs%c3z*(v2(i+3,j,k)*u(i+3,j,k,2)+v2(i-3,j,k)*u(i-3,j,k,2))+ &
             &      coefs%c4z*(v2(i+4,j,k)*u(i+4,j,k,2)+v2(i-4,j,k)*u(i-4,j,k,2))
             tmpxx= coefs%c0x*(v2(i  ,j,k)*u(i,j  ,k,2)) +&
             &      coefs%c1x*(v2(i,j+1,k)*u(i,j+1,k,2)+v2(i,j-1,k)*u(i,j-1,k,2))+ &
             &      coefs%c2x*(v2(i,j+2,k)*u(i,j+2,k,2)+v2(i,j-2,k)*u(i,j-2,k,2))+ &
             &      coefs%c3x*(v2(i,j+3,k)*u(i,j+3,k,2)+v2(i,j-3,k)*u(i,j-3,k,2))+ &
             &      coefs%c4x*(v2(i,j+4,k)*u(i,j+4,k,2)+v2(i,j-4,k)*u(i,j-4,k,2))
             tmpyy= coefs%c0y*(v2(i  ,j,k)*u(i,j,k  ,2)) +&
             &      coefs%c1y*(v2(i,j,k+1)*u(i,j,k+1,2)+v2(i,j,k-1)*u(i,j,k-1,2))+ &
             &      coefs%c2y*(v2(i,j,k+2)*u(i,j,k+2,2)+v2(i,j,k-2)*u(i,j,k-2,2))+ &
             &      coefs%c3y*(v2(i,j,k+3)*u(i,j,k+3,2)+v2(i,j,k-3)*u(i,j,k-3,2))+ &
             &      coefs%c4y*(v2(i,j,k+4)*u(i,j,k+4,2)+v2(i,j,k-4)*u(i,j,k-4,2))
             u(i,j,k,3)=tmpxx+tmpzz+tmpyy
          end do
       end do
    end do
  end subroutine FD_2nd_3D_derivatives_scalar_adjoint
  
  subroutine FD_2nd_2D_derivatives_scalar_adjoint(bounds,u,v2)
    type(FDbounds)       ::                       bounds
    real                 :: u(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3,-1,3)
    
    real                 ::v2(bounds%nmin1:bounds%nmax1,bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)
    
    real                 :: tmpzz, tmpxx
    
    do j=bounds%nmin2,bounds%nmax2
       do i=bounds%nmin1,bounds%nmax1
          tmpzz= coefs%c0z*(v2(i  ,j,1)*u(i  ,j,1,2)) +&
          &      coefs%c1z*(v2(i+1,j,1)*u(i+1,j,1,2)+v2(i-1,j,1)*u(i-1,j,1,2))+ &
          &      coefs%c2z*(v2(i+2,j,1)*u(i+2,j,1,2)+v2(i-2,j,1)*u(i-2,j,1,2))+ &
          &      coefs%c3z*(v2(i+3,j,1)*u(i+3,j,1,2)+v2(i-3,j,1)*u(i-3,j,1,2))+ &
          &      coefs%c4z*(v2(i+4,j,1)*u(i+4,j,1,2)+v2(i-4,j,1)*u(i-4,j,1,2))
          tmpxx= coefs%c0x*(v2(i  ,j,1)*u(i,j  ,1,2)) +&
          &      coefs%c1x*(v2(i,j+1,1)*u(i,j+1,1,2)+v2(i,j-1,1)*u(i,j-1,1,2))+ &
          &      coefs%c2x*(v2(i,j+2,1)*u(i,j+2,1,2)+v2(i,j-2,1)*u(i,j-2,1,2))+ &
          &      coefs%c3x*(v2(i,j+3,1)*u(i,j+3,1,2)+v2(i,j-3,1)*u(i,j-3,1,2))+ &
          &      coefs%c4x*(v2(i,j+4,1)*u(i,j+4,1,2)+v2(i,j-4,1)*u(i,j-4,1,2))
          u(i,j,1,3)=tmpxx+tmpzz
       end do
    end do

  end subroutine FD_2nd_2D_derivatives_scalar_adjoint
  
end module FD_derivatives
