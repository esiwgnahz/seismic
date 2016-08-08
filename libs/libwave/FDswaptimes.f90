module FDswaptime_mod

  use GeneralParam_types
  use FD_types
  implicit none
  contains
   
    subroutine FDswaptime_omp(genpar,bounds,u)
      type(GeneralParam)  ::  genpar
      type(FDbounds)     ::          bounds
      real               ::                 u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
      integer         :: k,j,i

      !$OMP PARALLEL DO PRIVATE(k,j,i)
      do k=bounds%nmin3-genpar%nbound,bounds%nmax3+genpar%nbound
         do j=bounds%nmin2-4,bounds%nmax2+4
            do i=bounds%nmin1-4,bounds%nmax1+4
               u(i,j,k,-1) = u(i,j,k,0)
               u(i,j,k,0)  = u(i,j,k,1)
               u(i,j,k,1)  = u(i,j,k,2)
               u(i,j,k,2)  = u(i,j,k,3)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

    end subroutine FDswaptime_omp
    
    subroutine FDswaptime(genpar,bounds,u)
      type(GeneralParam)::genpar
      type(FDbounds)    ::       bounds
      real              ::              u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
      integer         :: k,j,i

      do k=bounds%nmin3-genpar%nbound,bounds%nmax3+genpar%nbound
         do j=bounds%nmin2-4,bounds%nmax2+4
            do i=bounds%nmin1-4,bounds%nmax1+4
               u(i,j,k,-1) = u(i,j,k,0)
               u(i,j,k,0)  = u(i,j,k,1)
               u(i,j,k,1)  = u(i,j,k,2)
               u(i,j,k,2)  = u(i,j,k,3)
            end do
         end do
      end do

    end subroutine FDswaptime
    

end module FDswaptime_mod
