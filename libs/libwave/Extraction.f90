module Extraction_mod

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod

  implicit none

contains

  subroutine Extraction_shot_simple(bounds,dat,elev,u,genpar,it)
    type(FDbounds)    ::            bounds
    type(DataSpace)   ::                   dat
    type(ModelSpace_elevation) ::              elev
    type(GeneralParam)::                              genpar
    real              ::                            u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    integer           :: it,k,j
    

    if (genpar%rec_type.eq.0) then
       do k=1,dat%ny
          do j=1,dat%nx
             dat%data(it,j,k)=u(elev%irec_z(j,k),j,k,2)
          end do
       end do
    else
       do k=1,dat%ny
          do j=1,dat%nx
             dat%data(it,j,k)=u(elev%irec_z(j,k),j,k,2)-u(-elev%irec_z(j,k),j,k,2)
          end do
       end do
    end if

  end subroutine Extraction_shot_simple

  subroutine Extraction_shot_sinc(bounds,dat,elev,u,genpar,it)
    type(FDbounds)    ::          bounds
    type(DataSpace)   ::                 dat
    type(ModelSpace_elevation) ::              elev
    type(GeneralParam)::                            genpar
    real              ::                     u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                          bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    integer           :: it,k,i,j,l
    real              :: fdum
    real,allocatable  :: sinc(:)
    
    allocate(sinc(genpar%lsinc))

    if (genpar%rec_type.eq.0) then
       do k=1,dat%ny
          do j=1,dat%nx
             call mksinc(sinc,genpar%lsinc,elev%drec_z(j,k)/genpar%dz)
             fdum = 0.
             do l=-genpar%lsinc/2,genpar%lsinc/2
                fdum = fdum + sinc(genpar%lsinc/2+1+l)* &
             &  u(elev%irec_z(j,k)+l,j,k,2)
             end do
             dat%data(it,j,k) = fdum
          end do
       end do
    else
       do k=1,dat%ny
          do j=1,dat%nx
             call mksinc(sinc,genpar%lsinc,elev%drec_z(j,k)/genpar%dz)
             fdum = 0.
             do l=-genpar%lsinc/2,genpar%lsinc/2
                fdum = fdum + sinc(genpar%lsinc/2+1+l)* &
             &  (u(elev%irec_z(j,k)+l,j,k,2) - u(-elev%irec_z(j,k)+l,j,k,2))
             end do
             dat%data(it,j,k) = fdum
          end do
       end do
    end if
       
    deallocate(sinc)

  end subroutine Extraction_shot_sinc

  subroutine Extraction_wavefield(bounds,model,dat,elev,u,genpar,it,counter)
    type(FDbounds)    ::          bounds
    type(ModelSpace)  ::                 model
    type(DataSpace)   ::                     dat
    type(ModelSpace_elevation) ::                elev
    type(GeneralParam)::                                genpar
    real              ::                              u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    integer           :: it,i,k,j,counter
    
    real, allocatable :: buffer_sou(:,:)

    allocate(buffer_sou(model%nz,model%nx))

    MODULO:if (mod(it-1,genpar%snapi).eq.0) then
       counter=counter+1
       if (genpar%surf_type.ne.0) then
          do k=1,model%ny
             buffer_sou = 0.
             do j=1,model%nx
                do i=elev%ielev_z(j,k),model%nz
                   buffer_sou(i,j) = u(i,j,k,2)
                end do
             end do
             dat%wave(:,:,k,counter)=buffer_sou
          end do
       else
          dat%wave(:,:,:,counter)=u(1:model%nz,1:model%nx,1:model%ny,2)         
       end if
    end if MODULO
    
    deallocate(buffer_sou)

  end subroutine Extraction_wavefield

!  subroutine Extraction_shot_simple(bounds,dat,u,irec_z,drec_z,genpar,it)
!    type(FDbounds)    ::            bounds
!    type(DataSpace)   ::                   dat
!    type(GeneralParam):: genpar
!    real              ::                       u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
!    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
!    real              :: drec_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)
!    integer           :: irec_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)
!    integer           :: it,k,j
!    
!
!    if (genpar%rec_type.eq.0) then
!       do k=1,dat%ny
!          do j=1,dat%nx
!             dat%data(it,j,k)=u(irec_z(j,k),j,k,2)
!          end do
!       end do
!    else
!       do k=1,dat%ny
!          do j=1,dat%nx
!             dat%data(it,j,k)=u(irec_z(j,k),j,k,2)-u(-irec_z(j,k),j,k,2)
!          end do
!       end do
!    end if
!
!  end subroutine Extraction_shot_simple
!
!  subroutine Extraction_shot_sinc(bounds,dat,u,irec_z,drec_z,genpar,it)
!    type(FDbounds)    ::          bounds
!    type(DataSpace)   ::                 dat
!    type(GeneralParam):: genpar
!    real              ::                     u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
!    &                                          bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
!    integer           :: irec_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)
!    real              :: drec_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3)
!    integer           :: it,k,i,j,l
!    real              :: fdum
!    real,allocatable  :: sinc(:)
!    
!    allocate(sinc(genpar%lsinc))
!
!    if (genpar%rec_type.eq.0) then
!       do k=1,dat%ny
!          do j=1,dat%nx
!             call mksinc(sinc,genpar%lsinc,drec_z(j,k)/genpar%dz)
!             fdum = 0.
!             do l=-genpar%lsinc/2,genpar%lsinc/2
!                fdum = fdum + sinc(genpar%lsinc/2+1+l)* &
!             &  u(irec_z(j,k)+l,j,k,2)
!             end do
!             dat%data(it,j,k) = fdum
!          end do
!       end do
!    else
!       do k=1,dat%ny
!          do j=1,dat%nx
!             call mksinc(sinc,genpar%lsinc,drec_z(j,k)/genpar%dz)
!             fdum = 0.
!             do l=-genpar%lsinc/2,genpar%lsinc/2
!                fdum = fdum + sinc(genpar%lsinc/2+1+l)* &
!             &  (u(irec_z(j,k)+l,j,k,2) - u(-irec_z(j,k)+l,j,k,2))
!             end do
!             dat%data(it,j,k) = fdum
!          end do
!       end do
!    end if
!       
!    deallocate(sinc)
!
!  end subroutine Extraction_shot_sinc

end module Extraction_mod
