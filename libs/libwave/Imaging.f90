module Imaging_mod
  
  use sep

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod

  use ExtractPadModel_mod

  implicit none

contains

  subroutine Imaging_condition_from_disk(model,genpar)
    type(ModelSpace) :: model
    type(GeneralParam) :: genpar
    integer :: i,j,k,l,counter
    real, dimension (:,:,:), allocatable :: fwd,bwd,tmpim,tmpil
 
    real    :: maxillu,taper
    integer :: ierr,blocksize,index

    allocate(tmpim(model%nz,model%nxw,model%nyw))
    allocate(tmpil(model%nz,model%nxw,model%nyw))

    tmpim=0.
    tmpil=0.

    allocate(fwd(model%nz,model%nxw,model%nyw))
    allocate(bwd(model%nz,model%nxw,model%nyw))

    blocksize=4*model%nz*model%nxw*model%nyw
    call sseek(model%waFtag,0,0)
    call sseek(model%waBtag,0,0)
    do l=1,genpar%ntsnap

       fwd=0.
       bwd=0.

       index=genpar%ntsnap-l
       call sseek_block(model%waFtag,index,blocksize,0)
       call sreed(model%waBtag,bwd,blocksize)
       call sreed(model%waFtag,fwd,blocksize)
       
       !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
       do k=1,model%nyw
          do j=1,model%nxw
             taper=model%taperx(j)*model%tapery(k)
             do i=1,model%nz
                tmpim(i,j,k)=tmpim(i,j,k)+fwd(i,j,k)*bwd(i,j,k)*taper
                tmpil(i,j,k)=tmpil(i,j,k)+fwd(i,j,k)*fwd(i,j,k)*taper
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end do

    maxillu=maxval(tmpil)
    tmpil=(tmpil+maxillu*1e-4)/rms(tmpil,model%nz,model%nxw,model%nyw)/(model%nz*model%nxw*model%nyw)

    deallocate(fwd,bwd)
    
    allocate(model%image(model%nz,model%nx,model%ny))
    allocate(model%illum(model%nz,model%nx,model%ny))
    
    model%image=0.
    model%illum=0.

    call mod_window_pad(.false.,model%image,tmpim,model)
    call mod_window_pad(.false.,model%illum,tmpil,model)

    deallocate(tmpim,tmpil)

  end subroutine Imaging_condition_from_disk

  subroutine scale_illumination(mod)
    type(ModelSpace) ::  mod
    real :: maxillu

    maxillu=maxval(mod%illumsmall)
    mod%illumsmall=(mod%illumsmall+maxillu*1e-4)/rms(mod%illumsmall,mod%nz,mod%nxw,mod%nyw)/(mod%nz*mod%nxw*mod%nyw)
  
  end subroutine scale_illumination

  real function rms(array,nz,nx,ny)
    real, dimension(nz,nx,ny) :: array
    integer :: nx,nz,ny
    integer :: i,j,k
    
    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=1,ny
       do j=1,nx
          do i=1,nz
             rms=rms+array(i,j,k)**2      
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    rms=sqrt(rms)

  end function rms

  subroutine Imaging_condition_sourceonly_from_disk(bounds,model,elev,u,genpar,it)
    type(FDbounds)    ::                            bounds
    type(ModelSpace)  ::                                   model
    type(ModelSpace_elevation) ::                                elev
    type(GeneralParam)::                                                genpar
    real              ::                                             u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer :: i,j,k,counter
    real, dimension (:,:,:), allocatable :: fwd

    integer :: ierr,blocksize,index
    real    :: taper

    MODULO:if (mod(it,genpar%snapi).eq.0) then
       model%counter=model%counter+1
       ! First read the source wavefield from disk
       allocate(fwd(model%nz,model%nxw,model%nyw))
       blocksize=4*model%nz*model%nxw*model%nyw
       call sseek(model%waFtag,0,0)
       fwd=0.
       index=genpar%ntsnap-model%counter
       call sseek_block(model%waFtag,index,blocksize,0)
       call sreed(model%waFtag,fwd,blocksize)

       if (genpar%surf_type.ne.0) then
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=elev%ielev_z(j,k),model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+u(i,j,k)*fwd(i,j,k)*taper
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+fwd(i,j,k)*fwd(i,j,k)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=1,model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+u(i,j,k)*fwd(i,j,k)*taper
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+fwd(i,j,k)*fwd(i,j,k)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if                  

       deallocate(fwd)

    end if MODULO

  end subroutine Imaging_condition_sourceonly_from_disk

  subroutine Imaging_condition_LSRTM_sourceonly_from_disk(bounds,model,elev,u,genpar,it)
    type(FDbounds)    ::                                  bounds
    type(ModelSpace)  ::                                         model
    type(ModelSpace_elevation) ::                                      elev
    type(GeneralParam)::                                                      genpar
    real              ::                                                   u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer :: i,j,k,counter
    real, dimension (:,:,:), allocatable :: fwd

    integer :: ierr,blocksize,index
    real    :: taper

    MODULO:if (mod(it,genpar%snapi).eq.0) then
       model%counter=model%counter+1
       ! First read the source wavefield from disk
       allocate(fwd(model%nz,model%nxw,model%nyw))
       blocksize=4*model%nz*model%nxw*model%nyw
       call sseek(model%waFtag,0,0)
       fwd=0.
       index=genpar%ntsnap-model%counter
       call sseek_block(model%waFtag,index,blocksize,0)
       call sreed(model%waFtag,fwd,blocksize)

       if (genpar%surf_type.ne.0) then
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=elev%ielev_z(j,k),model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+2*u(i,j,k)*fwd(i,j,k)*taper/model%vel(i,j,k)
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+fwd(i,j,k)*fwd(i,j,k)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=1,model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+2*u(i,j,k)*fwd(i,j,k)*taper/model%vel(i,j,k)
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+fwd(i,j,k)*fwd(i,j,k)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if                  

       deallocate(fwd)

    end if MODULO

  end subroutine Imaging_condition_LSRTM_sourceonly_from_disk

  subroutine Imaging_condition_sourceonly_from_memory(bounds,model,elev,u,genpar,it)
    type(FDbounds)    ::                              bounds
    type(ModelSpace)  ::                                     model
    type(ModelSpace_elevation) ::                                  elev
    type(GeneralParam)::                                                  genpar
    real              ::                                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer :: i,j,k,counter
    real, dimension (:,:,:), pointer :: bwd

    integer :: ierr,blocksize,index
    real    :: taper
    
    MODULO:if (mod(it,genpar%snapi).eq.0) then
       
       model%counter=model%counter+1
       index=genpar%ntsnap-model%counter+1

       if (genpar%surf_type.ne.0) then
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=elev%ielev_z(j,k),model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+u(i,j,k)*model%wvfld%wave(i,j,k,index,1)*taper
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+model%wvfld%wave(i,j,k,index,1)*model%wvfld%wave(i,j,k,index,1)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=1,model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+u(i,j,k)*model%wvfld%wave(i,j,k,index,1)*taper
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+model%wvfld%wave(i,j,k,index,1)*model%wvfld%wave(i,j,k,index,1)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if                  

    end if MODULO

  end subroutine Imaging_condition_sourceonly_from_memory

  subroutine Imaging_condition_LSRTM_sourceonly_from_memory(bounds,model,elev,u,genpar,it)
    type(FDbounds)    ::                                   bounds
    type(ModelSpace)  ::                                          model
    type(ModelSpace_elevation) ::                                      elev
    type(GeneralParam)::                                                      genpar
    real              ::                                                    u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, & 
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer :: i,j,k,counter
    real, dimension (:,:,:), pointer :: bwd

    integer :: ierr,blocksize,index
    real    :: taper
    
    MODULO:if (mod(it,genpar%snapi).eq.0) then
       
       model%counter=model%counter+1
       index=genpar%ntsnap-model%counter+1

       if (genpar%surf_type.ne.0) then
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=elev%ielev_z(j,k),model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+2*taper*u(i,j,k)*model%wvfld%wave(i,j,k,index,1)/model%vel(i,j,k)
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+model%wvfld%wave(i,j,k,index,1)*model%wvfld%wave(i,j,k,index,1)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=1,model%nz
                   model%imagesmall(i,j,k)= model%imagesmall(i,j,k)+2*taper*u(i,j,k)*model%wvfld%wave(i,j,k,index,1)/model%vel(i,j,k)
                   model%illumsmall(i,j,k)= model%illumsmall(i,j,k)+model%wvfld%wave(i,j,k,index,1)*model%wvfld%wave(i,j,k,index,1)*taper
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if                  

    end if MODULO

  end subroutine Imaging_condition_LSRTM_sourceonly_from_memory

  subroutine Injection_Born(bounds,model,elev,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(ModelSpace_elevation) ::        elev
    type(GeneralParam)::                         genpar 
    real              ::                      u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it,i,j,k   
    real              :: taper

    if (mod(it,genpar%snapi).eq.0) then

       model%counter=model%counter+1

       if (genpar%surf_type.ne.0) then
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=elev%ielev_z(j,k),model%nz
                   u(i,j,k)=u(i,j,k)+2*taper*model%image(i,j,k)*model%wvfld%wave(i,j,k,model%counter,1)/model%vel(i,j,k)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=1,model%nz
                   u(i,j,k)=u(i,j,k)+2*taper*model%image(i,j,k)*model%wvfld%wave(i,j,k,model%counter,1)/model%vel(i,j,k)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if                  

    end if 
    
  end subroutine Injection_Born

  subroutine Injection_Born_from_disk(bounds,model,elev,u,genpar,it)
    
    type(FDbounds)    ::              bounds
    type(ModelSpace)  ::                     model
    type(ModelSpace_elevation) ::                  elev
    type(GeneralParam)::                                  genpar 
    real              ::                                u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it   
    real, dimension (:,:,:), allocatable :: fwd
    integer           :: ierr,i,j,k,blocksize
    real              :: taper

    allocate(fwd(model%nz,model%nxw,model%nyw))
    fwd=0.

    blocksize=4*model%nz*model%nxw*model%nyw

    if (mod(it,genpar%snapi).eq.0) then

       call sreed(model%waFtag,fwd,blocksize)       

       if (genpar%surf_type.ne.0) then
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=elev%ielev_z(j,k),model%nz
                   u(i,j,k)=u(i,j,k)+2*taper*model%image(i,j,k)*fwd(i,j,k)/model%vel(i,j,k)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO PRIVATE(k,j,taper,i)
          do k=1,model%nyw
             do j=1,model%nxw
                taper=model%taperx(j)*model%tapery(k)
                do i=1,model%nz
                   u(i,j,k)=u(i,j,k)+2*taper*model%image(i,j,k)*fwd(i,j,k)/model%vel(i,j,k)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if                  

    end if
    
    deallocate(fwd)

  end subroutine Injection_Born_from_disk

end module Imaging_mod
