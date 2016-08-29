module Extraction_mod

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod

  implicit none

contains

  ! Nearest neighbor
  subroutine Extraction_array_simple(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::             bounds
    type(ModelSpace)  ::                    model
    type(TraceSpace), dimension(:) ::             data
    type(GeneralParam)::                                 genpar 
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    
    integer           :: i
    
    real, allocatable :: deltai(:)


    allocate(deltai(3))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)
    deltai(3)=1./genpar%delta(3)

    if (genpar%twoD) then
       if (genpar%rec_type.eq.0) then
          do i=1,size(data)
             data(i)%trace(it,1)=data(i)%trace(it,1)+model%vel(data(i)%icoord(1),data(i)%icoord(2),1)**2*&
             &                   u(data(i)%icoord(1),data(i)%icoord(2),1)*product(deltai(1:2))
          end do
       else
          do i=1,size(data)
             data(i)%trace(it,1)=data(i)%trace(it,1)+model%vel(data(i)%icoord(1),data(i)%icoord(2),1)**2*&
             &                   (u(data(i)%icoord(1),data(i)%icoord(2),1)-u(-data(i)%icoord(1),data(i)%icoord(2),1))*product(deltai(1:2))
          end do
       end if
    else
       if (genpar%rec_type.eq.0) then
          !$OMP PARALLEL DO PRIVATE(i)
          do i=1,size(data)
             data(i)%trace(it,1)=data(i)%trace(it,1)+model%vel(data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3))**2*&
             &                   u(data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3))*product(deltai)  
          end do
          !$OMP END PARALLEL DO    
       else
          !$OMP PARALLEL DO PRIVATE(i)
          do i=1,size(data)
             data(i)%trace(it,1)=data(i)%trace(it,1)+model%vel(data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3))**2*&
             &                   (u(data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3))-u(-data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3)))*product(deltai) 
          end do
          !$OMP END PARALLEL DO 
       end if
    end if
    deallocate(deltai)

  end subroutine Extraction_array_simple

  ! Sinc interpolation
  subroutine Extraction_array_sinc(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::           bounds
    type(ModelSpace)  ::                  model
    type(TraceSpace), dimension(:) ::           data
    type(GeneralParam)::                               genpar 
    real              ::                             u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    
    if (genpar%twoD) then
       do i=1,size(data)
          call Extraction_1trace_sinc_xy(bounds,model%vel(data(i)%icoord(1),data(i)%icoord(2),1),data(i),u,genpar,it)
       end do
    else      
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)       
          call Extraction_1trace_sinc_xyz(bounds,model%vel(data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3)),data(i),u,genpar,it)
       end do       
       !$OMP END PARALLEL DO
    end if

  end subroutine Extraction_array_sinc

  ! Sinc interpolation
  subroutine Extraction_array_lint3(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::            bounds
    type(ModelSpace)  ::                   model
    type(TraceSpace), dimension(:) ::            data
    type(GeneralParam)::                                genpar 
    real              ::                              u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    
    if (genpar%twoD) then
       do i=1,size(data)
          call Extraction_1trace_lint2_xy(bounds,model%vel(data(i)%icoord(1),data(i)%icoord(2),1),data(i),u,genpar,it)
       end do
    else      
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)       
          call Extraction_1trace_lint3_xyz(bounds,model%vel(data(i)%icoord(1),data(i)%icoord(2),data(i)%icoord(3)),data(i),u,genpar,it)
       end do       
       !$OMP END PARALLEL DO
    end if

  end subroutine Extraction_array_lint3

  subroutine Extraction_1trace_sinc_xyz(bounds,vel,data,u,genpar,it)
    type(FDbounds)    ::                bounds
    real              ::                       vel
    type(TraceSpace)  ::                           data
    type(GeneralParam)::                                    genpar 
    real              ::                                  u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j,k
    integer           :: minx,maxx,miny,maxy,minz,maxz

    allocate(sinc(genpar%lsinc,3))
    allocate(deltai(3))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)
    deltai(3)=1./genpar%delta(3)

    minx=-genpar%lsinc*0.5
    minz=minx
    miny=minx
    
    maxx=genpar%lsinc*0.5
    maxz=maxx
    maxy=maxx

    do i=1,3
       call mksinc(sinc(:,i),genpar%lsinc,data%dcoord(i)*deltai(i))
    end do
    
    if (genpar%rec_type.eq.0) then
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                data%trace(it,1)=data%trace(it,1)+ &
                &   vel**2*u(i+data%icoord(1),j+data%icoord(2),k+data%icoord(3))*&
                &   sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*product(deltai) 
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    
    else
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                data%trace(it,1)=data%trace(it,1)+ &
                &   vel**2* &
                &   (u(i+data%icoord(1),j+data%icoord(2),k+data%icoord(3))-u(-i-data%icoord(1),j+data%icoord(2),k+data%icoord(3)))*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*product(deltai) 
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    deallocate(sinc,deltai)

  end subroutine Extraction_1trace_sinc_xyz

  subroutine Extraction_1trace_lint3_xyz(bounds,vel,data,u,genpar,it)
    type(FDbounds)    ::                 bounds
    real              ::                        vel
    type(TraceSpace)  ::                            data
    type(GeneralParam)::                                    genpar 
    real              ::                                 u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    integer           :: i,j

    real              :: f,utmp,utmpm
    real              :: fx(3),gx(3)
    integer           :: ix(3)

    allocate(deltai(3))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)
    deltai(3)=1./genpar%delta(3)
    
    do i=1,3
       f=(data%coord(i)-genpar%omodel(i))*deltai(i)
       j=f
       ix(i)=1+j
       fx(i)=f-j
       gx(i)=1.-fx(i)
    end do
    
    utmp=               gx(1)*gx(2)*gx(3)*u(ix(1)  ,ix(2)  ,ix(3)  )+&
       &                gx(1)*gx(2)*fx(3)*u(ix(1)  ,ix(2)  ,ix(3)+1)+&
       &                gx(1)*fx(2)*gx(3)*u(ix(1)  ,ix(2)+1,ix(3)  )+&
       &                gx(1)*fx(2)*fx(3)*u(ix(1)  ,ix(2)+1,ix(3)+1)+&
       &                fx(1)*gx(2)*gx(3)*u(ix(1)+1,ix(2)  ,ix(3)  )+&
       &                fx(1)*gx(2)*fx(3)*u(ix(1)+1,ix(2)  ,ix(3)+1)+&
       &                fx(1)*fx(2)*gx(3)*u(ix(1)+1,ix(2)+1,ix(3)  )+&
       &                fx(1)*fx(2)*fx(3)*u(ix(1)+1,ix(2)+1,ix(3)+1)

    if (genpar%rec_type.eq.0) then

       data%trace(it,1)=data%trace(it,1)+vel**2*(utmp)*product(deltai)

    else

       utmpm=           gx(1)*gx(2)*gx(3)*u(-ix(1)  ,ix(2)  ,ix(3)  )+&
       &                gx(1)*gx(2)*fx(3)*u(-ix(1)  ,ix(2)  ,ix(3)+1)+&
       &                gx(1)*fx(2)*gx(3)*u(-ix(1)  ,ix(2)+1,ix(3)  )+&
       &                gx(1)*fx(2)*fx(3)*u(-ix(1)  ,ix(2)+1,ix(3)+1)+&
       &                fx(1)*gx(2)*gx(3)*u(-ix(1)+1,ix(2)  ,ix(3)  )+&
       &                fx(1)*gx(2)*fx(3)*u(-ix(1)+1,ix(2)  ,ix(3)+1)+&
       &                fx(1)*fx(2)*gx(3)*u(-ix(1)+1,ix(2)+1,ix(3)  )+&
       &                fx(1)*fx(2)*fx(3)*u(-ix(1)+1,ix(2)+1,ix(3)+1)

       data%trace(it,1)=data%trace(it,1)+vel**2*(utmp-utmpm)*product(deltai)

    end if    

    deallocate(deltai)

  end subroutine Extraction_1trace_lint3_xyz

  subroutine Extraction_1trace_lint2_xy(bounds,vel,data,u,genpar,it)
    type(FDbounds)    ::                 bounds
    real              ::                        vel
    type(TraceSpace)  ::                            data
    type(GeneralParam)::                                    genpar 
    real              ::                                 u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    integer           :: i,j

    real              :: f,utmp,utmpm
    real              :: fx(2),gx(2)
    integer           :: ix(2)

    allocate(deltai(2))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)
    
    do i=1,2
       f=(data%coord(i)-genpar%omodel(i))*deltai(i)
       j=f
       ix(i)=1+j
       fx(i)=f-j
       gx(i)=1.-fx(i)
    end do
    
    utmp=               gx(1)*gx(2)*u(ix(1)  ,ix(2)  ,1)+&
       &                gx(1)*fx(2)*u(ix(1)  ,ix(2)+1,1)+&
       &                fx(1)*gx(2)*u(ix(1)+1,ix(2)  ,1)+&
       &                fx(1)*fx(2)*u(ix(1)+1,ix(2)+1,1)

    if (genpar%rec_type.eq.0) then

       data%trace(it,1)=data%trace(it,1)+vel**2*utmp*product(deltai)

    else     
         
       utmpm=           gx(1)*gx(2)*u(-ix(1)  ,ix(2)  ,1)+&
       &                gx(1)*fx(2)*u(-ix(1)  ,ix(2)+1,1)+&
       &                fx(1)*gx(2)*u(-ix(1)+1,ix(2)  ,1)+&
       &                fx(1)*fx(2)*u(-ix(1)+1,ix(2)+1,1)

       data%trace(it,1)=data%trace(it,1)+vel**2*(utmp-utmpm)*product(deltai)

    end if    

    deallocate(deltai)

  end subroutine Extraction_1trace_lint2_xy

  subroutine Extraction_1trace_sinc_xy(bounds,vel,data,u,genpar,it)
    type(FDbounds)    ::               bounds
    real              ::                      vel
    type(TraceSpace)  ::                           data
    type(GeneralParam)::                                  genpar 
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j,k
    integer           :: minx,maxx,minz,maxz

    allocate(sinc(genpar%lsinc,2))
    allocate(deltai(2))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)

    minx=-genpar%lsinc*0.5
    minz=minx
    
    maxx=genpar%lsinc*0.5
    maxz=maxx

    do i=1,2
       call mksinc(sinc(:,i),genpar%lsinc,data%dcoord(i)*deltai(i))
    end do
    
    if (genpar%rec_type.eq.0) then
       do j=minx,maxx
          do i=minz,maxz
             data%trace(it,1)=data%trace(it,1)+ &
             &   vel**2* &
             &   u(i+data%icoord(1),j+data%icoord(2),1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*product(deltai) 
          end do
       end do
    else
       do j=minx,maxx
          do i=minz,maxz
             data%trace(it,1)=data%trace(it,1)+ &
             &   vel**2* &
             &   (u(i+data%icoord(1),j+data%icoord(2),1)-u(-i-data%icoord(1),j+data%icoord(2),1))*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*product(deltai) 
          end do
       end do
    end if
    
    deallocate(sinc,deltai)
    
  end subroutine Extraction_1trace_sinc_xy
  
  subroutine Extraction_wavefield(bounds,model,elev,u,genpar,it,dat)
    optional          ::                                        dat
    type(FDbounds)    ::          bounds
    type(ModelSpace)  ::                 model
    type(WaveSpace)   ::                       dat
    type(ModelSpace_elevation) ::                  elev
    type(GeneralParam)::                                  genpar
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i,k,j
    
    real, allocatable :: buffer_sou(:,:)

    allocate(buffer_sou(model%nz,model%nxw))

    MODULO:if (mod(it,genpar%snapi).eq.0) then
       dat%counter=dat%counter+1
       if (genpar%surf_type.ne.0) then
          do k=1,model%nyw
             buffer_sou = 0.
             do j=1,model%nxw
                do i=elev%ielev_z(j,k),model%nz
                   buffer_sou(i,j) = u(i,j,k)
                end do
             end do
             dat%wave(:,:,k,dat%counter,1)=buffer_sou
          end do
       else
          dat%wave(:,:,:,dat%counter,1)=u(1:model%nz,1:model%nxw,1:model%nyw)         
       end if
    end if MODULO
    
    deallocate(buffer_sou)

  end subroutine Extraction_wavefield

  subroutine Extraction_wavefield_copy_to_disk(bounds,model,elev,u,genpar,it,dat)
    optional          ::                                                     dat
    type(FDbounds)    ::                       bounds
    type(ModelSpace)  ::                              model
    type(WaveSpace)   ::                                    dat
    type(ModelSpace_elevation) ::                               elev
    type(GeneralParam)::                                               genpar
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i,k,j
    
    real, allocatable :: buffer_sou(:,:)
    character(len=8)  :: wavetag

    if (genpar%tmin.eq.1) then
       wavetag=model%waFtag
    else
       wavetag=model%waBtag
    end if

    allocate(buffer_sou(model%nz,model%nxw))

    MODULO:if (mod(it,genpar%snapi).eq.0) then
       if (genpar%surf_type.ne.0) then
          do k=1,model%nyw
             buffer_sou = 0.
             do j=1,model%nxw
                do i=elev%ielev_z(j,k),model%nz
                   buffer_sou(i,j) = u(i,j,k)
                end do
             end do
             call srite(wavetag,buffer_sou,4*model%nxw*model%nz)
          end do
       else  
          call srite(wavetag,u(1:model%nz,1:model%nxw,1:model%nyw),4*model%nxw*model%nyw*model%nz)
       end if
    end if MODULO
    
    deallocate(buffer_sou)

  end subroutine Extraction_wavefield_copy_to_disk

end module Extraction_mod
