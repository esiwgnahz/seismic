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
    
    if (genpar%twoD) then
       do i=1,size(data)
          call Extraction_1trace_simple_xy(bounds,model,data(i),u,genpar,it)
       end do
    else
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)
          call Extraction_1trace_simple_xyz(bounds,model,data(i),u,genpar,it)
       end do  
       !$OMP END PARALLEL DO    
    end if

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
          call Extraction_1trace_sinc_xy(bounds,model,data(i),u,genpar,it)
       end do
    else      
       !$OMP PARALLEL DO PRIVATE(i)
       do i=1,size(data)       
          call Extraction_1trace_sinc_xyz(bounds,model,data(i),u,genpar,it)
       end do       
       !$OMP END PARALLEL DO
    end if

  end subroutine Extraction_array_sinc

  subroutine Extraction_1trace_sinc_xyz(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::                bounds
    type(ModelSpace)  ::                       model
    type(TraceSpace)  ::                             data
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
    deltai(1)=1./genpar%dz
    deltai(2)=1./genpar%dx
    deltai(3)=1./genpar%dy

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
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                data%trace(it,1)=data%trace(it,1)+ &
                &   model%vel(data%icoord(1),data%icoord(2),data%icoord(3))**2* &
                &   u(i+data%icoord(1),j+data%icoord(2),k+data%icoord(3))*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*product(deltai) 
             end do
          end do
       end do
    
    else
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                data%trace(it,1)=data%trace(it,1)+ &
                &   model%vel(data%icoord(1),data%icoord(2),data%icoord(3))**2* &
                &   (u(i+data%icoord(1),j+data%icoord(2),k+data%icoord(3))-u(i-data%icoord(1),j+data%icoord(2),k+data%icoord(3)))*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*product(deltai) 
             end do
          end do
       end do
    end if

    deallocate(sinc,deltai)

  end subroutine Extraction_1trace_sinc_xyz

  subroutine Extraction_1trace_sinc_xy(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::               bounds
    type(ModelSpace)  ::                      model
    type(TraceSpace)  ::                            data
    type(GeneralParam)::                                   genpar 
    real              ::                                 u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j,k
    integer           :: minx,maxx,minz,maxz

    allocate(sinc(genpar%lsinc,2))
    allocate(deltai(2))
    deltai(1)=1./genpar%dz
    deltai(2)=1./genpar%dx

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
             &   model%vel(data%icoord(1),data%icoord(2),1)**2* &
             &   u(i+data%icoord(1),j+data%icoord(2),1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*product(deltai) 
          end do
       end do
    else
       do j=minx,maxx
          do i=minz,maxz
             data%trace(it,1)=data%trace(it,1)+ &
             &   model%vel(data%icoord(1),data%icoord(2),1)**2* &
             &   (u(i+data%icoord(1),j+data%icoord(2),1)-u(i-data%icoord(1),j+data%icoord(2),1))*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*product(deltai) 
          end do
       end do
    end if
    
    deallocate(sinc,deltai)
    
  end subroutine Extraction_1trace_sinc_xy
  
  subroutine Extraction_1trace_simple_xyz(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::                  bounds
    type(ModelSpace)  ::                         model
    type(TraceSpace)  ::                               data
    type(GeneralParam)::                                      genpar 
    real              ::                                    u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    integer           :: i,j,k

    allocate(deltai(3))
    deltai(1)=1./genpar%dz
    deltai(2)=1./genpar%dx
    deltai(3)=1./genpar%dy

    
    if (genpar%rec_type.eq.0) then
       data%trace(it,1)=data%trace(it,1)+ &
       &   model%vel(data%icoord(1),data%icoord(2),data%icoord(3))**2* &
       &   u(data%icoord(1),data%icoord(2),data%icoord(3))*product(deltai) 
    else     
       data%trace(it,1)=data%trace(it,1)+ &
       &   model%vel(data%icoord(1),data%icoord(2),data%icoord(3))**2* &
       &   (u(data%icoord(1),data%icoord(2),data%icoord(3))-u(-data%icoord(1),data%icoord(2),data%icoord(3)))*product(deltai)
    end if

    deallocate(deltai)

  end subroutine Extraction_1trace_simple_xyz

  subroutine Extraction_1trace_simple_xy(bounds,model,data,u,genpar,it)
    type(FDbounds)    ::                 bounds
    type(ModelSpace)  ::                        model
    type(TraceSpace)  ::                              data
    type(GeneralParam)::                                     genpar 
    real              ::                                   u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    integer           :: i,j

    allocate(deltai(2))
    deltai(1)=1./genpar%dz
    deltai(2)=1./genpar%dx

    
    if (genpar%rec_type.eq.0) then
       data%trace(it,1)=data%trace(it,1)+ &
       &   model%vel(data%icoord(1),data%icoord(2),1)**2* &
       &   u(data%icoord(1),data%icoord(2),1)*product(deltai) 
    else     
       data%trace(it,1)=data%trace(it,1)+ &
       &   model%vel(data%icoord(1),data%icoord(2),1)**2* &
       &   (u(data%icoord(1),data%icoord(2),1)-u(-data%icoord(1),data%icoord(2),1))*product(deltai)
    end if

    deallocate(deltai)

  end subroutine Extraction_1trace_simple_xy

  subroutine Extraction_wavefield(bounds,model,dat,elev,u,genpar,counter,it)
    type(FDbounds)    ::          bounds
    type(ModelSpace)  ::                 model
    type(WaveSpace)   ::                       dat
    type(ModelSpace_elevation) ::                  elev
    type(GeneralParam)::                                  genpar
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
    integer           :: it
    integer           :: i,k,j,counter
    
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
             dat%wave(:,:,k,counter,1)=buffer_sou
          end do
       else
          dat%wave(:,:,:,counter,1)=u(1:model%nz,1:model%nx,1:model%ny,2)         
       end if
    end if MODULO
    
    deallocate(buffer_sou)

  end subroutine Extraction_wavefield

end module Extraction_mod
