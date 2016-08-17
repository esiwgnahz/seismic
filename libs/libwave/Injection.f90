module Injection_mod

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod
  use FD_types

  implicit none

contains

  subroutine Injection_sinc(bounds,model,tracevec,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(TraceSpace),dimension(:)  ::    tracevec
    type(GeneralParam)::                            genpar 
    real              ::                          u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    
    if (genpar%twoD) then
       do i=1,size(tracevec)
          call Injection_source_sinc_xz(bounds,model,tracevec(i),u,genpar,it)
       end do
    else
       do i=1,size(tracevec)
          call Injection_source_sinc_xyz(bounds,model,tracevec(i),u,genpar,it) 
       end do
    end if
    
  end subroutine Injection_sinc

  subroutine Injection_rho_sinc(bounds,model,tracevec,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(TraceSpace),dimension(:)  ::    tracevec
    type(GeneralParam)::                            genpar 
    real              ::                          u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    
    if (genpar%twoD) then
       do i=1,size(tracevec)
          call Injection_source_rho_sinc_xz(bounds,model,tracevec(i),u,genpar,it)
       end do
    else
       do i=1,size(tracevec)
          call Injection_source_rho_sinc_xyz(bounds,model,tracevec(i),u,genpar,it) 
       end do
    end if
    
  end subroutine Injection_rho_sinc

  subroutine Injection_source_sinc_xyz(bounds,model,sou,u,genpar,it)
    type(FDbounds)    ::               bounds
    type(ModelSpace)  ::                      model
    type(TraceSpace)  ::                            sou
    type(GeneralParam)::                                  genpar 
    real              ::                                u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
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
       call mksinc(sinc(:,i),genpar%lsinc,sou%dcoord(i)*deltai(i))
    end do
    
    if (genpar%shot_type.eq.0) then
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))+&
                &   model%vel(sou%icoord(1),sou%icoord(2),sou%icoord(3))**2* &
                &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai)            
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))+&
                &   model%vel(sou%icoord(1),sou%icoord(2),sou%icoord(3))**2* &
                &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai) 

                u(i-sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i-sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))-&
                &   model%vel(sou%icoord(1),sou%icoord(2),sou%icoord(3))**2* &
                &   sou%trace(it,1)*sinc(maxz+1-i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai)            
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    deallocate(sinc,deltai)

  end subroutine Injection_source_sinc_xyz

  subroutine Injection_source_rho_sinc_xyz(bounds,model,sou,u,genpar,it)
    type(FDbounds)    ::                   bounds
    type(ModelSpace)  ::                          model
    type(TraceSpace)  ::                                sou
    type(GeneralParam)::                                      genpar 
    real              ::                                   u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
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
       call mksinc(sinc(:,i),genpar%lsinc,sou%dcoord(i)*deltai(i))
    end do
    
    if (genpar%shot_type.eq.0) then
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))+&
                &   model%rho(sou%icoord(1),sou%icoord(2),sou%icoord(3))*&
                &   model%vel(sou%icoord(1),sou%icoord(2),sou%icoord(3))**2* &
                &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai)            
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))+&
                &   model%rho(sou%icoord(1),sou%icoord(2),sou%icoord(3))*&
                &   model%vel(sou%icoord(1),sou%icoord(2),sou%icoord(3))**2* &
                &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai) 

                u(i-sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i-sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))-&
                &   model%rho(sou%icoord(1),sou%icoord(2),sou%icoord(3))*&
                &   model%vel(sou%icoord(1),sou%icoord(2),sou%icoord(3))**2* &
                &   sou%trace(it,1)*sinc(maxz+1-i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai)            
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    deallocate(sinc,deltai)

  end subroutine Injection_source_rho_sinc_xyz

  subroutine Injection_source_sinc_xz(bounds,model,sou,u,genpar,it)
    type(FDbounds)    ::              bounds
    type(ModelSpace)  ::                     model
    type(TraceSpace)  ::                           sou
    type(GeneralParam)::                                 genpar 
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j
    integer           :: minx,maxx,minz,maxz

    allocate(sinc(genpar%lsinc,2))
    allocate(deltai(3))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)

    minx=-genpar%lsinc*0.5
    minz=minx
    
    maxx=genpar%lsinc*0.5
    maxz=maxx

    do i=1,2
       call mksinc(sinc(:,i),genpar%lsinc,sou%dcoord(i)*deltai(i))
    end do
    
    if (genpar%shot_type.eq.0) then
       do j=minx,maxx
          do i=minz,maxz            
             u(i+sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i+sou%icoord(1),j+sou%icoord(2),1)+&
             &   model%vel(sou%icoord(1),sou%icoord(2),1)**2* &
             &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)
             
          end do
       end do
    else
       do j=minx,maxx
          do i=minz,maxz            
             u(i+sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i+sou%icoord(1),j+sou%icoord(2),1)+&
             &   model%vel(sou%icoord(1),sou%icoord(2),1)**2* &
             &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)          
             u(i-sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i-sou%icoord(1),j+sou%icoord(2),1)-&
             &   model%vel(sou%icoord(1),sou%icoord(2),1)**2* &
             &   sou%trace(it,1)*sinc(maxz+1-i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)
          end do
       end do
    end if

    deallocate(sinc,deltai)

  end subroutine Injection_source_sinc_xz

  subroutine Injection_source_rho_sinc_xz(bounds,model,sou,u,genpar,it)
    type(FDbounds)    ::              bounds
    type(ModelSpace)  ::                     model
    type(TraceSpace)  ::                           sou
    type(GeneralParam)::                                 genpar 
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j
    integer           :: minx,maxx,minz,maxz

    allocate(sinc(genpar%lsinc,2))
    allocate(deltai(3))
    deltai(1)=1./genpar%delta(1)
    deltai(2)=1./genpar%delta(2)

    minx=-genpar%lsinc*0.5
    minz=minx
    
    maxx=genpar%lsinc*0.5
    maxz=maxx

    do i=1,2
       call mksinc(sinc(:,i),genpar%lsinc,sou%dcoord(i)*deltai(i))
    end do
    
    if (genpar%shot_type.eq.0) then
       do j=minx,maxx
          do i=minz,maxz            
             u(i+sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i+sou%icoord(1),j+sou%icoord(2),1)+&
             &   model%rho(sou%icoord(1),sou%icoord(2),1)*&
             &   model%vel(sou%icoord(1),sou%icoord(2),1)**2* &
             &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)
             
          end do
       end do
    else
       do j=minx,maxx
          do i=minz,maxz            
             u(i+sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i+sou%icoord(1),j+sou%icoord(2),1)+&
             &   model%rho(sou%icoord(1),sou%icoord(2),1)*&
             &   model%vel(sou%icoord(1),sou%icoord(2),1)**2* &
             &   sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)          
             u(i-sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i-sou%icoord(1),j+sou%icoord(2),1)-&
             &   model%rho(sou%icoord(1),sou%icoord(2),1)*&
             &   model%vel(sou%icoord(1),sou%icoord(2),1)**2* &
             &   sou%trace(it,1)*sinc(maxz+1-i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)
          end do
       end do
    end if

    deallocate(sinc,deltai)

  end subroutine Injection_source_rho_sinc_xz

  subroutine Injection_source_simple_xyz(bounds,model,sou,u,genpar,it)
    type(FDbounds)    ::                 bounds
    type(ModelSpace)  ::                        model
    type(TraceSpace)  ::                              sou
    type(GeneralParam)::                                    genpar 
    real              ::                                  u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    real              :: deltai
    integer           :: i,j,k,iz,ix,iy

    iz=sou%icoord(1)
    ix=sou%icoord(2)

    if(genpar%twoD) then
       iy=1
       deltai=1./(genpar%delta(1)*genpar%delta(2))
    else
       iy=sou%icoord(3)
       deltai=1./(genpar%delta(1)*genpar%delta(2)*genpar%delta(3))
    end if

    if (genpar%shot_type.eq.0) then
       u( iz,ix,iy)=u( iz,ix,iy)+model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai  
    else
       u( iz,ix,iy)=u( iz,ix,iy)+model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai      
       u(-iz,ix,iy)=u(-iz,ix,iy)-model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai
    end if

  end subroutine Injection_source_simple_xyz

  subroutine Injection_source_rho_simple_xyz(bounds,model,sou,u,genpar,it)
    type(FDbounds)    ::                 bounds
    type(ModelSpace)  ::                        model
    type(TraceSpace)  ::                              sou
    type(GeneralParam)::                                    genpar 
    real              ::                                  u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    real              :: deltai
    integer           :: i,j,k,iz,ix,iy

    iz=sou%icoord(1)
    ix=sou%icoord(2)

    if(genpar%twoD) then
       iy=1
       deltai=1./(genpar%delta(1)*genpar%delta(2))
    else
       iy=sou%icoord(3)
       deltai=1./(genpar%delta(1)*genpar%delta(2)*genpar%delta(3))
    end if

    if (genpar%shot_type.eq.0) then
       u( iz,ix,iy)=u( iz,ix,iy)+model%rho(iz,ix,iy)*model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai  
    else
       u( iz,ix,iy)=u( iz,ix,iy)+model%rho(iz,ix,iy)*model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai      
       u(-iz,ix,iy)=u(-iz,ix,iy)-model%rho(iz,ix,iy)*model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai
    end if

  end subroutine Injection_source_rho_simple_xyz

  subroutine Injection_Born(bounds,model,tracevec,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(TraceSpace),dimension(:)  ::    tracevec
    type(GeneralParam)::                            genpar 
    real              ::                          u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it   

    if (mod(it,genpar%snapi).eq.0) then
       model%wvfld%counter=model%wvfld%counter+1
       u(1:model%nz,1:model%nx,1:model%ny)=u(1:model%nz,1:model%nx,1:model%ny)+2*dprod(model%image(1:model%nz,1:model%nx,1:model%ny),model%wvfld%wave(1:model%nz,1:model%nx,1:model%ny,model%wvfld%counter,1))/model%vel(1:model%nz,1:model%nx,1:model%ny)
    end if
    
  end subroutine Injection_Born

end module Injection_mod
