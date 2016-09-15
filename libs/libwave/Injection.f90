module Injection_mod

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod
  use FD_types

  implicit none

contains

  subroutine Injection_lint(bounds,model,tracevec,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(TraceSpace),dimension(:)  ::    tracevec
    type(GeneralParam)::                            genpar 
    real              ::                          u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    integer           :: type_inject
    real              :: v2r

    if (genpar%tmax.eq.1) then
       type_inject=genpar%rec_type
    else
       type_inject=genpar%shot_type
    end if

    if (genpar%twoD) then
       do i=1,size(tracevec)
          if (genpar%withRho) then
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),1)**2* &
             &   model%rho(tracevec(i)%icoord(1),tracevec(i)%icoord(2),1)
          else
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),1)**2
          end if
          call Injection_source_lint3_xz(bounds,v2r,tracevec(i),u,genpar,it,type_inject)
       end do
    else
       !$OMP PARALLEL DO PRIVATE(i,v2r)
       do i=1,size(tracevec)
          if (genpar%withRho) then
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))**2* &
             &   model%rho(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))
          else
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))**2
          end if
          call Injection_source_lint3_xyz(bounds,v2r,tracevec(i),u,genpar,it,type_inject) 
       end do    
       !$OMP END PARALLEL DO
    end if
    
  end subroutine Injection_lint

  subroutine Injection_sinc(bounds,model,tracevec,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(TraceSpace),dimension(:)  ::    tracevec
    type(GeneralParam)::                            genpar 
    real              ::                          u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    integer           :: type_inject
    real              :: v2r

    if (genpar%tmax.eq.1) then
       type_inject=genpar%rec_type
    else
       type_inject=genpar%shot_type
    end if

    if (genpar%twoD) then
       do i=1,size(tracevec)
          if (genpar%withRho) then
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))**2* &
             &   model%rho(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))
          else
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))**2
          end if
          call Injection_source_sinc_xz(bounds,v2r,tracevec(i),u,genpar,it,type_inject)
       end do
    else
       !$OMP PARALLEL DO PRIVATE(i,v2r)
       do i=1,size(tracevec)
          if (genpar%withRho) then
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))**2* &
             &   model%rho(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))
          else
             v2r=model%vel(tracevec(i)%icoord(1),tracevec(i)%icoord(2),tracevec(i)%icoord(3))**2
          end if
          call Injection_source_sinc_xyz(bounds,v2r,tracevec(i),u,genpar,it,type_inject) 
       end do    
       !$OMP END PARALLEL DO
    end if
    
  end subroutine Injection_sinc

  subroutine Injection_LSRTM_sinc(bounds,model,tracevec,u,genpar,it)
    
    type(FDbounds)    ::    bounds
    type(ModelSpace)  ::           model
    type(TraceSpace),dimension(:)  ::    tracevec
    type(GeneralParam)::                            genpar 
    real              ::                          u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: i
    integer           :: type_inject
    real              :: v2r

    if (genpar%tmax.eq.1) then
       type_inject=genpar%rec_type
    else
       type_inject=genpar%shot_type
    end if

    if (genpar%twoD) then
       do i=1,size(tracevec)
          v2r=1/product(genpar%delta(1:2))
          call Injection_source_sinc_xz(bounds,v2r,tracevec(i),u,genpar,it,type_inject)
       end do
    else
       !$OMP PARALLEL DO PRIVATE(i,v2r)
       do i=1,size(tracevec)
          v2r=1/product(genpar%delta(1:3))
          call Injection_source_sinc_xyz(bounds,v2r,tracevec(i),u,genpar,it,type_inject) 
       end do
       !$OMP END PARALLEL DO
    end if
    
  end subroutine Injection_LSRTM_sinc

  subroutine Injection_source_sinc_xyz(bounds,v2r,sou,u,genpar,it,type_inject)
    type(FDbounds)    ::               bounds
    real              ::                      v2r
    type(TraceSpace)  ::                          sou
    type(GeneralParam)::                                genpar 
    real              ::                              u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: type_inject

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j,k
    integer           :: minx,maxx,miny,maxy,minz,maxz

    allocate(sinc(genpar%lsinc,3))
    allocate(deltai(3))
    deltai(1)=genpar%delta(1)
    deltai(2)=genpar%delta(2)
    deltai(3)=genpar%delta(3)

    minx=-genpar%lsinc*0.5
    minz=minx
    miny=minx
    
    maxx=genpar%lsinc*0.5
    maxz=maxx
    maxy=maxx

    do i=1,3
       call mksinc(sinc(:,i),genpar%lsinc,sou%dcoord(i)/deltai(i))
    end do
    
    if (type_inject.eq.0) then
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=miny,maxy
          do j=minx,maxx
             do i=minz,maxz
                u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i+sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))+&
                &   v2r*sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
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
                &   v2r*sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai) 

                u(i-sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))=&
                &   u(i-sou%icoord(1),j+sou%icoord(2),k+sou%icoord(3))-&
                &   v2r*sou%trace(it,1)*sinc(maxz+1-i,1)*sinc(maxx+1+j,2)*sinc(maxy+1+k,2)*&
                &   product(deltai)            
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    deallocate(sinc,deltai)

  end subroutine Injection_source_sinc_xyz

  subroutine Injection_source_sinc_xz(bounds,v2r,sou,u,genpar,it,type_inject)
    type(FDbounds)    ::              bounds
    real              ::                     v2r
    type(TraceSpace)  ::                         sou
    type(GeneralParam)::                               genpar 
    real              ::                             u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: type_inject

    real, allocatable :: deltai(:)
    real, allocatable :: sinc(:,:)
    integer           :: i,j
    integer           :: minx,maxx,minz,maxz

    allocate(sinc(genpar%lsinc,2))
    allocate(deltai(3))
    deltai(1)=genpar%delta(1)
    deltai(2)=genpar%delta(2)

    minx=-genpar%lsinc*0.5
    minz=minx
    
    maxx=genpar%lsinc*0.5
    maxz=maxx

    do i=1,2
       call mksinc(sinc(:,i),genpar%lsinc,sou%dcoord(i)/deltai(i))
    end do
    
    if (genpar%shot_type.eq.0) then
       do j=minx,maxx
          do i=minz,maxz            
             u(i+sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i+sou%icoord(1),j+sou%icoord(2),1)+&
             &   v2r*sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)
             
          end do
       end do
    else
       do j=minx,maxx
          do i=minz,maxz            
             u(i+sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i+sou%icoord(1),j+sou%icoord(2),1)+&
             &   v2r*sou%trace(it,1)*sinc(maxz+1+i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)          
             u(i-sou%icoord(1),j+sou%icoord(2),1)=&
             &   u(i-sou%icoord(1),j+sou%icoord(2),1)-&
             &   v2r*sou%trace(it,1)*sinc(maxz+1-i,1)*sinc(maxx+1+j,2)*&
             &   deltai(1)*deltai(2)
          end do
       end do
    end if

    deallocate(sinc,deltai)

  end subroutine Injection_source_sinc_xz

  subroutine Injection_source_lint3_xyz(bounds,v2r,sou,u,genpar,it,type_inject)
    type(FDbounds)    ::                bounds
    real              ::                       v2r
    type(TraceSpace)  ::                           sou
    type(GeneralParam)::                                 genpar 
    real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: type_inject

    real, allocatable :: deltai(:)
    integer           :: i,j,k

    real              :: f,utmp,utmpm
    real              :: fx(3),gx(3)
    integer           :: ix(3)

    allocate(deltai(3))
    deltai(1)=genpar%delta(1)
    deltai(2)=genpar%delta(2)
    deltai(3)=genpar%delta(3)

    do i=1,3
       f=(sou%coord(i)-genpar%omodel(i))/deltai(i)
       j=f
       ix(i)=1+j
       fx(i)=f-j
       gx(i)=1.-fx(i)
    end do

    u(ix(1)  ,ix(2)  ,ix(3)  )= u(ix(1)  ,ix(2)  ,ix(3)  ) + gx(1)*gx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)  ,ix(2)  ,ix(3)+1)= u(ix(1)  ,ix(2)  ,ix(3)+1) + gx(1)*gx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)  ,ix(2)+1,ix(3)  )= u(ix(1)  ,ix(2)+1,ix(3)  ) + gx(1)*fx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)  ,ix(2)+1,ix(3)+1)= u(ix(1)  ,ix(2)+1,ix(3)+1) + gx(1)*fx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)+1,ix(2)  ,ix(3)  )= u(ix(1)+1,ix(2)  ,ix(3)  ) + fx(1)*gx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)+1,ix(2)  ,ix(3)+1)= u(ix(1)+1,ix(2)  ,ix(3)+1) + fx(1)*gx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)+1,ix(2)+1,ix(3)  )= u(ix(1)+1,ix(2)+1,ix(3)  ) + fx(1)*fx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)+1,ix(2)+1,ix(3)+1)= u(ix(1)+1,ix(2)+1,ix(3)+1) + fx(1)*fx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)

    if (type_inject.eq.0) then
       
       u(-ix(1)  ,ix(2)  ,ix(3)  )= u(-ix(1)  ,ix(2)  ,ix(3)  ) + gx(1)*gx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)  ,ix(2)  ,ix(3)+1)= u(-ix(1)  ,ix(2)  ,ix(3)+1) + gx(1)*gx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)  ,ix(2)+1,ix(3)  )= u(-ix(1)  ,ix(2)+1,ix(3)  ) + gx(1)*fx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)  ,ix(2)+1,ix(3)+1)= u(-ix(1)  ,ix(2)+1,ix(3)+1) + gx(1)*fx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)+1,ix(2)  ,ix(3)  )= u(-ix(1)+1,ix(2)  ,ix(3)  ) + fx(1)*gx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)+1,ix(2)  ,ix(3)+1)= u(-ix(1)+1,ix(2)  ,ix(3)+1) + fx(1)*gx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)+1,ix(2)+1,ix(3)  )= u(-ix(1)+1,ix(2)+1,ix(3)  ) + fx(1)*fx(2)*gx(3)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)+1,ix(2)+1,ix(3)+1)= u(-ix(1)+1,ix(2)+1,ix(3)+1) + fx(1)*fx(2)*fx(3)*v2r*sou%trace(it,1)*product(deltai)
       
    end if

    deallocate(deltai)

  end subroutine Injection_source_lint3_xyz

  subroutine Injection_source_lint3_xz(bounds,v2r,sou,u,genpar,it,type_inject)
    type(FDbounds)    ::               bounds
    real              ::                      v2r
    type(TraceSpace)  ::                          sou
    type(GeneralParam)::                                genpar 
    real              ::                              u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    integer           :: type_inject

    real, allocatable :: deltai(:)
    integer           :: i,j,k

    real              :: f,utmp,utmpm
    real              :: fx(2),gx(2)
    integer           :: ix(2)

    allocate(deltai(2))
    deltai(1)=genpar%delta(1)
    deltai(2)=genpar%delta(2)

    do i=1,2
       f=(sou%coord(i)-genpar%omodel(i))/deltai(i)
       j=f
       ix(i)=1+j
       fx(i)=f-j
       gx(i)=1.-fx(i)
    end do

    u(ix(1)  ,ix(2)  ,1)= u(ix(1)  ,ix(2)  ,1) + gx(1)*gx(2)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)  ,ix(2)+1,1)= u(ix(1)  ,ix(2)+1,1) + gx(1)*fx(2)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)+1,ix(2)  ,1)= u(ix(1)+1,ix(2)  ,1) + fx(1)*gx(2)*v2r*sou%trace(it,1)*product(deltai)
    u(ix(1)+1,ix(2)+1,1)= u(ix(1)+1,ix(2)+1,1) + fx(1)*fx(2)*v2r*sou%trace(it,1)*product(deltai)

    if (type_inject.eq.0) then
       
       u(-ix(1)  ,ix(2)  ,1)= u(-ix(1)  ,ix(2)  ,1) + gx(1)*gx(2)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)  ,ix(2)+1,1)= u(-ix(1)  ,ix(2)+1,1) + gx(1)*fx(2)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)+1,ix(2)  ,1)= u(-ix(1)+1,ix(2)  ,1) + fx(1)*gx(2)*v2r*sou%trace(it,1)*product(deltai)
       u(-ix(1)+1,ix(2)+1,1)= u(-ix(1)+1,ix(2)+1,1) + fx(1)*fx(2)*v2r*sou%trace(it,1)*product(deltai)
       
    end if

    deallocate(deltai)

  end subroutine Injection_source_lint3_xz

  subroutine Injection_source_simple_xyz(bounds,model,sou,u,genpar,it,type_inject)
    type(FDbounds)    ::                 bounds
    type(ModelSpace)  ::                        model
    type(TraceSpace)  ::                            sou
    type(GeneralParam)::                                  genpar 
    real              ::                                u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    real              :: deltai
    integer           :: i,j,k,iz,ix,iy
    integer           :: type_inject

    iz=sou%icoord(1)
    ix=sou%icoord(2)

    if(genpar%twoD) then
       iy=1
       deltai=(genpar%delta(1)*genpar%delta(2))
    else
       iy=sou%icoord(3)
       deltai=(genpar%delta(1)*genpar%delta(2)*genpar%delta(3))
    end if

    if (type_inject.eq.0) then
       u( iz,ix,iy)=u( iz,ix,iy)+model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai  
    else
       u( iz,ix,iy)=u( iz,ix,iy)+model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai      
       u(-iz,ix,iy)=u(-iz,ix,iy)-model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai
    end if

  end subroutine Injection_source_simple_xyz

  subroutine Injection_source_rho_simple_xyz(bounds,model,sou,u,genpar,it,type_inject)
    type(FDbounds)    ::                 bounds
    type(ModelSpace)  ::                        model
    type(TraceSpace)  ::                              sou
    type(GeneralParam)::                                    genpar 
    real              ::                                  u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
    &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer           :: it
    real              :: deltai
    integer           :: i,j,k,iz,ix,iy
    integer           :: type_inject

    iz=sou%icoord(1)
    ix=sou%icoord(2)

    if(genpar%twoD) then
       iy=1
       deltai=(genpar%delta(1)*genpar%delta(2))
    else
       iy=sou%icoord(3)
       deltai=(genpar%delta(1)*genpar%delta(2)*genpar%delta(3))
    end if

    if (type_inject.eq.0) then
       u( iz,ix,iy)=u( iz,ix,iy)+model%rho(iz,ix,iy)*model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai  
    else
       u( iz,ix,iy)=u( iz,ix,iy)+model%rho(iz,ix,iy)*model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai      
       u(-iz,ix,iy)=u(-iz,ix,iy)-model%rho(iz,ix,iy)*model%vel(iz,ix,iy)**2*sou%trace(it,1)*deltai
    end if

  end subroutine Injection_source_rho_simple_xyz

end module Injection_mod
