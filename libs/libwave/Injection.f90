module Injection_mod

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod
  use FD_types

  implicit none

contains

  subroutine Injection_source(bounds,model,dat,elev,u,genpar)
    type(FDbounds)    ::      bounds
    type(ModelSpace)  ::             model
    type(DataSpace)   ::                   dat
    type(ModelSpace_elevation)   ::            elev
    real              ::                            u(bounds%nmin1-4:bounds%nmax1+4)
    type(GeneralParam):: genpar 
    real              :: dxi,dzi,dyi
    real, allocatable :: sinc(:)
    integer           :: i

    allocate(sinc(genpar%lsinc))

    dxi=1./genpar%dx
    dzi=1./genpar%dz
    dyi=1./genpar%dy
    
    if (genpar%shot_type.eq.0) then
       if (elev%dshot_z.eq.0) then
          u(elev%ishot_z) = u(elev%ishot_z)+ &
          &  model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dzi*dyi
       else
          call mksinc(sinc,genpar%lsinc,elev%dshot_z*dzi)
          do i=-genpar%lsinc*0.5,genpar%lsinc*0.5
             u(i+elev%ishot_z) = u(i+elev%ishot_z)+ &
             &  model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+1+i)* &
             &  dxi*dzi*dyi
          end do
       end if
    else
       if (elev%dshot_z.eq.0.) then
          u(elev%ishot_z) = &
          &  u(elev%ishot_z)+ &
          &  model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dyi*dzi
          u(elev%ishot_z-2) = &
          &  u(elev%ishot_z-2)- &
          &  model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dyi*dzi
       else
          call mksinc(sinc,genpar%lsinc,elev%dshot_z*dzi)
          do i=-genpar%lsinc*0.5,genpar%lsinc*0.5
             u( i+elev%ishot_z) = u( i+elev%ishot_z)+ &
             &  model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+i+1)* &
             &  dxi*dyi*dzi
             u(-i+elev%ishot_z-2) = u(-i+elev%ishot_z-2)- &
             &  model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+1-i)* &
             &  dxi*dyi*dzi
          end do
       endif
    end if
    deallocate(sinc)

  end subroutine Injection_source

  subroutine Injection_source_rho(bounds,model,dat,elev,u,genpar)
    type(FDbounds)    ::          bounds
    type(ModelSpace)  ::                 model
    type(DataSpace)   ::                       dat
    type(ModelSpace_elevation)   ::                elev
    real              ::                                u(bounds%nmin1-4:bounds%nmax1+4)
    type(GeneralParam):: genpar 
    real              :: dxi,dzi,dyi
    real, allocatable :: sinc(:)
    integer           :: i

    allocate(sinc(genpar%lsinc))

    dxi=1./genpar%dx
    dzi=1./genpar%dz
    dyi=1./genpar%dy
    
    if (genpar%shot_type.eq.0) then
       if (elev%dshot_z.eq.0) then
          u(elev%ishot_z) = u(elev%ishot_z)+ &
          &  model%rho(elev%ishot_z,elev%ishot_x,elev%ishot_y)*model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dzi*dyi
       else
          call mksinc(sinc,genpar%lsinc,elev%dshot_z*dzi)
          do i=-genpar%lsinc*0.5,genpar%lsinc*0.5
             u(i+elev%ishot_z) = u(i+elev%ishot_z)+ &
             &  model%rho(elev%ishot_z,elev%ishot_x,elev%ishot_y)*model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+1+i)* &
             &  dxi*dzi*dyi
          end do
       end if
    else
       if (elev%dshot_z.eq.0.) then
          u(elev%ishot_z) = &
          &  u(elev%ishot_z)+ &
          &  model%rho(elev%ishot_z,elev%ishot_x,elev%ishot_y)*model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dyi*dzi
          u(elev%ishot_z-2) = &
          &  u(elev%ishot_z-2)- &
          &  model%rho(elev%ishot_z,elev%ishot_x,elev%ishot_y)*model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dyi*dzi
       else
          call mksinc(sinc,genpar%lsinc,elev%dshot_z*dzi)
          do i=-genpar%lsinc*0.5,genpar%lsinc*0.5
             u( i+elev%ishot_z) = u( i+elev%ishot_z)+ &
             &  model%rho(elev%ishot_z,elev%ishot_x,elev%ishot_y)*model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+i+1)* &
             &  dxi*dyi*dzi
             u(-i+elev%ishot_z-2) = u(-i+elev%ishot_z-2)- &
             &  model%rho(elev%ishot_z,elev%ishot_x,elev%ishot_y)*model%vel(elev%ishot_z,elev%ishot_x,elev%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+1-i)* &
             &  dxi*dyi*dzi
          end do
       endif
    end if

    deallocate(sinc)

  end subroutine Injection_source_rho

end module Injection_mod
