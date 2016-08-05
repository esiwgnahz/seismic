module Injection_mod

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod

  implicit none

contains

  subroutine Injection_source(bounds,mod,dat,u)
    type(FDbounds)    ::      bounds
    type(ModelSpace)  ::             mod
    type(DataSpace)   ::                 dat
    real              ::                     u(bounds%nmin1-4:bounds%nmax1+4)
    type(GeneralParam):: genpar
    real              :: dxi,dzi,dyi
    real, allocatable :: sinc(:)

    allocate(sinc(genpar%lsinc))

    dxi=1./mod%dx
    dzi=1./mod%dz
    dyi=1./mod%dy
    
    if (mod%shot_type.eq.0) then
       if (mod%dshot_z.eq.0) then
          u(mod%ishot_z) = u(mod%ishot_z)+ &
          &  mod%vel(mod%ishot_z,mod%ishot_x,mod%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dzi*dyi
       else
          call mksinc(sinc,genpar%lsinc,dshot_z*dzi)
          do i=-genpar%lsinc*0.5,genpar%lsinc*0.5
             u(i+mod%ishot_z) = u(i+mod%ishot_z)+ &
             &  mod%vel(mod%ishot_z,mod%ishot_x,mod%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+1+i)* &
             &  dxi*dzi*dyi
          end do
       end if
    else
       if (mod%dshot_z.eq.0.) then
          u(mod%ishot_z) = &
          &  u(mod%ishot_z)+ &
          &  mod%vel(mod%ishot_z,mod%ishot_x,mod%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dyi*dzi
          u(mod%ishot_z-2) = &
          &  u(mod%ishot_z-2)- &
          &  mod%vel(mod%ishot_z,ishot_x,mod%ishot_y)**2*dat%source(dat%it)* &
          &  dxi*dyi*dzi
       else
          call mksinc(sinc,genpar%lsinc,dshot_z*dzi)
          do i=-genpar%lsinc*0.5,genpar%lsinc*0.5
             u( i+mod%ishot_z) = u( i+mod%ishot_z)+ &
             &  mod%vel(mod%ishot_z,mod%ishot_x,mod%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+i+1)* &
             &  dxi*dyi*dzi
             u(-i+mod%ishot_z-2) = u(-i+mod%ishot_z-2)- &
             &  mod%vel(mod%ishot_z,mod%ishot_x,mod%ishot_y)**2* &
             &  dat%source(dat%it)*sinc(genpar%lsinc*0.5+1-i)* &
             &  dxi*dyi*dzi
          end do
       endif
    end if

    deallocate(sinc)

  end subroutine Injection_source

end module Injection_mod
