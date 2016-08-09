module Propagator_mod
  
  ! Types
  use FD_types
  use ModelSpace_types
  use DataSpace_types
  use GeneralParam_types

  ! Computing
  use Boundary_mod
  use Injection_mod
  use FDcoefs_assign
  use Extraction_mod
  use FD_derivatives
  use FDswaptime_mod

  implicit none

contains

!  subroutine propagator_acoustic(FD_coefs,FD_scheme,bounds,mod,elev,genpar,nbound)
  subroutine propagator_acoustic(FD_coefs,FD_scheme,InjSrc,TimeDer,TimeSwap,bounds,model,dat,elev,genpar)

    interface
       subroutine FD_coefs      (coef)
         use FD_types
         type(UnscaledFDcoefs)   ::      coef        
       end subroutine FD_coefs

       subroutine FD_scheme     (genpar,bounds,u,model)
         use ModelSpace_types
         use GeneralParam_types
         use FD_types
         type(GeneralParam)   :: genpar
         type(ModelSpace)     ::                 model
         type(FDbounds)       ::        bounds
         real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
       end subroutine FD_scheme
      
       subroutine InjSrc(bounds,model,dat,elev,u,genpar)
         use GeneralParam_types
         use ModelSpace_types
         use DataSpace_types
         use FD_types
         type(FDbounds)    ::      bounds
         type(ModelSpace)  ::             model
         type(DataSpace)   ::                   dat
         type(ModelSpace_elevation)   ::            elev
         real              ::                  u(bounds%nmin1-4:bounds%nmax1+4)
         type(GeneralParam):: genpar 
       end subroutine InjSrc

       subroutine TimeDer(genpar,bounds,u)
         use GeneralParam_types
         use FD_types
         type(GeneralParam) :: genpar
         type(FDbounds)     ::        bounds
         real               ::               u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
       end subroutine TimeDer
       
       subroutine TimeSwap(genpar,bounds,u)
         use GeneralParam_types
         use FD_types
         type(GeneralParam) :: genpar
         type(FDbounds)     ::        bounds
         real               ::               u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3)
       end subroutine TimeSwap
       
      end interface

    type(GeneralParam)         :: genpar
    type(ModelSpace)           :: model
    type(DataSpace)            :: dat
    type(FDbounds)             :: bounds
    type(ModelSpace_elevation) :: elev
    type(HigdonParam)          :: hig
    type(ScaledFDcoefs)        :: scaled_fdcoefs
    type(UnscaledFDcoefs)      :: fdcoefs
    integer                    :: it,counter
    real, allocatable :: u(:,:,:,:)

    allocate(u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3))
    
    allocate(elev%irec_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
    allocate(elev%drec_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
    if (genpar%surf_type.ne.0) then
       allocate(elev%ielev_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
       allocate(elev%delev_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
    end if

    allocate(hig%gx(16, bounds%nmin1:bounds%nmax1, bounds%nmin3:bounds%nmax3))
    allocate(hig%gy(16, bounds%nmin1:bounds%nmax1, bounds%nmin2:bounds%nmax2))
    allocate(hig%gz(16, bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))

    u=0.

    call FD_coefs(fdcoefs)
    call FD_types_assign_scaled_coefs(fdcoefs,scaled_fdcoefs,genpar)
    call FD_derivatives_coef_init(scaled_fdcoefs)

    call Higdon(genpar%dt,model,bounds,hig)
    write(0,*) 'INFO: done with Higdon'
    call  ModelSpace_compute_shotz_positions(elev)
    write(0,*) 'INFO: done with source elevation'
    write(0,*) elev%ishot_z,elev%dshot_z
    call ModelSpace_elevation_parameters(elev,bounds,genpar)
    write(0,*) 'INFO: done with receiver elevation'

    write(0,*) elev%ishot_x,elev%ishot_y,genpar%dt
    counter=0
    TIME_LOOPS:do it=1,dat%nt     
       dat%it=it
       if (mod(it,50).eq.0) write (0,*) "Step",it," of ",dat%nt,"time steps (Source Wavefield)"
!       write(0,*) minval(u),maxval(u)
       call Extraction_shot_simple(bounds,dat,elev,u,genpar,it)
       call Extraction_wavefield(bounds,model,dat,elev,u,genpar,it,counter)
       call Boundary_set_free_surface(bounds,model,elev,u,genpar)
       call FD_scheme(genpar,bounds,u,model)
       call InjSrc(bounds,model,dat,elev,u(:,elev%ishot_x,elev%ishot_y,3),genpar)
       call TimeDer(genpar,bounds,u)
!       if (genpar%shot_type.eq.0) then
!          call Boundary0(genpar,bounds,u,model,hig)
!       else
!          call Boundary1(genpar,bounds,u,model,hig)
!       end if

       call TimeSwap(genpar,bounds,u)
    end do TIME_LOOPS
   
    deallocate(u)
    call deallocateHigdonParam(hig)
    call deallocateModelSpace_elev(elev)

  end subroutine propagator_acoustic
end module Propagator_mod
