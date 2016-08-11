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
  subroutine propagator_acoustic(FD_coefs,FD_scheme,ExtractData,InjSrc,TimeDer,TimeSwap,bounds,model,sou,datavec,wfld,elev,genpar)

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

       subroutine ExtractData(bounds,model,datavec,u,genpar,it)
         use GeneralParam_types
         use ModelSpace_types
         use DataSpace_types
         
         type(FDbounds)    ::             bounds
         type(ModelSpace)  ::                    model
         type(TraceSpace), dimension(:) ::             datavec
         type(GeneralParam)::                                 genpar 
         real              ::                               u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
         &                 bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
         integer           :: it
         
       end subroutine ExtractData

       subroutine InjSrc(bounds,model,sou,u,genpar,it)
         use GeneralParam_types
         use ModelSpace_types
         use DataSpace_types
         use FD_types
         type(FDbounds)    ::      bounds
         type(ModelSpace)  ::             model
         type(TraceSpace)   ::                   sou
         type(GeneralParam):: genpar 
         real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
         integer           :: it
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
    type(TraceSpace), dimension(:) :: datavec
    type(TraceSpace)           :: sou
    type(WaveSpace)            :: wfld
    type(FDbounds)             :: bounds
    type(ModelSpace_elevation) :: elev
    type(HigdonParam)          :: hig
    type(ScaledFDcoefs)        :: scaled_fdcoefs
    type(UnscaledFDcoefs)      :: fdcoefs
    integer                    :: it,counter
    real, allocatable :: u(:,:,:,:)

    allocate(u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound,-1:3))
    
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
    call ModelSpace_compute_xyz_positions(elev,sou)
    call ModelSpace_compute_array_xyz_position(elev,datavec)
    call ModelSpace_elevation_parameters(elev,bounds,genpar)

    counter=0
    TIME_LOOPS:do it=1,sou%dimt%nt     

       if (mod(it,50).eq.0) write (0,*) "Step",it," of ",sou%dimt%nt,"time steps (Source Wavefield)"

       call ExtractData(bounds,model,datavec,u(:,:,:,2),genpar,it)     
       call Extraction_wavefield(bounds,model,wfld,elev,u,genpar,counter,it)
       call Boundary_set_free_surface(bounds,model,elev,u,genpar)
       call FD_scheme(genpar,bounds,u,model)
       call InjSrc(bounds,model,sou,u(:,:,:,3),genpar,it)
       call TimeDer(genpar,bounds,u)
       if (genpar%shot_type.eq.0) then
          call Boundary0(genpar,bounds,u,model,hig)
       else
          call Boundary1(genpar,bounds,u,model,hig)
       end if

       call TimeSwap(genpar,bounds,u)

    end do TIME_LOOPS
   
    deallocate(u)
    call deallocateHigdonParam(hig)
    call deallocateModelSpace_elev(elev)

  end subroutine propagator_acoustic
end module Propagator_mod
