! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Propagator_mod
  
  ! Types
  use omp_lib
  use FD_types
  use ModelSpace_types
  use DataSpace_types
  use GeneralParam_types

  ! Computing
  use Imaging_mod
  use Boundary_mod
  use Injection_mod
  use FDcoefs_assign
  use Extraction_mod
  use FD_derivatives
  use FDswaptime_mod

  implicit none

contains

  subroutine propagator_acoustic(FD_coefs,FD_scheme,Injection,TimeDer,TimeSwap,bounds,model,elev,genpar,&
  &                      sou,wfld,datavec,ExtractData,ExtractWave,ImagingCondition)
    optional             sou,wfld,datavec,ExtractData,ExtractWave,ImagingCondition
    interface
       subroutine FD_coefs      (coef)
         use FD_types
         type(UnscaledFDcoefs)   ::      coef        
       end subroutine FD_coefs

       subroutine ImagingCondition(bounds,model,elev,u,genpar,it)
         use GeneralParam_types
         use ModelSpace_types
         use DataSpace_types

         type(FDbounds)    ::      bounds
         type(ModelSpace)  ::            model
         type(ModelSpace_elevation) ::          elev
         type(GeneralParam)::                         genpar
         real              ::                        u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
         &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
         integer           :: it
       end subroutine ImagingCondition

       subroutine FD_scheme     (genpar,bounds,u2,u3,model)
         use ModelSpace_types
         use GeneralParam_types
         use FD_types
         type(GeneralParam)   :: genpar
         type(ModelSpace)     ::                 model
         type(FDbounds)       ::        bounds
         real                 :: u2(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
         real                 :: u3(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
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

       subroutine ExtractWave(bounds,model,elev,u,genpar,it,dat)
         use GeneralParam_types
         use ModelSpace_types
         use DataSpace_types

         optional :: dat

         type(FDbounds)    :: bounds
         type(ModelSpace)  ::        model
         type(WaveSpace)   ::              dat
         type(ModelSpace_elevation) ::        elev
         type(GeneralParam)::                         genpar
         real              ::                       u(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, &
         &                                            bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
         integer           :: it
       end subroutine ExtractWave

       subroutine Injection(bounds,model,sou,u,genpar,it)
         use GeneralParam_types
         use ModelSpace_types
         use DataSpace_types
         use FD_types
         type(FDbounds)    ::      bounds
         type(ModelSpace)  ::             model
         type(TraceSpace), dimension(:)   ::                   sou
         type(GeneralParam):: genpar 
         real                 :: u(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
         integer           :: it
       end subroutine Injection

       subroutine TimeDer(genpar,bounds,grid)
         use GeneralParam_types
         use DataSpace_types
         use FD_types
         type(GeneralParam) :: genpar
         type(FDbounds)     ::        bounds
         type(USpace)       :: grid
       end subroutine TimeDer

       subroutine TimeSwap(grid)
         use DataSpace_types
         type(USpace)   :: grid
       end subroutine TimeSwap

    end interface

    type(GeneralParam)         :: genpar
    type(ModelSpace)           :: model
    type(TraceSpace), dimension(:) :: datavec
    type(TraceSpace), dimension(:) :: sou
    type(WaveSpace)            :: wfld
    type(FDbounds)             :: bounds
    type(ModelSpace_elevation) :: elev
    type(HigdonParam)          :: hig
    type(ScaledFDcoefs)        :: scaled_fdcoefs
    type(UnscaledFDcoefs)      :: fdcoefs
    integer                    :: it
    type(USpace) :: grid

    integer :: i,counting(10),count_rate,count_max,tmin,tmax,tstep
    real    :: totcount(9)
    logical :: verb,optim,testnan

    verb=genpar%verbose
    optim=genpar%optim

    testnan=.false.

    ! Allocate wavefields and set them to zero first
    ! It is important to set them to zero as small rounding errors creep in otherwise
    call allocateUSpace(grid,genpar,bounds)
    call zeroUSpace(grid)

    if (genpar%surf_type.ne.0) then
       if (.not.allocated(elev%ielev_z)) allocate(elev%ielev_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
       if (.not.allocated(elev%delev_z)) allocate(elev%delev_z(bounds%nmin2:bounds%nmax2,bounds%nmin3:bounds%nmax3))
    end if

    allocate(hig%gx(16, bounds%nmin1:bounds%nmax1, bounds%nmin3:bounds%nmax3))
    allocate(hig%gy(16, bounds%nmin1:bounds%nmax1, bounds%nmin2:bounds%nmax2))
    allocate(hig%gz(16, bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))

    call FD_coefs(fdcoefs)
    call FD_types_assign_scaled_coefs(fdcoefs,scaled_fdcoefs,genpar)
    call FD_derivatives_coef_init(scaled_fdcoefs)

    call Higdon(genpar%dt,model,bounds,hig)
    if (present(sou))     call ModelSpace_compute_array_xyz_positions(genpar,sou)
    if (present(datavec)) call ModelSpace_compute_array_xyz_positions(genpar,datavec)
    call ModelSpace_elevation_parameters(elev,bounds,genpar)
    !
    totcount=0.
    if (present(wfld)) wfld%counter=0

    TIME_LOOPS:do it=genpar%tmin,genpar%tmax,genpar%tstep    

       counting=0.

       if (verb.and.(mod(it,100).eq.0)) write (0,*) "INFO: Step",it," of ",max(genpar%tmin,genpar%tmax),"time steps"
 
!       testnan=there_is_nan_in_array2(grid%u2,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here1 testnan, it=',it
!          exit
!       end if
       
       !       write(0,*) 'here1',it
       if (verb) call system_clock(counting(1),count_rate,count_max)
       if (present(ExtractData)) then
          if (present(datavec)) call ExtractData(bounds,model,datavec,grid%u2,genpar,it) 
       end if

!       testnan=there_is_nan_in_array2(grid%u2,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here2 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here2',it
       if (verb) call system_clock(counting(2),count_rate,count_max)    
       if (present(ExtractWave)) then
          if (present(wfld)) then
             call ExtractWave(bounds,model,elev,grid%u2,genpar,it,dat=wfld)
          else
             call ExtractWave(bounds,model,elev,grid%u2,genpar,it)
          end if
       end if

!       testnan=there_is_nan_in_array2(grid%u2,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here3 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here3',it
       if (verb) call system_clock(counting(3),count_rate,count_max)
       if (present(ImagingCondition).and.(genpar%tstep.eq.-1)) then
          call ImagingCondition(bounds,model,elev,grid%u2,genpar,it)
       end if

!       testnan=there_is_nan_in_array2(grid%u2,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here4 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here4',it
       if (verb) call system_clock(counting(4),count_rate,count_max)      
       call Boundary_set_free_surface_grid(bounds,model,elev,grid,genpar)

!       testnan=there_is_nan_in_array2(grid%u2,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here5 u2 testnan, it=',it
!          exit
!       end if
!       testnan=there_is_nan_in_array2(grid%u3,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here5 u3 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here5',it
       if (verb) call system_clock(counting(5),count_rate,count_max)
       call FD_scheme(genpar,bounds,grid%u2,grid%u3,model)

       !       write(0,*) 'here6'
       if (verb) call system_clock(counting(6),count_rate,count_max)

!       testnan=there_is_nan_in_array2(grid%u3,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here6 u3 testnan, it=',it
!          exit
!       end if
       ! This is regular injection for FD modeling/receiver injection
       if (present(sou)) &
       & call Injection(bounds,model,sou,grid%u3(:,:,:),genpar,it)

!       testnan=there_is_nan_in_array2(grid%u3,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here7 u3 testnan, it=',it
!          exit
!       end if
       ! This is Born modeling, with genpar%tstep=+1
       if (present(ImagingCondition).and.(genpar%tstep.eq.1)) then
          call ImagingCondition(bounds,model,elev,grid%u3(:,:,:),genpar,it)
       end if

!       testnan=there_is_nan_in_array2(grid%u3,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here8 u3 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here7'
       if (verb) call system_clock(counting(7),count_rate,count_max)
       call TimeDer(genpar,bounds,grid)

!       testnan=there_is_nan_in_array2(grid%u3,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here9 u3 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here8'
       if (verb) call system_clock(counting(8),count_rate,count_max)
       if (optim) then
          call Boundary_3d(genpar,bounds,grid,model,hig)
          !call Boundary0_opt_grid(genpar,bounds,grid,model,hig)
       else
          call Boundary_3d_noomp(genpar,bounds,grid,model,hig)
          !call Boundary0_opt_grid_noomp(genpar,bounds,grid,model,hig)
       end if
    
!       testnan=there_is_nan_in_array2(grid%u3,bounds,genpar)
!       if (testnan) then
!          write(0,*) 'found here10 u3 testnan, it=',it
!          exit
!       end if
       !       write(0,*) 'here9'
       if (verb) call system_clock(counting(9),count_rate,count_max)
       call TimeSwap(grid)

       !       write(0,*) 'here10'
       if (verb) call system_clock(counting(10),count_rate,count_max)

       if (verb) then
          do i=1,9
             totcount(i)=totcount(i)+float(counting(i+1)-counting(i))/float(count_rate)/60.
          end do
       end if

    end do TIME_LOOPS

    if (verb) then
       write(0,*) 'INFO ---------------------------'
       write(0,*) 'INFO Total time               = ',sum(totcount),'mn'
       write(0,*) 'INFO ---------------------------'
       write(0,*) 'INFO  * Extract Data          = ',100*totcount(1)/sum(totcount),'%',totcount(1),'mn'
       write(0,*) 'INFO  * Extract wave          = ',100*totcount(2)/sum(totcount),'%',totcount(2),'mn'
       write(0,*) 'INFO  * Imaging condition     = ',100*totcount(3)/sum(totcount),'%',totcount(3),'mn'
       write(0,*) 'INFO  * Boundary free surface = ',100*totcount(4)/sum(totcount),'%',totcount(4),'mn'
       write(0,*) 'INFO  * FD stencil            = ',100*totcount(5)/sum(totcount),'%',totcount(5),'mn'
       write(0,*) 'INFO  * Inject Src            = ',100*totcount(6)/sum(totcount),'%',totcount(6),'mn'
       write(0,*) 'INFO  * Time derivative       = ',100*totcount(7)/sum(totcount),'%',totcount(7),'mn'
       write(0,*) 'INFO  * Absorbing boundaries  = ',100*totcount(8)/sum(totcount),'%',totcount(8),'mn'
       write(0,*) 'INFO  * Time Swap             = ',100*totcount(9)/sum(totcount),'%',totcount(9),'mn'
       write(0,*) 'INFO ---------------------------'
    end if

    call deallocateHigdonParam(hig)
    call deallocateUSpace(grid)

    if (genpar%surf_type.ne.0) then
       if (allocated(elev%ielev_z)) deallocate(elev%ielev_z)
       if (allocated(elev%delev_z)) deallocate(elev%delev_z)
    end if


  end subroutine propagator_acoustic

 
  function there_is_nan_in_array2(array,bounds,genpar) result(log)
    type(FDbounds)       ::             bounds
    type(GeneralParam)   ::                     genpar
    logical :: log
    integer :: n,j,k,i
    real    :: array(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)

    log=.false.

    !$OMP PARALLEL DO PRIVATE(k,j,i)
    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          do i=bounds%nmin1,bounds%nmax1
             if (isnan(array(i,j,k))) then
                log=.true.
                exit
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end function there_is_nan_in_array2

end module Propagator_mod
