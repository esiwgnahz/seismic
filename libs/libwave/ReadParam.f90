module ReadParam_mod

  use sep
  use GeneralParam_types

  implicit none

contains
  
  subroutine read_3D_params(genpar)
    type(GeneralParam)  ::  genpar

    call from_param('fmax',genpar%fmax,30.)
    call from_param('ntaper',genpar%ntaper,20)
    call from_param('snapi',genpar%snapi,4)
    call from_param('aperture_x',genpar%aperture(1))
    call from_param('aperture_y',genpar%aperture(2))
    call from_param('num_threads',genpar%nthreads,4)
    call from_param('lsinc',genpar%lsinc,7)
    call from_param('maxsize',genpar%max_memory,16.)

    call from_param('rec_type',genpar%rec_type,0)
    call from_param('shot_type',genpar%shot_type,0)
    call from_param('surf_type',genpar%surf_type,0)

    call from_param('twoD',genpar%twoD,.false.)
    call from_param('withRho',genpar%withRho,.false.)
    call from_param('LSRTM',genpar%LSRTM,.false.)
    call from_param('Born',genpar%Born,.false.)
    call from_param('verb',genpar%verbose,.true.)
    call from_param('Optim',genpar%optim,.true.)
    call from_param('Task',genpar%task,'MOD')

    if (genpar%LSRTM) call from_param('CHEAPLSRTM',genpar%CHEAPLSRTM,.false.)
       
    call from_param('write_forward_wavefield',genpar%WriteFwdWvfld,.false.)

    if (.not.genpar%twoD) then     
       genpar%nbound=4
    else
       genpar%nbound=0
    end if

    if (genpar%withRho) then
       genpar%coefpower=1
    else
       genpar%coefpower=2
    end if

    write(0,*) 'INFO: -------- Parameters ---------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO: task  =',genpar%task
    write(0,*) 'INFO: fmax  =',genpar%fmax,'Hz'
    write(0,*) 'INFO: ntaper=',genpar%ntaper,'points'
    write(0,*) 'INFO: snapi =',genpar%snapi,'time slices'
    write(0,*) 'INFO: aperture_x=',genpar%aperture(1)
    write(0,*) 'INFO: aperture_y=',genpar%aperture(2)
    write(0,*) 'INFO:'
    write(0,*) 'INFO: num_threads=',genpar%nthreads
    write(0,*) 'INFO: maxsize    =',genpar%max_memory,'Gb'
    write(0,*) 'INFO: lsinc      =',genpar%lsinc,'points'
    write(0,*) 'INFO:'
    write(0,*) 'INFO: twoD                   =',genpar%twoD
    write(0,*) 'INFO: withRho                =',genpar%withRho
    write(0,*) 'INFO: write_forward_wavefield=',genpar%WriteFwdWvfld
    write(0,*) 'INFO: LSRTM                  =',genpar%LSRTM
    write(0,*) 'INFO: Born                   =',genpar%Born
    if (genpar%LSRTM) write(0,*) 'INFO: CHEAPLSRTM             =',genpar%CHEAPLSRTM
    write(0,*) 'INFO: ----------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO: -------- Source/Receiver parameters ------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO: rec_type =',genpar%rec_type
    write(0,*) 'INFO: shot_type=',genpar%shot_type
    write(0,*) 'INFO: surf_type=',genpar%surf_type
    write(0,*) 'INFO:'
    write(0,*) 'INFO: ------------------------------------------'

  end subroutine read_3D_params

end module ReadParam_mod
