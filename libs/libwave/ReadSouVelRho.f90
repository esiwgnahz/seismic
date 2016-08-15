module Readsouvelrho_mod
  
  use sep

  use DataSpace_types
  use GeneralParam_types

  implicit none

contains

  subroutine readsou(source,genpar)
    type(TraceSpace)::source
    type(GeneralParam) ::   genpar
    real            ::sinc(genpar%lsinc)
    real            ::fdum,fscale
    real,dimension(:)  ,allocatable :: trace
    integer         :: it,l

    call from_history('n1',sourcevec%dimt%nt)
    call from_history('d1',sourcevec%dimt%dt)
    allocate(sourcevec%trace(sourcevec%dimt%nt,1))
    call sreed('in',sourcevec%trace(:,1),4*sourcevec%dimt%nt)
    
    genpar%dt=sourcevec%dimt%dt
    genpar%nt=sourcevec%dimt%nt

    call from_param('withRho',genpar%withRho,.false.)

    if (genpar%withRho) then
       allocate(trace(sourcevec%dimt%nt))
       ! Shift by 1/2 time step forward(reason: staggered grid)
       call mksinc(sinc,lsinc,+0.5)
       trace = sou(1:sourcevec%dimt%nt,1)
       do it=1,sourcevec%dimt%nt
          fscale = 1.
          fdum = 0.
          do l=-lsinc/2,lsinc/2
             if (it+l.ge.1 .and. it+l.le.sourcevec%dimt%nt) then
                fdum = fdum + sinc(lsinc/2+1+l)*trace(it+l)
             else
                fscale = fscale - sinc(lsinc/2+1+l)
             endif
          end do
          sourcevec%trace(it,1)= fdum / fscale
       end do
       deallocate(trace)
    end if

  end subroutine readsou

  subroutine readsoucoord(source,mod)    
    type(TraceSpace) ::   source
    type(ModelSpace) ::          mod
    real             :: mino(3),maxo(3),mid(3)
    integer          :: i

    mino(1)=mod%oz
    mino(2)=mod%ox
    mino(3)=mod%oy

    maxo(1)=(mod%nz-1)*mod%dz+mod%oz
    maxo(2)=(mod%nx-1)*mod%dx+mod%ox
    maxo(3)=(mod%ny-1)*mod%dy+mod%oy

    mid=(maxo-mino)/2

    call from_param('source_z',source%coord(1),mido(1))
    call from_param('source_x',source%coord(2),mido(2))
    call from_param('source_y',source%coord(3),mido(3))

    do i=1,3
       if ((source%coord(i).lt.mino(i)).or.(source%coord(i).gt.maxo(i))) then
          write(0,*) 'ERROR:  Problem with axis',i
          call erexit('ERROR: source_coord is not within model bounds, exit now')
       end if
    end do

    write(0,*) 'source coordinates =',source%coord

  end subroutine readsoucoord

  subroutine readvel(mod,genpar)
    type(ModelSpace)   ::mod
    type(GeneralParam) ::genpar
    
    if (auxpar('n1','i',mod%nz,mod%veltag).eq.0)  & 
    &    call erexit('need n1:nz')
    call putch('From aux(vel): nz','i',mod%nz)
    if (auxpar('o1','r',genpar%omodel(1),mod%veltag).eq.0.)  & 
    &    call erexit('need o1:omodel(1)')
    call putch('From aux(vel): oz','r',genpar%omodel(1))
    if (auxpar('d1','r',genpar%delta(1),mod%veltag).eq.0.)  & 
    &    call erexit('need d1:delta(1)')
    call putch('From aux(vel): dz','r',genpar%delta(1))
    
    if (auxpar('n2','i',mod%nx,mod%veltag).eq.0)  & 
    &    call erexit('need n2:nx')
    call putch('From aux(vel): nx','i',mod%nx)
    if (auxpar('o2','r',genpar%omodel(2),mod%veltag).eq.0.)  & 
    &    call erexit('need o2:omodel(2)')
    call putch('From aux(vel): oz','r',genpar%omodel(2))
    if (auxpar('d2','r',genpar%delta(2),mod%veltag).eq.0.)  & 
    &    call erexit('need d2:delta(2)')
    call putch('From aux(vel): dx','r',genpar%delta(2))

    if (genpar%twoD) then
       mod%ny=1
       genpar%delta(3)=1
       genpar%omodel(3)=0.
    else
       if (auxpar('n3','i',mod%ny,mod%veltag).eq.0)  & 
       &    call erexit('need n3:ny')
       call putch('From aux(vel): ny','i',mod%ny)
       if (auxpar('o3','r',genpar%omodel(3),mod%veltag).eq.0.)  & 
       &    call erexit('need o3:omodel(3)')
       call putch('From aux(vel): oz','r',genpar%omodel(3))
       if (auxpar('d3','r',genpar%delta(3),mod%veltag).eq.0.)  & 
       &    call erexit('need d3:delta(3)')
       call putch('From aux(vel): dz','r',genpar%delta(3))
    end if
       
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO:  Velocity Model Dimensions  '
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO: nz  =',mod%nz
    write(0,*) 'INFO: nx  =',mod%nx
    write(0,*) 'INFO: ny  =',mod%ny
    write(0,*) 'INFO:'
    write(0,*) 'INFO: dz  =',genpar%delta(1)
    write(0,*) 'INFO: dx  =',genpar%delta(2)
    write(0,*) 'INFO: dy  =',genpar%delta(3)
    write(0,*) 'INFO:'
    write(0,*) 'INFO: oz  =',genpar%omodel(1)
    write(0,*) 'INFO: ox  =',genpar%omodel(2)
    write(0,*) 'INFO: oy  =',genpar%omodel(3)
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO:'

    mod%oz=genpar%omodel(1)
    mod%ox=genpar%omodel(2)
    mod%oy=genpar%omodel(3)

    mod%dz=genpar%delta(1)
    mod%dx=genpar%delta(2)
    mod%dy=genpar%delta(3)

    allocate(mod%vel(mod%nz,mod%nx,mod%ny)
    if (.not.exist_file(mod%veltag)) then
       call erexit('ERROR: Need velocity file, exit now')
    else
       call sreed(mod%veltag,vel,4*mod%nz*mod%nx*mod%ny)
       call auxclose(mod%veltag)
       call vel_check(mod,genpar)
    end if
    
    if(genpar%withRho) then
       if (.not.exist_file(mod%rhotag)) call erexit('ERROR: Need rho file, exit now')
    else
       call dim_consistency_chech(mod%rhotag,mod%veltag,genpar%twoD)
       call sreed(mod%rhotag,rho,4*mod%nz*mod%nx*mod%ny)
       call auxclose(mod%rhotag)
    end if
  end subroutine readvel
   !
  ! Subroutine to check the parameters supplied for stability and accuracy
  ! If the computations are too inaccurate or unstable, then the program
  ! will exit
  !
  !
  subroutine vel_check(mod,genpar)
    type(ModelSpace)   :: mod
    type(GeneralParam) ::genpar

    implicit none
    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    !
    ! Variable dictionary
    !
    !  fmax		real		maximum frequency of the source wavelet
    !  fmax_new	real		Suggested fmax if the wavefield is aliased
    !				given the parameters supplied  
    !  stab		real		Courant stability factor
    !  dtn		real		Suggested dt if the computations are unstable
    !  samplerate	real		number of samples per smallest wavelength
    !  vmin		real		minumum velocity (but bigger than 0) used
    !				for determining smallest wavelength with
    !				fmax
    !  lambdamin	real		computed minimum wavelength
    !
    integer :: i, j
    real :: vmax, vmin, stab, samplerate,dtn
    real :: lambdamin, fmax_new
    !
    ! check the velocity model for its maximum value
    !
    vmax=maxval(mod%vel)
    vmin=minval(mod%vel)

    if (genpar%twoD) then
       stab=vmax*genpar%dt/(min(genpar%delta(1),genpar%delta(2)))
    else
       stab=vmax*genpar%dt/(min(genpar%delta(1),genpar%delta(2),genpar%delta(3)))
    end if

    call putlin('     ')
    call putlin('The Courant number vmax*dt/min(dx,dz)') 
    call putlin('which controls')
    call putlin('Stability and accuracy') 
    call putlin('of the modeled wavefield is:')
    call putch('stab = Courant_number','r',stab)
    !
    ! if the computations are unstable, exit
    !
    if (stab.gt.0.45) then
       call putlin(							&
       &      'The Courant number is >0.45 so calculations for this')
       call putlin(							&
       &      'set of vmax, dt, and dx is unstable')
       call putlin(							&
       &      'You should be able to make the computations stable')
       call putlin(							&
       &      'by decreasing the time step size of the source wave')
       call putlin(							&
       &      'For the parameters given the time step size should be')
       call putlin(							&
       &      'set to no greater than')
       !
       ! compute the largest stable time step size and write it out to the header
       !
       dtn=0.45/stab*genpar%dt
       !
       call putch('dtn = new_dt','r',dtn)
       call putlin('     ')
       call erexit('ERROR: computations will be unstable')
       !
    else
       call putlin(							&
       &      'The Courant number is <= 0.45, calculations will procede')
       call putlin(' ')
    end if
    !
    ! now find the minumum wavelength and see if it is adequately sampled
    !
    lambdamin=vmin/genpar%fmax
    !
    ! compute its sample rate for the model given
    !
    if (genpar%twoD) then
       samplerate=lambdamin/(max(genpar%delta(1),genpar%delta(2)))
    else
       samplerate=lambdamin/(max(genpar%delta(1),genpar%delta(2),genpar%delta(3)))
    end if
    !
    ! if we are adequately sampled put this on the header file for
    ! the output seismograms
    !
    if (samplerate.ge.3.0) then
       call putch('samplerate = min_sample_rate','r',samplerate)
       call putlin('The wavefield is adequately sampled')
       call putlin(' ')
    end if
    !
    ! if we are pushing the accuracy put that out on the seismogram header
    !
    if (2.0.le.samplerate.and.samplerate.lt.3.0) then
       call putch('samplerate = min_sample_rate','r',samplerate)
       call putlin(							&
       &      'The wavefield might not be sampled well enough at high')
       call putlin(							&
       &      'frequencies or for long propagation distances.')
       call putlin(							&
       &      'The computations will procede, hope for the best')
       call putlin(							&
       &      'You may see some dispersion in the result.')
       call putlin(' ')
    end if
    !
    ! stop if the computations are aliased or nearly so
    !
    if (samplerate.lt.2.0) then
       call putch('samplerate = min_sample_rate','r',samplerate)
       call putlin(							&
       &      'The wavefield will be severely undersampled or aliased')
       call putlin(							&
       &      'wherever there are low velocities.  You should lower the')
       call putlin(							&
       &      'frequency content of your source, or use lower values of')
       call putlin(							&
       &      'dx, dz, and dt if you wish to keep the same frequencies')
       call putlin(							&
       &      'in your source wave.')
       call putlin(							&
       &      'The wavefield will be adequately sampled if you lower')
       call putlin(							&
       &      'the maximum frequency to:')
       !
       ! compute the suggested maximum frequency
       !
       fmax_new=samplerate/3.0*fmax
       !
       call putch('fmax_new = fmax_new','r',fmax_new)
       call putlin(' ')
       call putlin('And keep dx, dz, and dt the same. ')
       call erexit('ERROR: wavefield not adequately sampled')
    end if
    !
    ! if we made if this far, we can procede with the computations
    !
  end subroutine vel_check

  subroutine dim_consistency_chech(tag1,tag2,twoD)
    character(len=3) :: tag1,tag2

    integer :: n1,n2
    real    :: o1,o2
    real    :: d1,d2
    logical :: stat,twoD

    stat=auxpar('n1','i',n1,tag1)
    stat=auxpar('n1','i',n2,tag2)
    if (n1.ne.n2) then
       write(0,*) 'ERROR:',tag1,tag2,'not consistent with n1'
       call erexit('ERROR: leaving now')
    end if
    stat=auxpar('o1','r',o1,tag1)
    stat=auxpar('o1','r',o2,tag2)
    if (o1.ne.o2) then
       write(0,*) 'ERROR:',tag1,tag2,'not consistent with o1'
       call erexit('ERROR: leaving now')
    end if
    stat=auxpar('d1','r',d1,tag1)
    stat=auxpar('d1','r',d2,tag2)
    if (d1.ne.d2) then
       write(0,*) 'ERROR:',tag1,tag2,'not consistent with d1'
       call erexit('ERROR: leaving now')
    end if
    stat=auxpar('n2','i',n1,tag1)
    stat=auxpar('n2','i',n2,tag2)
    if (n1.ne.n2) then
       write(0,*) 'ERROR:',tag1,tag2,'not consistent with n2'
       call erexit('ERROR: leaving now')
    end if
    stat=auxpar('o2','r',o1,tag1)
    stat=auxpar('o2','r',o2,tag2)
    if (n1.ne.n2) then
       write(0,*) 'ERROR:',tag1,tag2,'not consistent with o2'
       call erexit('ERROR: leaving now')
    end if
    stat=auxpar('d2','r',d1,tag1)
    stat=auxpar('d2','r',d2,tag2)
    if (d1.ne.d2) then
       write(0,*) 'ERROR:',tag1,tag2,'not consistent with d2'
       call erexit('ERROR: leaving now')
    end if

    if (.not.twoD) then
       stat=auxpar('n3','i',n1,tag1)
       stat=auxpar('n3','i',n2,tag2)
       if (n1.ne.n2) then
          write(0,*) 'ERROR:',tag1,tag2,'not consistent with n3'
          call erexit('ERROR: leaving now')
       end if
       stat=auxpar('o3','r',o1,tag1)
       stat=auxpar('o3','r',o2,tag2)
       if (n1.ne.n2) then
          write(0,*) 'ERROR:',tag1,tag2,'not consistent with o3'
          call erexit('ERROR: leaving now')
       end if
       stat=auxpar('d3','r',d1,tag1)
       stat=auxpar('d3','r',d2,tag2)
       if (d1.ne.d2) then
          write(0,*) 'ERROR:',tag1,tag2,'not consistent with d3'
          call erexit('ERROR: leaving now')
       end if
    end if

  end subroutine dim_consistency_chech

end module readsouvelrho_mod
