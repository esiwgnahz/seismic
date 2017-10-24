! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module ExtractPadModel_mod

  use sep
  use Readsouvelrho_mod

  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  implicit none
  
contains
  
  subroutine extract_coord_source_receiver_patch(datavec,sourcevec,mod,genpar)
    type(TraceSpace), dimension(:), allocatable::datavec
    type(ModelSpace)                           ::                  mod
    type(TraceSpace)             ::        sourcevec
    type(GeneralParam)                         ::                      genpar

    integer :: i
    real    :: sx,sy
    real    :: gx1,gxn,gy1,gyn
    real    :: minx,miny,maxx,maxy
    real    :: mx1,mxn,my1,myn
    real    :: centerx,centery
    
    sx=sourcevec%coord(2)
    sy=sourcevec%coord(3)
    gx1=datavec(1)%coord(2)
    gy1=datavec(1)%coord(3)
    gxn=datavec(1)%coord(2)
    gyn=datavec(1)%coord(3)
    ! First find the limits of the receiver box
    do i=2,size(datavec)
       gx1=min(gx1,datavec(i)%coord(2))
       gxn=max(gxn,datavec(i)%coord(2))
       gy1=min(gy1,datavec(i)%coord(3))
       gyn=max(gyn,datavec(i)%coord(3))
    end do

    minx = min(sx, gx1)
    miny = min(sy, gy1)
    maxx = max(sx, gxn)
    maxy = max(sy, gyn)

    mx1 = (gx1+sx)/2                                  ! midpoint location for first receiver in x
    mxn = (gxn+sx)/2                                  ! midpoint location for last receiver in x
    my1 = (gy1+sy)/2                                  ! midpoint location for first receiver in y
    myn = (gyn+sy)/2                                  ! midpoint location for last receiver in y
    centerx = mx1 + abs(mxn-mx1)/2                    ! Middle between midpoints in x 
    centery = my1 + abs(myn-my1)/2                    ! Middle between midpoints in y 

    mod%oxw=min(minx,centerx-genpar%aperture(1)/2)
    mod%endx=max(maxx,centerx+genpar%aperture(1)/2)
    mod%oyw=min(miny,centery-genpar%aperture(2)/2)
    mod%endy=max(maxy,centery+genpar%aperture(2)/2)
    
    if (genpar%verbose) then
       write(0,*) 'INFO:-------------------------------------------------------------------------'
       write(0,*) 'INFO: Coordinates model according to aperture and source/receiver positions is'
       write(0,*) 'INFO:-------------------------------------------------------------------------'
       write(0,*) 'INFO:'
       write(0,*) 'INFO: endy,oxw ------------------ endy,endx'
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    '
       write(0,*) 'INFO:     |                           |    ' 
       write(0,*) 'INFO:  oyw,oxw ------------------- oyw,endx'
       write(0,*) 'INFO:'
       write(0,*) 'INFO: oxw =',mod%oxw,'  oyw =',mod%oyw
       write(0,*) 'INFO: endx=',mod%endx,' endy=',mod%endy
       write(0,*) 'INFO:'
       write(0,*) 'INFO:'
    end if

  end subroutine extract_coord_source_receiver_patch

  subroutine read_window_vel(mod,genpar,bounds)
    type(ModelSpace)  ::     mod
    type(GeneralParam)::         genpar
    type(FDbounds)    ::                bounds
    real,allocatable,dimension(:,:,:) :: tmpbig,tmpsmall

    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    integer  :: nout,i
    
    if (auxpar('n1','i',mod%nz,mod%veltag).eq.0)  & 
    &    call erexit('need n1:nz')
    call putch('From aux(vel): nz','i',mod%nz)
    if (auxpar('o1','r',mod%oz,mod%veltag).eq.0.)  & 
    &    call erexit('need o1:omodel(1)')
    call putch('From aux(vel): oz','r',mod%oz)
    if (auxpar('d1','r',mod%dz,mod%veltag).eq.0.)  & 
    &    call erexit('need d1:delta(1)')
    call putch('From aux(vel): dz','r',mod%dz)
    
    if (auxpar('n2','i',mod%nx,mod%veltag).eq.0)  & 
    &    call erexit('need n2:nx')
    call putch('From aux(vel): nx','i',mod%nx)
    if (auxpar('o2','r',mod%ox,mod%veltag).eq.0.)  & 
    &    call erexit('need o2:omodel(2)')
    call putch('From aux(vel): oz','r',mod%ox)
    if (auxpar('d2','r',mod%dx,mod%veltag).eq.0.)  & 
    &    call erexit('need d2:delta(2)')
    call putch('From aux(vel): dx','r',mod%dx)

    if (genpar%twoD) then
       mod%ny=1
       mod%dy=1.
       mod%oy=0.
    else
       if (auxpar('n3','i',mod%ny,mod%veltag).eq.0)  & 
       &    call erexit('need n3:ny')
       call putch('From aux(vel): ny','i',mod%ny)
       if (auxpar('o3','r',mod%oy,mod%veltag).eq.0.)  & 
       &    call erexit('need o3:omodel(3)')
       call putch('From aux(vel): oz','r',mod%oy)
       if (auxpar('d3','r',mod%dy,mod%veltag).eq.0.)  & 
       &    call erexit('need d3:delta(3)')
       call putch('From aux(vel): dz','r',mod%dy)
    end if
       
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO:  Velocity Model Dimensions  '
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO: nz  =',mod%nz
    write(0,*) 'INFO: nx  =',mod%nx
    write(0,*) 'INFO: ny  =',mod%ny
    write(0,*) 'INFO:'
    write(0,*) 'INFO: dz  =',mod%dz
    write(0,*) 'INFO: dx  =',mod%dx
    write(0,*) 'INFO: dy  =',mod%dy
    write(0,*) 'INFO:'
    write(0,*) 'INFO: oz  =',mod%oz
    write(0,*) 'INFO: ox  =',mod%ox
    write(0,*) 'INFO: oy  =',mod%oy
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO:'

    nout=0
    if ((mod%oxw.lt.mod%ox).or.(mod%oxw.gt.(mod%ox+(mod%nx-1)*mod%dx))) nout=nout+1
    if ((mod%endx.lt.mod%ox).or.(mod%endx.gt.(mod%ox+(mod%nx-1)*mod%dx))) nout=nout+1
    if ((mod%oyw.lt.mod%oy).or.(mod%oyw.gt.(mod%oy+(mod%ny-1)*mod%dy))) nout=nout+1
    if ((mod%endy.lt.mod%oy).or.(mod%endy.gt.(mod%oy+(mod%ny-1)*mod%dy))) nout=nout+1
    if ((mod%oxw.lt.mod%ox).and.(mod%endx.gt.(mod%ox+(mod%nx-1)*mod%dx))) nout=nout-1
    if ((mod%oyw.lt.mod%oy).and.(mod%endy.gt.(mod%oy+(mod%ny-1)*mod%dy))) nout=nout-1

    if (nout.eq.4) then
       write(0,*) 'ERROR:'
       write(0,*) 'ERROR: The source/receiver patch, with aperture, is outside velocity model'
       write(0,*) 'ERROR:'
       write(0,*) 'ERROR: Source/receiver patch coordinates'
       write(0,*) 'ERROR: mod%oxw=',mod%oxw
       write(0,*) 'ERROR: mod%oyw=',mod%oyw
       write(0,*) 'ERROR: mod%endx=',mod%endx
       write(0,*) 'ERROR: mod%endy=',mod%endy
       write(0,*) 'ERROR:'
       write(0,*) 'ERROR: Model dimensions coordinates'
       write(0,*) 'ERROR: begw=',mod%ox
       write(0,*) 'ERROR: begw=',mod%oy
       write(0,*) 'ERROR: endx=',mod%ox+(mod%nx-1)*mod%dx
       write(0,*) 'ERROR: endy=',mod%oy+(mod%ny-1)*mod%dy
       call erexit('ERROR: leaving now due to improper dimensions of velocity vs. receiver-source patch')
    end if

    mod%nxw=nint((mod%endx-mod%oxw)/mod%dx)+1
    mod%nyw=nint((mod%endy-mod%oyw)/mod%dy)+1
    if (genpar%twoD) then
       mod%nyw=1
       mod%oyw=mod%oy
    end if

    write(0,*) 'INFO:'
    write(0,*) 'INFO:--------------------------------------'
    write(0,*) 'INFO:  Windowed velocity Model Dimensions  '
    write(0,*) 'INFO:--------------------------------------'
    write(0,*) 'INFO: nz  =',mod%nz
    write(0,*) 'INFO: nx  =',mod%nxw
    write(0,*) 'INFO: ny  =',mod%nyw
    write(0,*) 'INFO:'
    write(0,*) 'INFO: dz  =',mod%dz
    write(0,*) 'INFO: dx  =',mod%dx
    write(0,*) 'INFO: dy  =',mod%dy
    write(0,*) 'INFO:'
    write(0,*) 'INFO: oz  =',mod%oz
    write(0,*) 'INFO: ox  =',mod%oxw
    write(0,*) 'INFO: oy  =',mod%oyw
    write(0,*) 'INFO:-----------------------------------'
    write(0,*) 'INFO:'

    allocate(tmpbig(  mod%nz,mod%nx, mod%ny))
    allocate(tmpsmall(mod%nz,mod%nxw,mod%nyw))

    tmpbig=0.
    tmpsmall=0.

    genpar%delta(1)=mod%dz
    genpar%delta(2)=mod%dx
    genpar%delta(3)=mod%dy

    if (.not.exist_file(mod%veltag)) then
       call erexit('ERROR: Need velocity file, exit now')
    else
       do i=1,mod%ny
          call sreed(mod%veltag,tmpbig(:,:,i),4*mod%nz*mod%nx)
       end do
       call auxclose(mod%veltag)
       call vel_check(tmpbig,genpar)
    end if

    call mod_window_pad(.true.,tmpbig,tmpsmall,mod)

    genpar%omodel(1)=mod%oz
    genpar%omodel(2)=mod%oxw
    genpar%omodel(3)=mod%oyw

    if (genpar%shot_type.gt.0 .or. genpar%surf_type.gt.0) then
       if (genpar%shot_type.gt.0) bounds%nmin1 = -(genpar%ntaper+genpar%lsinc/2+2)+1
       if (genpar%surf_type.gt.0) bounds%nmin1 = -genpar%lsinc+1
    else
       bounds%nmin1 = -genpar%ntaper+1
    endif
    bounds%nmax1 =  mod%nz+genpar%ntaper
    bounds%nmin2 = -genpar%ntaper+1
    bounds%nmax2 =  mod%nxw+genpar%ntaper
    if (.not. genpar%twoD) then
       bounds%nmin3 = -genpar%ntaper+1
       bounds%nmax3 =  mod%nyw+genpar%ntaper
    else
       bounds%nmin3 = 1
       bounds%nmax3 = 1
    end if
  
    genpar%delta(1)=mod%dz
    genpar%delta(2)=mod%dx
    genpar%delta(3)=mod%dy

    allocate(mod%vel(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)) 
   
    call model_pad(tmpsmall,mod%vel,bounds,mod%nz,mod%nxw,mod%nyw,genpar)
   
    tmpsmall=0.
    tmpbig=0.
    !call srite('tmpvel',mod%vel,4*(bounds%nmax1-bounds%nmin1+9)*(bounds%nmax2-bounds%nmin2+9)*(bounds%nmax3-bounds%nmin3+2*genpar%nbound+1))
    !call to_history('n1',bounds%nmax1-bounds%nmin1+9,'tmpvel')
    !call to_history('n2',bounds%nmax2-bounds%nmin2+9,'tmpvel')
    !call to_history('n3',bounds%nmax3-bounds%nmin3+2*genpar%nbound+1,'tmpvel')
    if(genpar%withRho) then
       if (.not.exist_file(mod%rhotag)) then
          call erexit('ERROR: Need rho file, exit now')
       else
          call dim_consistency_check(mod%rhotag,mod%veltag,genpar%twoD)
          allocate(mod%rho(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
          allocate(mod%rho2(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
          call sreed(mod%rhotag,tmpbig,4*mod%nz*mod%nx*mod%ny)
          call auxclose(mod%rhotag)
          call mod_window_pad(.true.,tmpbig,tmpsmall,mod)
          call model_pad(tmpsmall,mod%rho,bounds,mod%nz,mod%nxw,mod%nyw,genpar)
          call Interpolate(mod,bounds,genpar)
       end if
    end if

    if(mod%exist_vel2) then
       call dim_consistency_check(mod%vel2tag,mod%veltag,genpar%twoD)
       allocate(mod%vel2(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
       call sreed(mod%vel2tag,tmpbig,4*mod%nz*mod%nx*mod%ny)
       call auxclose(mod%vel2tag)
       call vel_check(tmpbig,genpar)
       call mod_window_pad(.true.,tmpbig,tmpsmall,mod)
       call model_pad(tmpsmall,mod%vel2,bounds,mod%nz,mod%nxw,mod%nyw,genpar)
    end if

    deallocate(tmpbig,tmpsmall)
  end subroutine read_window_vel

  subroutine read_vel(   mod,genpar)
    type(ModelSpace)  :: mod
    type(GeneralParam)::     genpar

    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    
    integer :: i

    if (auxpar('n1','i',mod%nz,mod%veltag).eq.0)  & 
    &    call erexit('need n1:nz')
    call putch('From aux(vel): nz','i',mod%nz)
    if (auxpar('o1','r',mod%oz,mod%veltag).eq.0.)  & 
    &    call erexit('need o1:omodel(1)')
    call putch('From aux(vel): oz','r',mod%oz)
    if (auxpar('d1','r',mod%dz,mod%veltag).eq.0.)  & 
    &    call erexit('need d1:delta(1)')
    call putch('From aux(vel): dz','r',mod%dz)
    
    if (auxpar('n2','i',mod%nx,mod%veltag).eq.0)  & 
    &    call erexit('need n2:nx')
    call putch('From aux(vel): nx','i',mod%nx)
    if (auxpar('o2','r',mod%ox,mod%veltag).eq.0.)  & 
    &    call erexit('need o2:omodel(2)')
    call putch('From aux(vel): oz','r',mod%ox)
    if (auxpar('d2','r',mod%dx,mod%veltag).eq.0.)  & 
    &    call erexit('need d2:delta(2)')
    call putch('From aux(vel): dx','r',mod%dx)

    if (genpar%twoD) then
       mod%ny=1
       mod%dy=1.
       mod%oy=0.
    else
       if (auxpar('n3','i',mod%ny,mod%veltag).eq.0)  & 
       &    call erexit('need n3:ny')
       call putch('From aux(vel): ny','i',mod%ny)
       if (auxpar('o3','r',mod%oy,mod%veltag).eq.0.)  & 
       &    call erexit('need o3:omodel(3)')
       call putch('From aux(vel): oz','r',mod%oy)
       if (auxpar('d3','r',mod%dy,mod%veltag).eq.0.)  & 
       &    call erexit('need d3:delta(3)')
       call putch('From aux(vel): dz','r',mod%dy)
    end if
       
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO:  Velocity Model Dimensions  '
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO: nz  =',mod%nz
    write(0,*) 'INFO: nx  =',mod%nx
    write(0,*) 'INFO: ny  =',mod%ny
    write(0,*) 'INFO:'
    write(0,*) 'INFO: dz  =',mod%dz
    write(0,*) 'INFO: dx  =',mod%dx
    write(0,*) 'INFO: dy  =',mod%dy
    write(0,*) 'INFO:'
    write(0,*) 'INFO: oz  =',mod%oz
    write(0,*) 'INFO: ox  =',mod%ox
    write(0,*) 'INFO: oy  =',mod%oy
    write(0,*) 'INFO:-----------------------------'
    write(0,*) 'INFO:'

    allocate(mod%vel(  mod%nz,mod%nx, mod%ny))

    genpar%delta(1)=mod%dz
    genpar%delta(2)=mod%dx
    genpar%delta(3)=mod%dy

    if (.not.exist_file(mod%veltag)) then
       call erexit('ERROR: Need velocity file, exit now')
    else
       do i=1,mod%ny
          call sreed(mod%veltag,mod%vel(:,:,i),4*mod%nz*mod%nx)
       end do
       call auxclose(mod%veltag)
       call vel_check(mod%vel,genpar)
    end if

  end subroutine read_vel

  subroutine read_window_ref( mod,genpar,bounds)
    type(ModelSpace)  ::      mod
    type(GeneralParam)::          genpar
    type(FDbounds)    ::                 bounds
    real,allocatable,dimension(:,:,:) :: tmpbig,tmpsmall

    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    
    allocate(tmpbig(mod%nz,mod%nx,mod%ny))
    allocate(tmpsmall(mod%nz,mod%nxw,mod%nyw))

    tmpbig=0.
    tmpsmall=0.
    if (.not.exist_file(mod%reftag)) then
       call erexit('ERROR: Need reflectivity file, exit now')
    else
       call dim_consistency_check(mod%veltag,mod%reftag,genpar%twoD)
       allocate(mod%image(bounds%nmin1-4:bounds%nmax1+4,bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
       call sreed(mod%reftag,tmpbig,4*mod%nz*mod%nx*mod%ny)
       call auxclose(mod%reftag)
       call mod_window_pad(.true.,tmpbig,tmpsmall,mod)
       call model_pad(tmpsmall,mod%image,bounds,mod%nz,mod%nxw,mod%nyw,genpar)
    end if
    deallocate(tmpbig,tmpsmall)
  end subroutine read_window_ref

  subroutine copy_window_vel_gath(    mod,modgath,genpar,genpargath,boundsgath,ishot)
    type(ModelSpace)                ::mod
    type(ModelSpace)                ::    modgath
    type(GeneralParam)              ::            genpar
    type(GeneralParam)              ::                   genpargath
    type(FDbounds)                  ::                              boundsgath
    real,allocatable,dimension(:,:,:) :: tmpsmall

    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    integer  :: nout,nshots,ishot

    nout=0

    genpargath=genpar
    modgath%nx=mod%nx
    modgath%dx=mod%dx
    modgath%ox=mod%ox
    modgath%nz=mod%nz
    modgath%dz=mod%dz
    modgath%oz=mod%oz
    modgath%ny=mod%ny
    modgath%dy=mod%dy
    modgath%oy=mod%oy

    if ((modgath%oxw.lt.modgath%ox).or.(modgath%oxw.gt.(modgath%ox+(modgath%nx-1)*modgath%dx))) nout=nout+1
    if ((modgath%endx.lt.modgath%ox).or.(modgath%endx.gt.(modgath%ox+(modgath%nx-1)*modgath%dx))) nout=nout+1
    if ((modgath%oyw.lt.modgath%oy).or.(modgath%oyw.gt.(modgath%oy+(modgath%ny-1)*modgath%dy))) nout=nout+1
    if ((modgath%endy.lt.modgath%oy).or.(modgath%endy.gt.(modgath%oy+(modgath%ny-1)*modgath%dy))) nout=nout+1

    if ((modgath%oxw.lt.modgath%ox).and.(modgath%endx.gt.(modgath%ox+(modgath%nx-1)*modgath%dx))) nout=nout-1
    if ((modgath%oyw.lt.modgath%oy).and.(modgath%endy.gt.(modgath%oy+(modgath%ny-1)*modgath%dy))) nout=nout-1

    if (nout.eq.4) then
       write(0,*) 'ERROR:'
       write(0,*) 'ERROR: The source/receiver patch, with aperture, is outside velocity model'
       write(0,*) 'ERROR:'
       write(0,*) 'ERROR: Source/receiver patch coordinates for shot #',ishot
       write(0,*) 'ERROR: modgath%oxw=',modgath%oxw
       write(0,*) 'ERROR: modgath%oyw=',modgath%oyw
       write(0,*) 'ERROR: modgath%endx=',modgath%endx
       write(0,*) 'ERROR: modgath%endy=',modgath%endy
       write(0,*) 'ERROR:'
       write(0,*) 'ERROR: Model dimensions coordinates'
       write(0,*) 'ERROR: begw=',modgath%ox
       write(0,*) 'ERROR: begw=',modgath%oy
       write(0,*) 'ERROR: endx=',modgath%ox+(modgath%nx-1)*modgath%dx
       write(0,*) 'ERROR: endy=',modgath%oy+(modgath%ny-1)*modgath%dy
       call erexit('ERROR: leaving now due to improper dimensions of velocity vs. receiver-source patch')
    end if

    modgath%nxw=nint((modgath%endx-modgath%oxw)/modgath%dx)+1
    modgath%nyw=nint((modgath%endy-modgath%oyw)/modgath%dy)+1
    if (genpar%twoD) then
       modgath%nyw=1
       modgath%oyw=modgath%oy
    end if

    if (genpar%verbose) then
       write(0,*) 'INFO:'
       write(0,*) 'INFO:--------------------------------------------------'
       write(0,*) 'INFO:  Windowed velocity Model Dimensions for shot #',ishot
       write(0,*) 'INFO:--------------------------------------------------'
       write(0,*) 'INFO: nz  =',modgath%nz
       write(0,*) 'INFO: nx  =',modgath%nxw
       write(0,*) 'INFO: ny  =',modgath%nyw
       write(0,*) 'INFO:'
       write(0,*) 'INFO: dz  =',modgath%dz
       write(0,*) 'INFO: dx  =',modgath%dx
       write(0,*) 'INFO: dy  =',modgath%dy
       write(0,*) 'INFO:'
       write(0,*) 'INFO: oz  =',modgath%oz
       write(0,*) 'INFO: ox  =',modgath%oxw
       write(0,*) 'INFO: oy  =',modgath%oyw
       write(0,*) 'INFO:---------------------------------------------------'
       write(0,*) 'INFO:'
    end if

    allocate(tmpsmall(modgath%nz,modgath%nxw,modgath%nyw))

    tmpsmall=0.

    call mod_window_pad(.true.,mod%vel,tmpsmall,modgath)

    genpargath%omodel(1)=modgath%oz
    genpargath%omodel(2)=modgath%oxw
    genpargath%omodel(3)=modgath%oyw

    if (genpar%shot_type.gt.0 .or. genpar%surf_type.gt.0) then
       if (genpar%shot_type.gt.0) boundsgath%nmin1 = -(genpar%ntaper+genpar%lsinc/2+2)+1
       if (genpar%surf_type.gt.0) boundsgath%nmin1 = -genpar%lsinc+1
    else
       boundsgath%nmin1 = -genpar%ntaper+1
    endif
    boundsgath%nmax1 =  modgath%nz+genpar%ntaper
    boundsgath%nmin2 = -genpar%ntaper+1
    boundsgath%nmax2 =  modgath%nxw+genpar%ntaper
    if (.not. genpar%twoD) then
       boundsgath%nmin3 = -genpar%ntaper+1
       boundsgath%nmax3 =  modgath%nyw+genpar%ntaper
    else
       boundsgath%nmin3 = 1
       boundsgath%nmax3 = 1
    end if

    allocate(modgath%vel(boundsgath%nmin1-4:boundsgath%nmax1+4, boundsgath%nmin2-4:boundsgath%nmax2+4, boundsgath%nmin3-genpar%nbound:boundsgath%nmax3+genpar%nbound)) 

    call model_pad(tmpsmall,modgath%vel,boundsgath,mod%nz,modgath%nxw,modgath%nyw,genpar)

    deallocate(tmpsmall)


  end subroutine copy_window_vel_gath

  subroutine mod_window_pad(adj,big,small,mod)
    type(ModelSpace) ::                   mod
    logical         ::      adj
    real, dimension(:,:,:) ::   big,small
    
    integer :: i,j,ix,iy,ix1,iy1
    real    :: coord_x, coord_y

    if (adj) then
       !$OMP PARALLEL DO PRIVATE(j,coord_y,iy,i,coord_x,ix)
       do j=1,mod%nyw
          coord_y=(j-1)*mod%dy+mod%oyw
          iy=nint((coord_y-mod%oy)/mod%dy)+1
          iy1=min(mod%ny,max(1,iy))
          do i=1,mod%nxw
             coord_x=(i-1)*mod%dx+mod%oxw
             ix=nint((coord_x-mod%ox)/mod%dx)+1
             ix1=min(mod%nx,max(1,ix))
             
             small(:,i,j)=small(:,i,j)+big(:,ix1,iy1)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(j,coord_y,iy,i,coord_x,ix)
       do j=1,mod%nyw
          coord_y=(j-1)*mod%dy+mod%oyw
          iy=nint((coord_y-mod%oy)/mod%dy)+1
          if ((iy.lt.1).or.(iy.gt.mod%ny)) cycle
          do i=1,mod%nxw
             coord_x=(i-1)*mod%dx+mod%oxw
             ix=nint((coord_x-mod%ox)/mod%dx)+1
             if ((ix.lt.1).or.(ix.gt.mod%nx)) cycle

             big(:,ix,iy)=big(:,ix,iy)+small(:,i,j)
          end do
       end do       
       !$OMP END PARALLEL DO     
    end if

  end subroutine mod_window_pad

  subroutine mod_copy_image_to_disk(mod)
    type(ModelSpace) ::                   mod
    real, dimension(:,:,:), allocatable :: tmpim,tmpil

    integer :: i,j,ix,iy,ix1,iy1
    real    :: coord_x, coord_y

    allocate(tmpim(mod%nz,mod%nx,mod%ny))
    allocate(tmpil(mod%nz,mod%nx,mod%ny))
    tmpim=0.
    tmpil=0.

    !$OMP PARALLEL DO PRIVATE(j,coord_y,iy,i,coord_x,ix)
    do j=1,mod%nyw
       coord_y=(j-1)*mod%dy+mod%oyw
       iy=nint((coord_y-mod%oy)/mod%dy)+1
       if ((iy.lt.1).or.(iy.gt.mod%ny)) cycle

       do i=1,mod%nxw
          coord_x=(i-1)*mod%dx+mod%oxw
          ix=nint((coord_x-mod%ox)/mod%dx)+1
          if ((ix.lt.1).or.(ix.gt.mod%nx)) cycle
          
          tmpim(:,ix,iy)=tmpim(:,ix,iy)+mod%imagesmall(:,i,j)
          tmpil(:,ix,iy)=tmpil(:,ix,iy)+mod%illumsmall(:,i,j)
       end do
    end do
    !$OMP END PARALLEL DO
 
    do j=1,mod%ny
       call srite('image',tmpim(:,:,j),4*mod%nx*mod%nz)
    end do

    do j=1,mod%ny
       call srite('image',tmpil(:,:,j),4*mod%nx*mod%nz)
    end do

    deallocate(tmpim,tmpil) 

  end subroutine mod_copy_image_to_disk

  subroutine mod_copy_gradient_to_disk(mod,nparam)
    type(ModelSpace) ::                mod
    integer          :: nparam
    real, dimension(:,:,:,:), allocatable :: tmpim
    real, dimension(:,:,:),   allocatable :: tmpil

    integer :: i,j,k,ix,iy,ix1,iy1
    real    :: coord_x, coord_y

    allocate(tmpim(mod%nz,mod%nx,mod%ny,nparam))
    allocate(tmpil(mod%nz,mod%nx,mod%ny))
    tmpim=0.
    tmpil=0.

    if (nparam.eq.1) then
       !$OMP PARALLEL DO PRIVATE(j,coord_y,iy,i,coord_x,ix)
       do j=1,mod%nyw
          coord_y=(j-1)*mod%dy+mod%oyw
          iy=nint((coord_y-mod%oy)/mod%dy)+1
          if ((iy.lt.1).or.(iy.gt.mod%ny)) cycle
          
          do i=1,mod%nxw
             coord_x=(i-1)*mod%dx+mod%oxw
             ix=nint((coord_x-mod%ox)/mod%dx)+1
             if ((ix.lt.1).or.(ix.gt.mod%nx)) cycle
             
             tmpim(:,ix,iy,1)=tmpim(:,ix,iy,1)+mod%imagesmall(:,i,j)
             tmpil(:,ix,iy)=tmpil(:,ix,iy)+mod%illumsmall(:,i,j)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(j,coord_y,iy,i,coord_x,ix)
       do j=1,mod%nyw
          coord_y=(j-1)*mod%dy+mod%oyw
          iy=nint((coord_y-mod%oy)/mod%dy)+1
          if ((iy.lt.1).or.(iy.gt.mod%ny)) cycle
          
          do i=1,mod%nxw
             coord_x=(i-1)*mod%dx+mod%oxw
             ix=nint((coord_x-mod%ox)/mod%dx)+1
             if ((ix.lt.1).or.(ix.gt.mod%nx)) cycle
             
             tmpim(:,ix,iy,1)=tmpim(:,ix,iy,1)+mod%imagesmall_nparam(:,i,j,1)
             tmpim(:,ix,iy,2)=tmpim(:,ix,iy,2)+mod%imagesmall_nparam(:,i,j,2)
             tmpil(:,ix,iy)=tmpil(:,ix,iy)+mod%illumsmall(:,i,j)
          end do
       end do
       !$OMP END PARALLEL DO
    end if
       
    do i=1,nparam
       do j=1,mod%ny
          call srite('gradient',tmpim(:,:,j,i),4*mod%nx*mod%nz)
       end do
    end do

    do j=1,mod%ny
       call srite('gradient',tmpil(:,:,j),4*mod%nx*mod%nz)
    end do

    deallocate(tmpim,tmpil) 

  end subroutine mod_copy_gradient_to_disk

  subroutine mod_copy_image(mod,image,illum)
    type(ModelSpace) ::     mod
    real, dimension(:,:,:) :: image,illum
    real, dimension(:,:,:), allocatable :: tmpim,tmpil

    integer :: i,j,ix,iy,ix1,iy1
    real    :: coord_x, coord_y

    allocate(tmpim(mod%nz,mod%nx,mod%ny))
    allocate(tmpil(mod%nz,mod%nx,mod%ny))
    tmpim=0.
    tmpil=0.

    do j=1,mod%nyw
       coord_y=(j-1)*mod%dy+mod%oyw
       iy=nint((coord_y-mod%oy)/mod%dy)+1
       if ((iy.lt.1).or.(iy.gt.mod%ny)) cycle

       do i=1,mod%nxw
          coord_x=(i-1)*mod%dx+mod%oxw
          ix=nint((coord_x-mod%ox)/mod%dx)+1
          if ((ix.lt.1).or.(ix.gt.mod%nx)) cycle
          
          tmpim(1:mod%nz,ix,iy)=tmpim(1:mod%nz,ix,iy)+mod%imagesmall(1:mod%nz,i,j)
          tmpil(1:mod%nz,ix,iy)=tmpil(1:mod%nz,ix,iy)+mod%illumsmall(1:mod%nz,i,j)
       end do
    end do

    image=image+tmpim
    illum=illum+tmpil

    deallocate(tmpim,tmpil) 

  end subroutine mod_copy_image

end module ExtractPadModel_mod
