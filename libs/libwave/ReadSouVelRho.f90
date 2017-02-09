module Readsouvelrho_mod
  
  use sep

  use Interpolate_mod
  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  implicit none

contains

  logical function is_trace_within_model(data,mod)
    type(TraceSpace) ::                  data
    type(ModelSpace) ::                       mod
    real :: minb(3)
    real :: maxb(3)
    integer :: i
    is_trace_within_model=.true.

    minb(1)=mod%oz
    maxb(1)=mod%oz+mod%dz*(mod%nz-1)
    minb(2)=mod%oxw
    maxb(2)=mod%oxw+mod%dx*(mod%nxw-1)
    minb(3)=mod%oyw
    maxb(3)=mod%oyw+mod%dy*(mod%nyw-1)
    
    do i=1,3
       if ((data%coord(i).lt.minb(i)).or.(data%coord(i).gt.maxb(i))) then
          is_trace_within_model=.false.
          exit
       end if
    end do
  end function is_trace_within_model

  subroutine are_traces_within_model(datavec,mod)
    type(TraceSpace), dimension(:) ::datavec
    type(ModelSpace)               ::        mod
    integer                        :: i
                     
    do i=1,size(datavec)
       if(.not.is_trace_within_model(datavec(i),mod)) then
          write(0,*) 'ERROR: Trace(',i,') is not within model bounds'
          write(0,*) 'ERROR: Z=',datavec(i)%coord(1)
          write(0,*) 'ERROR: X=',datavec(i)%coord(2)
          write(0,*) 'ERROR: Y=',datavec(i)%coord(3)
          write(0,*) 'ERROR: Quit now'
          call erexit('ERROR: trace is not withing velocity model bounds')
       end if
    end do
  end subroutine are_traces_within_model

  subroutine readtraces(datavec,sourcevec,genpar)
    type(TraceSpace), dimension(:), allocatable :: datavec
    type(TraceSpace), dimension(:) :: sourcevec
    type(GeneralParam)             :: genpar
    integer                        :: ntraces,n1,i
    real                           :: d1,o1

    write(0,*) 'INFO:------------------------'
    write(0,*) 'INFO: Starting reading traces'

    if (.not.exist_file('traces')) call erexit('ERROR: need traces file')
    call from_aux('traces','n1',n1)
    call from_aux('traces','d1',d1)
    call from_aux('traces','o1',o1)
    call from_aux('traces','n2',ntraces)
    if (n1.ne.sourcevec(1)%dimt%nt) then
       call erexit('ERROR: nt traces and sources are different, exit now')
    end if
    allocate(datavec(ntraces))
    do i=1,ntraces
       allocate(datavec(i)%trace(n1,1))
       call sreed('traces',datavec(i)%trace(:,1),4*n1)
       datavec(i)%dimt%nt=n1
       datavec(i)%dimt%ot=o1
       datavec(i)%dimt%dt=d1     
    end do

    write(0,*) 'INFO: Finished reading traces'
    write(0,*) 'INFO: -----------------------'

  end subroutine readtraces

  subroutine readcoords(datavec,sourcevec,genpar)
    type(TraceSpace), dimension(:) :: datavec,sourcevec
    type(GeneralParam)             :: genpar
    integer :: index_gx,   index_gy,index_sx,index_sy
    integer :: index_selev,index_gelev
    integer :: nkeys,ntraces,i
    real, allocatable, dimension(:) :: tracekeys

    write(0,*) 'INFO:-------------------------------------------'
    write(0,*) 'INFO: Starting reading traces/source coordinates'

    call from_param('keygx',index_gx)
    call from_param('keysx',index_sx)
    call from_param('keygy',index_gy)
    call from_param('keysy',index_sy)
    call from_param('keyselev',index_selev)
    call from_param('keygelev',index_gelev)

    if (.not.exist_file('coordfile')) call erexit('ERROR: need coordfile')
    call from_aux('coordfile','n1',nkeys)
    call from_aux('coordfile','n2',ntraces)

    if(ntraces.ne.size(datavec)) then
       call erexit('ERROR: size datavec ne number of traces in coordfiles, exit')
    end if

    allocate(tracekeys(nkeys))
    do i=1,ntraces
       tracekeys=0.
       call sreed('coordfile',tracekeys,4*nkeys)      
       datavec(i)%coord(1)=tracekeys(index_gelev)    
       datavec(i)%coord(2)=tracekeys(index_gx)    
       datavec(i)%coord(3)=tracekeys(index_gy)
    end do

    call auxclose('coordfile')
    sourcevec(1)%coord(1)=tracekeys(index_selev)   
    sourcevec(1)%coord(2)=tracekeys(index_sx)    
    sourcevec(1)%coord(3)=tracekeys(index_sy)

    write(0,*) 'INFO: Source coordinates'
    write(0,*) 'INFO: ------------------'
    write(0,*) 'INFO: z=',sourcevec(1)%coord(1)
    write(0,*) 'INFO: x=',sourcevec(1)%coord(2)
    write(0,*) 'INFO: y=',sourcevec(1)%coord(3)
    write(0,*) 'INFO: '

    write(0,*) 'INFO: Finished reading traces/source coordinates'
    write(0,*) 'INFO:-------------------------------------------'
    deallocate(tracekeys)

  end subroutine readcoords

  subroutine readgathercoords(shotgath,sourcevec,genpar)
    type(GatherSpace), dimension(:), allocatable :: shotgath
    type(TraceSpace), dimension(:), allocatable :: sourcevec
    type(GeneralParam)             :: genpar
    integer :: index_gx,   index_gy,index_sx,index_sy
    integer :: index_selev,index_gelev
    integer :: nkeys,ntraces,i,j,k,ngathers
    integer :: begi,endi
    real,    allocatable, dimension(:)  :: tracekeys
    integer, allocatable, dimension(:)  :: sindex
    integer, allocatable, dimension(:)  :: bindex
    real,    allocatable, dimension(:,:):: coord

    write(0,*) 'INFO:-------------------------------------------'
    write(0,*) 'INFO: Starting reading traces/source coordinates'

    call from_param('keygx',index_gx)
    call from_param('keysx',index_sx)
    call from_param('keygy',index_gy)
    call from_param('keysy',index_sy)
    call from_param('keyselev',index_selev)
    call from_param('keygelev',index_gelev)

    if (.not.exist_file('coordfile')) call erexit('ERROR: need coordfile')
    call from_aux('coordfile','n1',nkeys)
    call from_aux('coordfile','n2',ntraces)
    
    allocate(coord(nkeys,ntraces))
    call sreed('coordfile',coord,4*nkeys*ntraces)

    ngathers=number_of_shots(index_sx,index_sy,coord)
    allocate(shotgath(ngathers))
    allocate(sourcevec(ngathers))
    allocate(sindex(ngathers))
    allocate(bindex(ngathers))
    call number_of_traces(index_sx,index_sy,ngathers,sindex,bindex,coord)

    do i=1,ngathers
       shotgath(i)%ntraces=sindex(i)
       shotgath(i)%begi=bindex(i)
       allocate(shotgath(i)%gathtrace(sindex(i)))
    end do
    
    allocate(tracekeys(nkeys))
    begi=0
    endi=0
    k=0

    ! Now read the coordinates for each gathers
    write(0,*) 'INFO: ngathers',ngathers,'nkeys',nkeys,'ntraces',ntraces
    do i=1,ngathers
       begi=endi+1
       endi=begi+shotgath(i)%ntraces-1
       k=0
       do j=begi,endi
          tracekeys=0.  
          tracekeys=coord(:,j)
          k=k+1
          shotgath(i)%gathtrace(k)%coord(1)=tracekeys(index_gelev)    
          shotgath(i)%gathtrace(k)%coord(2)=tracekeys(index_gx)    
          shotgath(i)%gathtrace(k)%coord(3)=tracekeys(index_gy)
       end do
       sourcevec(i)%coord(1)=tracekeys(index_selev)   
       sourcevec(i)%coord(2)=tracekeys(index_sx)    
       sourcevec(i)%coord(3)=tracekeys(index_sy)

       if(genpar%verbose) write(0,*) 'INFO: Source coordinates for gather',i,'with ',shotgath(i)%ntraces,'traces are'
       if(genpar%verbose) write(0,*) 'INFO: ------------------------------------------------------------------------'
       if(genpar%verbose) write(0,*) 'INFO: x,y,z=',sourcevec(i)%coord(2),sourcevec(i)%coord(3),sourcevec(i)%coord(1)
       if(genpar%verbose) write(0,*) 'INFO: '
    end do

    write(0,*) 'INFO: Finished reading traces/source coordinates'
    write(0,*) 'INFO:-------------------------------------------'

    deallocate(tracekeys,sindex,bindex,coord)
  end subroutine readgathercoords

  subroutine readgathertraces(shotgath,source)
    type(GatherSpace), dimension(:):: shotgath
    type(TraceSpace),  dimension(:):: source
    integer                        :: nshots,ntraces,nfiletraces,nreadtraces,n1,j,i
    real                           :: d1,o1

    if (.not.exist_file('traces')) call erexit('ERROR: need traces file')
    call from_aux('traces','n1',n1)
    call from_aux('traces','d1',d1)
    call from_aux('traces','o1',o1)
    call from_aux('traces','n2',nfiletraces)
    if (n1.ne.source(1)%dimt%nt) then
       call erexit('ERROR: nt traces and sources are different, exit now')
    end if

    nreadtraces=0
    nshots=size(shotgath)
    write(0,*) 'INFO:-------------------------------------------'
    do i=1,nshots

       ntraces=shotgath(i)%ntraces
       if (modulo(i,10).eq.0) write(0,*) 'INFO: Starting reading',ntraces,'traces for shot',i

       do j=1,ntraces
          nreadtraces=nreadtraces+1
          allocate(shotgath(i)%gathtrace(j)%trace(n1,1))
          call sreed('traces',shotgath(i)%gathtrace(j)%trace(:,1),4*n1)
          shotgath(i)%gathtrace(j)%dimt%nt=n1
          shotgath(i)%gathtrace(j)%dimt%ot=o1
          shotgath(i)%gathtrace(j)%dimt%dt=d1     
       end do
    end do
    write(0,*) 'INFO:-------------------------------------------'

    write(0,*) 'INFO: Read a total of ',nreadtraces
    if(nreadtraces.ne.nfiletraces) call erexit('ERROR: number of read traces is different from n2 in traces')

    write(0,*) 'INFO: Finished reading traces'
    write(0,*) 'INFO: -----------------------'

  end subroutine readgathertraces

  subroutine readsou(source,genpar)
    type(TraceSpace), dimension(:), allocatable ::source
    type(GeneralParam) ::   genpar
    real            ::sinc(genpar%lsinc)
    real            ::fdum,fscale
    real,dimension(:)  ,allocatable :: trace
    integer         :: it,l

    allocate(source(1))
    call from_history('n1',source(1)%dimt%nt)
    call from_history('d1',source(1)%dimt%dt)
    call from_history('o1',source(1)%dimt%ot)
    allocate(source(1)%trace(source(1)%dimt%nt,1))
    call sreed('in',source(1)%trace(:,1),4*source(1)%dimt%nt)
    
    genpar%dt=source(1)%dimt%dt
    genpar%nt=source(1)%dimt%nt
    genpar%t0=source(1)%dimt%ot
    
    if (genpar%Born.or.genpar%LSRTM) then
       allocate(trace(source(1)%dimt%nt))
       trace=source(1)%trace(:,1)
       source(1)%trace(:,1)=0
       do it=2,source(1)%dimt%nt-1
          source(1)%trace(it,1)=source(1)%trace(it,1)+(-2*trace(it)+trace(it-1)+trace(it+1))/(source(1)%dimt%dt**2)
       end do
       deallocate(trace)
    end if

    if (genpar%withRho) then
       allocate(trace(source(1)%dimt%nt))
       ! Shift by 1/2 time step forward(reason: staggered grid)
       call mksinc(sinc,genpar%lsinc,+0.5)
       trace = source(1)%trace(1:source(1)%dimt%nt,1)
       do it=1,source(1)%dimt%nt
          fscale = 1.
          fdum = 0.
          do l=-genpar%lsinc/2,genpar%lsinc/2
             if (it+l.ge.1 .and. it+l.le.source(1)%dimt%nt) then
                fdum = fdum + sinc(genpar%lsinc/2+1+l)*trace(it+l)
             else
                fscale = fscale - sinc(genpar%lsinc/2+1+l)
             endif
          end do
          source(1)%trace(it,1)= fdum / fscale
       end do
       deallocate(trace)
    end if

  end subroutine readsou

  subroutine copysou2sougath(sou,sougath)
    type(TraceSpace), dimension(:) ::sou,sougath
    integer :: i,nt
    nt=sou(1)%dimt%nt
    do i=1,size(sougath)
       allocate(sougath(i)%trace(nt,1))
       sougath(i)%trace=sou(1)%trace
       sougath(i)%dimt%nt=sou(1)%dimt%nt
       sougath(i)%dimt%dt=sou(1)%dimt%dt
    end do
  end subroutine copysou2sougath

  subroutine readsoucoord(source,mod)    
    type(TraceSpace), dimension(:) :: source
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

    call from_param('source_z',source(1)%coord(1),0.)
    call from_param('source_x',source(1)%coord(2),mid(2))
    call from_param('source_y',source(1)%coord(3),mid(3))

    do i=1,3
       if ((source(1)%coord(i).lt.mino(i)).or.(source(1)%coord(i).gt.maxo(i))) then
          write(0,*) 'ERROR:  Problem with axis',i
          call erexit('ERROR: source_coord is not within model bounds, exit now')
       end if
    end do

    write(0,*) 'source coordinates =',source(1)%coord

  end subroutine readsoucoord

  subroutine readvel(    mod,genpar,bounds)
    type(ModelSpace)  :: mod
    type(GeneralParam)::     genpar
    type(FDbounds)    ::            bounds
    real,allocatable,dimension(:,:,:) :: tmp

    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    
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

    if (genpar%shot_type.gt.0 .or. genpar%surf_type.gt.0) then
       if (genpar%shot_type.gt.0) bounds%nmin1 = -(genpar%ntaper+genpar%lsinc/2+2)+1
       if (genpar%surf_type.gt.0) bounds%nmin1 = -genpar%lsinc+1
    else
       bounds%nmin1 = -genpar%ntaper+1
    endif
    bounds%nmax1 =  mod%nz+genpar%ntaper
    bounds%nmin2 = -genpar%ntaper+1
    bounds%nmax2 =  mod%nx+genpar%ntaper
    if (.not. genpar%twoD) then
       bounds%nmin3 = -genpar%ntaper+1
       bounds%nmax3 =  mod%ny+genpar%ntaper
    else
       bounds%nmin3 = 1
       bounds%nmax3 = 1
    end if
  
    mod%oz=genpar%omodel(1)
    mod%ox=genpar%omodel(2)
    mod%oy=genpar%omodel(3)

    mod%oxw=genpar%omodel(2)
    mod%oyw=genpar%omodel(3)

    mod%dz=genpar%delta(1)
    mod%dx=genpar%delta(2)
    mod%dy=genpar%delta(3)

    mod%nxw=mod%nx
    mod%nyw=mod%ny

    allocate(tmp(mod%nz,mod%nx,mod%ny))
    allocate(mod%vel(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))

    tmp=0.
    if (.not.exist_file(mod%veltag)) then
       call erexit('ERROR: Need velocity file, exit now')
    else
       call sreed(mod%veltag,tmp,4*mod%nz*mod%nx*mod%ny)
       call auxclose(mod%veltag)
       call vel_check(tmp,genpar)
       call model_pad(tmp,mod%vel,bounds,mod%nz,mod%nx,mod%ny,genpar)
    end if

    call srite('tmpvel',mod%vel,4*(bounds%nmax1-bounds%nmin1+1)*(bounds%nmax2-bounds%nmin2+1)*(bounds%nmax3-bounds%nmin3+1))
    call to_history('n1',bounds%nmax1-bounds%nmin1+1,'tmpvel')
    call to_history('n2',bounds%nmax2-bounds%nmin2+1,'tmpvel')
    call to_history('n3',bounds%nmax3-bounds%nmin3+1,'tmpvel')
    if(genpar%withRho) then
       if (.not.exist_file(mod%rhotag)) then
          call erexit('ERROR: Need rho file, exit now')
       else
          tmp=0.
          call dim_consistency_check(mod%rhotag,mod%veltag,genpar%twoD)
          allocate(mod%rho(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
          allocate(mod%rho2(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
          call sreed(mod%rhotag,tmp,4*mod%nz*mod%nx*mod%ny)
          call auxclose(mod%rhotag)
          call model_pad(tmp,mod%rho,bounds,mod%nz,mod%nx,mod%ny,genpar)
          call Interpolate(mod,bounds,genpar)
       end if
    end if
    deallocate(tmp)
  end subroutine readvel

  subroutine readref(         mod,genpar,bounds)
    type(ModelSpace)  ::      mod
    type(GeneralParam)::          genpar
    type(FDbounds)    ::                 bounds
    real,allocatable,dimension(:,:,:) :: tmp

    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar
    
    allocate(tmp(mod%nz,mod%nx,mod%ny))

    if (.not.exist_file(mod%reftag)) then
       call erexit('ERROR: Need reflectivity file, exit now')
    else
       tmp=0.
       call dim_consistency_check(mod%veltag,mod%reftag,genpar%twoD)
       allocate(mod%image(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))
       call sreed(mod%reftag,tmp,4*mod%nz*mod%nx*mod%ny)
       call auxclose(mod%reftag)
       call model_pad(tmp,mod%image,bounds,mod%nz,mod%nx,mod%ny,genpar)
    end if
    deallocate(tmp)
  end subroutine readref
  
  subroutine model_pad(tmp,field,bounds,nz,nx,ny,genpar)
    type(GeneralParam)     ::                    genpar
    type(FDbounds)         :: bounds
    real, dimension(:,:,:) :: tmp
    real                   :: field(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4, bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound)
    integer :: i,j,k,nx,ny,nz
    logical :: twoD

    field(1:nz,1:nx,1:ny)=tmp
    
    do k=1,ny
       do j=1,nx           
          do i=bounds%nmin1-4,0
             field(i,j,k)=field(1,j,k)
          end do         
          do i=nz+1,bounds%nmax1+4
             field(i,j,k)=field(nz,j,k)
          end do
       end do
    end do

    do k=1,ny
       do j=bounds%nmin2-4,0
          do i=bounds%nmin1-4,bounds%nmax1+4
             field(i,j,k)=field(i,1,k)
          end do         
       end do
       do j=nx+1,bounds%nmax2+4
          do i=bounds%nmin1-4,bounds%nmax1+4
             field(i,j,k)=field(i,nx,k)
          end do         
       end do
    end do

    if (.not.genpar%twoD) then
       do k=bounds%nmin3-genpar%nbound,0
          do j=bounds%nmin2-4,bounds%nmax2+4
             do i=bounds%nmin1-4,bounds%nmax1+4
                field(i,j,k)=field(i,j,1)
             end do
          end do
       end do
       do k=ny+1,bounds%nmax3+genpar%nbound
          do j=bounds%nmin2-4,bounds%nmax2+4
             do i=bounds%nmin1-4,bounds%nmax1+4
                field(i,j,k)=field(i,j,ny)
             end do
          end do
       end do
       
    end if
  end subroutine model_pad

  !
  ! Subroutine to check the parameters supplied for stability and accuracy
  ! If the computations are too inaccurate or unstable, then the program
  ! will exit
  !
  subroutine vel_check(vel,genpar)
    real, dimension(:,:,:) :: vel
    type(GeneralParam) ::genpar

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
    real :: lambdamin, fmax_new, stab_cond
    !
    ! check the velocity model for its maximum value
    !
    stab_cond=0.35
    stab_cond=0.45
    if (genpar%twoD) stab_cond=0.45

    vmax=maxval(vel)
    vmin=minval(vel)

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
    if (stab.gt.stab_cond) then
       if (genpar%twoD) then
          call putlin(							&
          &      'The Courant number is >0.45 so calculations for this')
       else
          call putlin(							&
          &      'The Courant number is >0.35 so calculations for this')
       end if
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
       dtn=stab_cond/stab*genpar%dt
       !
       call putch('dtn = new_dt','r',dtn)
       call putlin('     ')
       call erexit('ERROR: computations will be unstable')
       !
    else
       if (genpar%twoD) then
          call putlin(							&
          &      'The Courant number is <= 0.45, calculations will procede')
          call putlin(' ')
       else
          call putlin(							&
          &      'The Courant number is <= 0.35, calculations will procede')
          call putlin(' ')
       end if
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
       fmax_new=samplerate/3.0*genpar%fmax
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

  subroutine dim_consistency_check(tag1,tag2,twoD)
    character(len=3) :: tag1,tag2
    integer  :: fetch, hetch, tetch, getch, auxpar
    external :: fetch, hetch, tetch, getch, auxpar

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

  end subroutine dim_consistency_check

  integer function number_of_shots(index_sx,index_sy,coord)
    integer :: i,nkeys,index_sx,index_sy,ntraces
    real    :: sx,sy,sxtmp,sytmp
    real, allocatable, dimension(:)   :: tracekeys
    real,              dimension(:,:) :: coord

    ntraces=size(coord(1,:))
    nkeys=size(coord(:,1))
    allocate(tracekeys(nkeys))

    tracekeys=coord(:,1)
    sx=tracekeys(index_sx)
    sy=tracekeys(index_sy)
    number_of_shots=1
    do i=2,ntraces
       
       tracekeys=coord(:,i)
       sxtmp=tracekeys(index_sx)
       sytmp=tracekeys(index_sy)
       if ((sx.ne.sxtmp).or.(sy.ne.sytmp)) then
          number_of_shots=number_of_shots+1
          sx=sxtmp
          sy=sytmp
       end if
    end do

    deallocate(tracekeys)
  end function number_of_shots

  subroutine number_of_traces(index_sx,index_sy,nshots,sindex,bindex,coord)
    integer :: i,j,nshots,nkeys,index_sx,index_sy,ntraces
    real    :: sx,sy,sxtmp,sytmp
    real, allocatable, dimension(:)        :: tracekeys
    integer,           dimension(nshots)   :: sindex
    integer,           dimension(nshots)   :: bindex
    real,              dimension(:,:)      :: coord

    sindex=0
    ntraces=size(coord(1,:))
    nkeys=size(coord(:,1))
    allocate(tracekeys(nkeys))

    tracekeys=coord(:,1)
    sx=tracekeys(index_sx)
    sy=tracekeys(index_sy)
    j=1
    sindex(j)=sindex(j)+1
    bindex(j)=1
    do i=2,ntraces
       tracekeys=coord(:,i)
       sxtmp=tracekeys(index_sx)
       sytmp=tracekeys(index_sy)
       if ((sx.ne.sxtmp).or.(sy.ne.sytmp)) then
         j=j+1
         sx=sxtmp
         sy=sytmp
         bindex(j)=i
       endif
       sindex(j)=sindex(j)+1
    end do

    deallocate(tracekeys)
  end subroutine number_of_traces

end module readsouvelrho_mod
