! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Mute_gather

  use sep
  use DataSpace_types
  use GeneralParam_types

  implicit none

  type MuteParam
     real :: v0
     real :: t0
     real :: v1
     real :: t1
     real :: tramp

     type(GatherSpace), dimension(:), allocatable :: maskgath

  end type MuteParam

contains

  subroutine Init_MuteParam(mutepar)
    type(MuteParam) ::      mutepar

    call from_param('v0',mutepar%v0,1500.)
    call from_param('v1',mutepar%v1,10000.)
    call from_param('t0',mutepar%t0,0.)
    call from_param('t1',mutepar%t1,40.)
    call from_param('tramp',mutepar%tramp,0.5)

    write(0,*) 'INFO:-------------------------'
    write(0,*) 'INFO: Data muting parameters  '
    write(0,*) 'INFO:-------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  V0    = ',mutepar%v0
    write(0,*) 'INFO:  V1    = ',mutepar%v1
    write(0,*) 'INFO:  t0    = ',mutepar%t0
    write(0,*) 'INFO:  t1    = ',mutepar%t1
    write(0,*) 'INFO:  tramp = ',mutepar%tramp
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-------------------------'

  end subroutine Init_MuteParam

  subroutine MuteParam_compute_mask(mutepar,shotgath,sourcegath)
    type(MuteParam)    ::           mutepar
    type(GatherSpace), dimension(:) ::      shotgath
    type(TraceSpace),  dimension(:) ::               sourcegath

    real               :: hx,hy,sx,sy,gx,gy,dt,t0
    integer            :: i,j,nt

    real, dimension(:), allocatable :: tmp, tmp1

    dt=shotgath(1)%gathtrace(1)%dimt%dt
    t0=shotgath(1)%gathtrace(1)%dimt%ot
    nt=shotgath(1)%gathtrace(1)%dimt%nt

    allocate(tmp(nt),tmp1(nt))

    allocate(mutepar%maskgath(size(shotgath)))


    do i=1,size(shotgath)
       
       sx=sourcegath(i)%coord(2)
       sy=sourcegath(i)%coord(3)

       allocate(mutepar%maskgath(i)%gathtrace(shotgath(i)%ntraces))

       mutepar%maskgath(i)%ntraces=shotgath(i)%ntraces
       
       do j=1,shotgath(i)%ntraces

          allocate(mutepar%maskgath(i)%gathtrace(j)%trace(nt,1))

          gx=shotgath(i)%gathtrace(j)%coord(2)
          gy=shotgath(i)%gathtrace(j)%coord(3)

          hx=gx-sx
          hy=gy-sy
          
          call compute_top_mute(mutepar,dt,t0,hx,hy,tmp)
          call compute_bottom_mute(mutepar,dt,t0,hx,hy,tmp1)
          mutepar%maskgath(i)%gathtrace(j)%trace(:,1)=tmp*tmp1
          if (exist_file('mute_out')) call srite('mute',mutepar%maskgath(i)%gathtrace(j)%trace,4*nt)
       end do

    end do
    
    deallocate(tmp,tmp1)
 
  end subroutine MuteParam_compute_mask

  subroutine MuteParam_deallocate(mutepar)
    type(MuteParam) ::            mutepar
    integer :: i,j
    
    do j=1,size(mutepar%maskgath)
       do i=1,mutepar%maskgath(j)%ntraces
          deallocate(mutepar%maskgath(j)%gathtrace(i)%trace)
       end do
       deallocate(mutepar%maskgath(j)%gathtrace)
    end do
    deallocate(mutepar%maskgath)

  end subroutine MuteParam_deallocate

  subroutine compute_top_mute(mutepar,dt,t0,hx,hy,mute)
    type(MuteParam)    ::       mutepar
    real, dimension(:) ::                           mute
    real               ::               dt,t0,hx,hy

    integer            :: i,it,nt,nramp,beg,end
    real               :: t_init,tmute, pi=3.1415926535

    nt=size(mute)
    mute=0

    nramp=mutepar%tramp/dt+1.5
    t_init=mutepar%t0
    tmute=sqrt(t_init**2+((hx+hy)/mutepar%v0)**2)
    it=(tmute-t0)/dt+1.5

    beg=max(1,it-nramp)
    end=min(nt,it)
    do i=beg,end
       mute(i)=sin((i-it+nramp)*pi/(2*nramp))**2
    end do
    mute(it+1:)=1
    
  end subroutine compute_top_mute
    
  subroutine compute_bottom_mute(mutepar,dt,t0,hx,hy,mute)
    type(MuteParam)    ::       mutepar
    real, dimension(:) ::                       mute
    real               ::                 dt,t0,hx,hy

    integer            :: i,it,nt,nramp,beg,end
    real               :: t_init,tmute, pi=3.1415926535

    nt=size(mute)
    mute=1

    nramp=mutepar%tramp/dt+1.5
    t_init=mutepar%t1
    tmute=sqrt(t_init**2+((hx+hy)/mutepar%v1)**2)
    it=(tmute-t0)/dt+1.5

    beg=max(1,it+1)
    end=min(nt,it+1+nramp)
    do i=beg,end
       mute(i)=cos((i-it-1)*pi/(2*nramp))**2
    end do
    mute(it+1+nramp:)=0
    
  end subroutine compute_bottom_mute
    
    
end module Mute_gather
