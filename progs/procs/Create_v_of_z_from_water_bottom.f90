! Hang a 1D velocity profile from a picked surface
! v(z)=v0+coef*z^pow, where z is the distance from the picked surface
! Above the picked surface, the velocity is given by vwater.

program Create_v_of_z

  use sep

  implicit none

  real,    dimension(:), allocatable:: in,ou

  integer, dimension(:), allocatable:: n
  real,    dimension(:), allocatable:: o
  real,    dimension(:), allocatable:: d
  character(len=1024)               :: label

  integer :: ndim,i,j,i1,i2,i3,i4,i5,i6,nt
  real    :: v0,pow,coef,z0,z,z1,z_ref,vwater,dvc,max_depth,dr,max_z
  logical :: found1,found2,found3,found4,found5,found6,velo
  real :: pi=3.141592653589793238462643

  call sep_init(SOURCE)
  ndim=sep_dimension()
  allocate(n(ndim), o(ndim), d(ndim))
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do

  allocate(in(n(1)),ou(n(1)))

  call from_param("vwater",vwater,1500.)
  call from_param("v0",v0,vwater)
  call from_param("pow",pow,1.)
  call from_param("coef",coef,2.)

  max_depth=o(1)+(n(1)-1)*d(1)
  call from_param("depth_v_constant",dvc,max_depth)
  call from_param("depth_ref",dr,dvc/2)
  call from_param("max_depth_from_reflector",max_z,max_depth)
  call from_param("velo",velo,.true.)
  call from_param("nt",nt,2) ! minimum distance between itop and ibot

  do j=1,n(2)

     call sep_read(in)
     found1=.false.
     found2=.false.
     found3=.false.
     found4=.false.
     found5=.false.
     found6=.false.

     z=0.
     z1=0.

     do i=1,n(1)
        z0=o(1)+(i-1)*d(1)

        if (found1) then
           if (found2) then
              if (found3) then
                 if (found4) then
                    if (found5) then
                       if (.not.found6) then
                          if (in(i).ne.in(i-1)) then
                             if (in(i).ne.0) then
                                found6=.true.
                                i6=i
                             end if
                          end if
                       end if
                    end if
                 end if
              end if
           end if
        end if
        if (found1) then
           if (found2) then
              if (found3) then
                 if (found4) then
                    if (.not.found5) then
                       if (in(i).ne.in(i-1)) then
                          if (in(i).ne.0) then
                             found5=.true.
                             i5=i
                          end if
                       end if
                    end if
                 end if
              end if
           end if
        end if
        if (found1) then
           if (found2) then
              if (found3) then
                 if (.not.found4) then
                    if (in(i).ne.in(i-1)) then
                       if (in(i).ne.0) then
                          found4=.true.
                          i4=i
                       end if
                    end if
                 end if
              end if
           end if
        end if
        if (found1) then
           if (found2) then
              if (.not.found3) then
                 if (in(i).ne.in(i-1)) then
                    if (in(i).ne.0) then
                       found3=.true.
                       i3=i
                    end if
                 end if
              end if
           end if
        end if
        if (found1) then
           if (.not.found2) then
              if (in(i).ne.in(i-1)) then
                 if (in(i).ne.0) then
                    found2=.true.
                    i2=i
                 end if
              end if
           end if
        end if
        if(.not.found1) then
           if (in(i).ne.0) then
              found1=.true.
              i1=i
              z_ref=z0
           end if
        end if

        if (found1) then 
           if (.not.found2) then
              z=z+d(1)
              if (velo) then
                 if (z0.le.dvc) then
                    ou(i)=v0+coef*z**pow
                 else 
                    z1=z0-dr
                    ou(i)=v0+coef*z1**pow
                 end if
              else
                 ou(i)=v0+v0*(sin((z0-z_ref)*pi/(2*max_z)))**2
              end if
              if (z.gt.max_z+d(1)) then
                 ou(i)=ou(i-1)
              end if
           end if
        else
           ou(i)=vwater
        end if
     end do
     if (found2) ou(i2:)=vwater
     if (found3) ou(i3:)=v0
     if (found4) ou(i4:)=vwater
     if (found5) ou(i5:)=v0
     if (found6) ou(i6:)=vwater
     call sep_write(ou)
!     if (found1.and.found2) then
!        write(0,*) '1',ou(itop:ibot)
!        write(0,*) '2',ou(itop:)
!     end if
  end do

end program Create_v_of_z
