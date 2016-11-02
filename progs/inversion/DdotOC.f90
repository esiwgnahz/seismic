program Ddot_OC

  use sep
  use occgmod_mod
  use ReadData_mod
  use DataSpace_types_flt

  implicit none

  type(cube)       :: dat1,dat2
  double precision :: ddot
  logical          :: add
  integer :: ndim

  call sep_init(SOURCE)
  ndim=sep_dimension('file1')

  call ReadData_dim('file1',dat1,ndim)
  call ReadData_dim('file2',dat2,ndim)
  if (.not.ns_are_consistent(ndim,dat1,dat2)) call erexit('ERROR: dimensions file1 and file2 are not the same, exist now')

  call ReadData_cube('file1',dat1)
  call auxclose('file1')
  call ReadData_cube('file2',dat2)
  call auxclose('file2')

  call from_param('add',add,.false.)

  write(0,*) 'INFO: add=',add

  if(.not.exist_file('dot_product')) call erexit('ERROR: need dot_product file to save value')

  if (add) then
     call sreed('dot_product',ddot,8)
     call auxclose('dot_product')
  else
     ddot=0.d0
     call to_history('n1',1,'dot_product')
     call to_history('esize',8,'dot_product')
  end if

  write(0,*) 'INFO: ddot in',ddot
  ddot=ddot+dot_product(dat1%dat,dat2%dat)
  write(0,*) 'INFO: ddot out',ddot

  call srite('dot_product',ddot,8)

end program Ddot_OC
