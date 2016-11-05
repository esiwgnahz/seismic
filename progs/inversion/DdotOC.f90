program Ddot_OC

  use sep
  use occgmod_mod
  use ReadData_mod
  use DataSpace_types_flt

  implicit none

  type(cube)       :: dat1,dat2
  double precision :: ddotd
  real             :: ddots
  logical          :: add, double_precision
  integer :: ndim,esize

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
  call from_param('double_precision',double_precision,.false.)

  write(0,*) 'INFO: add=',add
  write(0,*) 'INFO: double precision',double_precision

  if(.not.exist_file('dot_product')) call erexit('ERROR: need dot_product file to save value')

  if (add) then
     if (double_precision) then
        call sreed('dot_product',ddotd,8)
     else
        call sreed('dot_product',ddots,4)
     end if
     call auxclose('dot_product')
  else
     ddotd=0.d0
     ddots=0.
  end if

  if (double_precision) then
     write(0,*) 'INFO: ddot in',ddotd
     ddotd=ddotd+dot_product(dat1%dat,dat2%dat)
     write(0,*) 'INFO: ddot out',ddotd
  else
     write(0,*) 'INFO: ddot in',ddots
     ddotd=dot_product(dat1%dat,dat2%dat)
     ddots=ddots+sngl(ddotd)
     write(0,*) 'INFO: ddot out',ddots
  end if

  esize=4
  if (double_precision) then
     esize=8
     call srite('dot_product',ddotd,8)
  else
     call srite('dot_product',ddots,4)
  end if
  call to_history('n1',1,'dot_product')
  call to_history('esize',esize,'dot_product')

end program Ddot_OC
