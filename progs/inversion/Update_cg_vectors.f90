! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Update_cg_vectors

  use sep
  use ReadData_mod
  use DataSpace_types_flt

  type(cube)       :: dat1,dat2 
  double precision :: alfa,beta
  integer          :: ndim,iter

  call sep_init(SOURCE)
  ndim=sep_dimension('file1')

  call from_param('iter',iter)

  write(0,*) 'INFO: iter=',iter

  call ReadData_dim('file1',dat1,ndim)
  call ReadData_dim('file2',dat2,ndim)

  if (.not.ns_are_consistent(ndim,dat1,dat2)) call erexit('ERROR: dimensions file1 and file2 are not the same, exist now')

  call ReadData_cube('file1',dat1)
  call auxclose('file1')
  call ReadData_cube('file2',dat2)
  call auxclose('file2')

  call sreed('alfabeta',alfa,8)
  call sreed('alfabeta',beta,8)

  write(0,*) 'INFO: updating file1=',beta,'*file1+',alfa,'*file2'
  dat1%dat=beta*dat1%dat+alfa*dat2%dat

  call WriteData_cube('file1',dat1)
  call WriteData_dim('file1',dat1,ndim)

end program Update_cg_vectors
