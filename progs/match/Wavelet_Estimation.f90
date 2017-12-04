program Wavelet_Estimation

  use sep
  use cgstep_mod

  use identity_mod
  use nchconest
  use nhelicon_mod
  use solver_smp_mod
  use BuildAdaptiveFilter_mod

  implicit none

  integer, external                    :: getch
  real,    dimension(:), allocatable   :: synt_data,obs_data,matched_data,weight
  real,    dimension(:), allocatable   :: in_wavelet, ou_wavelet
  type (filter)                        :: match_n,laplac_n, match_tmp

  integer :: ndim                        ! dimension of dataset
  integer :: ndim_domain                 ! dimension for the filtering (will be 2 here)
  integer :: niter

  character (len=20) :: targetname
  integer, dimension(:), pointer :: n                  ! data-size, dimension(ndim)
  real,    dimension(:), pointer :: d, o
  integer, dimension(:), pointer :: n_obs                  ! data-size, dimension(ndim)
  real,    dimension(:), pointer :: d_obs, o_obs

  integer :: i,stat,ilag,ipatch,idim,ipanels,npanels   ! junk variables
  integer :: n12
  integer :: verb
  character(len=1024)    :: label

  integer, dimension(:), allocatable :: status,nfilt_n

  call sep_init()
  call Wavelet_estimation_doc()

  write(0,*) 'INFO:'
  write(0,*) 'INFO:'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO: Wavelet estimation with time-domain '
  write(0,*) 'INFO:           matching filters          '
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:'

  ndim=sep_dimension("modeled_data")
  ndim_domain=2      ! Filters work on a 2D plane only
  allocate(nfilt_n(1))

  allocate(n(ndim), o(ndim), d(ndim))
  allocate(n_obs(ndim), o_obs(ndim), d_obs(ndim))
  do i=1,ndim
     call sep_get_data_axis_par("modeled_data",i,n(i),o(i),d(i),label)
     call sep_get_data_axis_par("observed_data",i,n_obs(i),o_obs(i),d_obs(i),label)
     if ((n(i).ne.n_obs(i))) &
     & call seperr('ERROR: n,d or o are not equal in modeled and observed data, quit now')
  end do
  if (exist_file("weight")) then
     do i=1,ndim
        call sep_get_data_axis_par("weight",i,n_obs(i),o_obs(i),d_obs(i),label)
        if ((n(i).ne.n_obs(i))) &
        & call seperr('ERROR: n,d or o are not equal in weigth and modeled data, quit now')    
     end do
  end if
  deallocate(n_obs,o_obs,d_obs)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:       Dimensions Observed Data      '
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:'

  do i=1,ndim
     write(0,*) 'INFO: - - - - - - - '
     write(0,*) 'INFO: n(',i,')=',n(i)
     write(0,*) 'INFO: o(',i,')=',o(i)
     write(0,*) 'INFO: d(',i,')=',d(i)
  end do

  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:'

  n12=n(1)*n(2)

  ! Read synt_data filter paramaters
  call from_param('nfilt',nfilt_n,(/21/))
  call from_param('niter',niter,product(nfilt_n)-1)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:        Inversion parameters         '
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO: Number of iterations        =',niter
  write(0,*) 'INFO: Size of filter in # samples =',nfilt_n(1)
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:'

  allocate(synt_data(n12),obs_data(n12),matched_data(n12),in_wavelet(n(1)),ou_wavelet(n(1)))

  call create_match_filter(nfilt_n,match_tmp,1,n)
     
  write(0,*) 'INFO:'
  write(0,*) 'INFO: --------------------'
  write(0,*) 'INFO: Starting processing '
  write(0,*) 'INFO: --------------------'
  write(0,*) 'INFO:'

  obs_data=0;synt_data=0;in_wavelet=0.;ou_wavelet=0.;matched_data=0.
  
  call sreed("modeled_data",synt_data,n12*4)
  call sreed("observed_data",obs_data,n12*4)
  call sreed("input_wavelet",in_wavelet,4*n(1))
  
  if(exist_file("weight")) then
     allocate(weight(n12))
     call sreed("weight",weight,n12*4)
  end if

  write(0,*) 'INFO: Finished reading model, observed and wavelet'

  ! PROCESSING PART
  call allocatehelix(match_n,product(nfilt_n))
  match_n%lag=match_tmp%lag

  if (verb>2) write(0,*) 'INFO:  -----------------------------------'
  if (verb>2) write(0,*) 'INFO:  Calling findmatch with no weight...'
  if (verb>2) write(0,*) 'INFO:  -----------------------------------'
  
  call nchconest_init(synt_data,match_n)
  if (exist_file("weight")) then
     call weightm_init(weight)
     call solver_smp(m=match_n%flt,d=obs_data,Fop=nchconest_lop,Wop=weightm_lop,stepper=cgstep,niter=niter,verb=.true.)
     call weightm_close()
  else
     call solver_smp(m=match_n%flt,d=obs_data,Fop=nchconest_lop,stepper=cgstep,niter=niter,verb=.true.)
  end if
  call nchconest_close()

  call nchelicon_init(match_n)
  stat=nchelicon_lop(.false.,.false.,in_wavelet,ou_wavelet)
  stat=nchelicon_lop(.false.,.false.,synt_data,matched_data)
  call srite("output_wavelet",ou_wavelet,n(1)*4)
  call srite("matched_data",matched_data,n12*4)
 
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO: Done with processing'
  write(0,*) 'INFO: ------------------------------------'

  call to_history('n1',n(1),'output_wavelet')
  call to_history('d1',d(1),'output_wavelet')
  call to_history('o1',o(1),'output_wavelet')
  
  do i=1,ndim
     call sep_put_data_axis_par("matched_data",i,n(i),o(i),d(i),label)
  end do

end program WAVELET_ESTIMATION

subroutine Wavelet_estimation_doc()
  call sep_add_doc_line("NAME")
  call sep_add_doc_line("    Wavelet_Estimation.x - Wavelet estimation for FWI ")
  call sep_add_doc_line("")
  call sep_add_doc_line("SYNOPSIS")
  call sep_add_doc_line("    Wavelet_Estimation input_wavelet=inw.H output_wavelet=ouw.H modeled_data=mod.H observed_data=obs.H > /dev/null ")
  call sep_add_doc_line("")
  call sep_add_doc_line("INPUT FILE")
  call sep_add_doc_line("")
  call sep_add_doc_line("    modeled_data - sepfile")
  call sep_add_doc_line("             mod.H (nt,ntraces) = input")
  call sep_add_doc_line("    observed_data - sepfile")
  call sep_add_doc_line("             obs.H (nt,ntraces) = input")
  call sep_add_doc_line("    input_wavelet - sepfile")
  call sep_add_doc_line("             inw.H (n)          = input")
  call sep_add_doc_line("")
  call sep_add_doc_line("INPUT PARAMETERS")
  call sep_add_doc_line("")
  call sep_add_doc_line("    niter - int")
  call sep_add_doc_line("               Number of iterations for 1D filter estimation")
  call sep_add_doc_line("")
  call sep_add_doc_line("    nfilt - int")
  call sep_add_doc_line("               Number of filter coefficients (time only)")
  call sep_add_doc_line("")
  call sep_add_doc_line("OPTIONAL INPUT FILE")
  call sep_add_doc_line("    weight - sepfile")
  call sep_add_doc_line("")   
  call sep_add_doc_line("OPTIONAL OUTPUT FILE")
  call sep_add_doc_line("    matched_data - sepfile")
  call sep_add_doc_line("             match.H (nt,ntraces) = output")
  call sep_add_doc_line("")   
  call sep_add_doc_line("OUTPUT PARAMETERS")
  call sep_add_doc_line("    output_wavelet - sepfile")
  call sep_add_doc_line("            img.H (nz,nx,ny,2) = output")
  call sep_add_doc_line("            n4=1 f1=0 RTM image, n4=1 f1=1 illumination")
  call sep_add_doc_line("")
  call sep_add_doc_line("DESCRIPTION")
  call sep_add_doc_line("      Compute an improved wavelet to match the modeled data to the observed data solving")
  call sep_add_doc_line("      f(m)=|W(Lm - d_obs)|_2 where L is the convolution with the modeled data")
  call sep_add_doc_line("                             and m the filter coefficients. W is the weighting operator")
  call sep_add_doc_line("                             which is optional")
  call sep_add_doc_line("COMMENTS")
  call sep_add_doc_line("")
  call doc('SOURCE')
end subroutine Wavelet_estimation_doc
