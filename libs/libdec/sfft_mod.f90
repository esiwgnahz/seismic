module sfft_mod
  use fftw_mod 
  use omp_lib

  implicit none

  integer, private :: nx,nw,ny,nz

  real :: pi=3.141592653589793238462643

contains

  ! Initializing dimensions
  ! FFTW for the 2nd and 3rd axes goes from 0 to pi -pi to 0
  ! so the last n/2+1 to n samples are the negative frequencies from
  ! -pi to 0.

  subroutine init_fft_dimension_1d_r2c_c2r(nx_in,dx_in,dw,ow)
    integer :: nx_in
    real    :: dx_in

    real    :: dw
    real    :: ow

    nx=nx_in
    nw=nx_in/2+1
    dw=2*pi/(nx_in*dx_in)
    ow=0.

  end subroutine init_fft_dimension_1d_r2c_c2r

  subroutine init_fft_dimension_2d(ny_in,nz_in,dy_in,dz_in,dky,dkz,oky,okz)
    integer :: ny_in,nz_in
    real    :: dy_in,dz_in

    real    :: dky,dkz
    real    :: oky,okz
  
    nz=nz_in
    ny=ny_in
    dkz=2*pi/(nz_in*dz_in)
    dky=2*pi/(ny_in*dy_in)
    okz=0.
    oky=0.
  
  end subroutine init_fft_dimension_2d

  ! Initializing plans
  subroutine initialize_fft1d_double_r2c_c2r(plan)
    complex        :: freq(nw)
    real           :: data(nx)
    integer(8), dimension(2) :: plan

    call sfftw_plan_dft_r2c_1d(plan(1),nx,data,freq,FFTW_ESTIMATE)
    call sfftw_plan_dft_c2r_1d(plan(2),nx,freq,data,FFTW_ESTIMATE)
  end subroutine initialize_fft1d_double_r2c_c2r

  ! Initializing plans
  subroutine initialize_fft1d_r2c_c2r(plan)
    complex        :: freq(nw)
    real           :: data(nx)
    integer(8), dimension(2) :: plan

    call sfftw_plan_dft_r2c_1d(plan(1),nx,data,freq,FFTW_ESTIMATE)
    call sfftw_plan_dft_c2r_1d(plan(2),nx,freq,data,FFTW_ESTIMATE)
  end subroutine initialize_fft1d_r2c_c2r

  subroutine fft1d_r2c_c2r_n23(dat,wld,n23,plan,back)    ! scaled 1D FFT T <-> W
    complex        :: wld(nw,n23)
    real           :: dat(nx,n23)
    logical, intent(IN)      :: back
    integer(8), dimension(2) :: plan
    integer                  :: i,n23

    if(back) then
       do i=1,n23
          call sfftw_execute_dft_c2r(plan(2),wld(:,i)/nx,dat(:,i))
       end do
    else
       do i=1,n23
          call sfftw_execute_dft_r2c(plan(1),dat(:,i),wld(:,i))
       end do
    endif
  end subroutine fft1d_r2c_c2r_n23

  subroutine fft1d_double_r2c_c2r_n23(dat,wld,n23,plan,back)    ! scaled 1D FFT T <-> W
    complex       :: wld(nw,n23)
    real          :: dat(nx,n23)
    logical, intent(IN)      :: back
    integer(8), dimension(2) :: plan
    integer                  :: i,n23

    if(back) then
       do i=1,n23
          call sfftw_execute_dft_c2r(plan(2),wld(:,i)/nx,dat(:,i))
       end do
    else
       do i=1,n23
          call sfftw_execute_dft_r2c(plan(1),dat(:,i),wld(:,i))
       end do
    endif
  end subroutine fft1d_double_r2c_c2r_n23


  subroutine destroy_fft(plan)
    integer(8), dimension(2) :: plan

    call sfftw_destroy_plan(plan(1))
    call sfftw_destroy_plan(plan(2))
  end subroutine destroy_fft

end module sfft_mod
