# 1 "<stdin>"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "<stdin>"
module hycdsolver_smp_mod
! 0 = W (F J m - d)
  use chain0_mod
  use hycdsolver_report_mod
  use hpf_mod
  implicit none
  logical, parameter, private :: T = .true., F = .false.
  contains
  subroutine hycdsolver_smp( m,d, Fop, Wop, stepper, niter , &
    & Jop,m0,err,resd,mmov,rmov,verb)
    optional :: Jop,m0,err,resd,mmov,rmov,verb
    interface
!-------------------------- begin definitions -----------
      integer function Fop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function Wop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function Jop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function stepper(sd,alfa,beta,ss,rd,gg,rm,eps)
        real, dimension(:) :: ss,rd,gg
        real, dimension (:), optional :: rm
        double precision :: alfa,beta
        logical :: sd
        real,optional :: eps
      end function
    end interface
    real, dimension(:), intent(in) :: d, m0
    integer, intent(in) :: niter
    logical, intent(in) :: verb
    real, dimension(:), intent(out) :: m,err, resd
    real, dimension(:,:), intent(out) :: rmov, mmov
    real, dimension(size(m)) :: g, s, mtmp
    real, dimension(size(d)), target :: rr, gg, ss
    real, dimension(size(d)+size(m)), target :: tt
    real, dimension(:), pointer :: rd, gd, td
    real, dimension(:), allocatable :: tmpd
    real, dimension(:), pointer :: rm, gm, tm
    integer :: iter, stat,nf_eval=0&
      &,ng_eval=0,icount,ncount=10
    logical :: sd
    double precision :: alfa, beta,of,ofp
    allocate(tmpd(size(d)))
    rd => rr(1:size(d))
    gd => gg(1:size(d))
    td => tt(1:size(d))
    tm => tt(1+size(d):)
    ss = 0. ! residual step
    s = 0. ! solution step
    stat=Wop(F,F,-d,rd)
    ! Rd = -W d
    if (present( m0)) then
      m=m0 ! m = m0
      call chain0(Wop,Fop,F,T,m,rd,td)
      ! Rd+= WF m0
      nf_eval=nf_eval+1
    else
      m=0
    end if
    sd = T
!-------------------------- begin iterations ------------
    do iter = 1,niter
      tmpd=hpfp(rd)
      ! First derivative hyperbolic norm
!Now computes \delta m
      call chain0(Wop,Fop,T,F,g,tmpd,td)
      ! g = (WF)'Rd
      if (present(Jop)) then
        tm=g
        stat = Jop(F,F,tm, g )
      end if
! g = J g
      ng_eval=ng_eval+1
!Now computes \delta r
      call chain0(Wop,Fop,F,F,g,gd,td)
      ! Gd = (WF) g
      nf_eval=nf_eval+1
!Now computes coefficients
      stat = stepper(sd,alfa,beta,ss,rr,gg)
      ! Computes alfa and beta for conj. directions
      if (stat.eq.1) then
        exit ! got stuck descending
      end if
      mtmp = m + alfa*g + beta*s
      rd=-d
      stat=Fop(F,T,mtmp,rd) ! rd=Fm-d
      stat=Wop(F,F,rd,td)
      rd=td ! rd=W(Fm-d)
      nf_eval=nf_eval+1
      of=hpf(rd)
! If we are not decreasing the value of the OF, we divide the steps by 2
      if (iter.gt.1) then
        icount=0
        do while ((of.gt.ofp).and.(icount.le.ncount))
          alfa=alfa/2
          beta=beta/2
          mtmp= m + alfa*g + beta*s
          rd=-d
          stat=Fop(F,T,mtmp,rd)
          ! rd=Fm-d
          stat=Wop(F,F,rd,td)
          rd=td ! rd=W(Fm-d)
          nf_eval=nf_eval+1
          icount=icount+1
          of=hpf(rd)
        end do
      end if
      if (icount.ge.ncount) then
        write(0,*) 'INFO: Got stuck, keep last good m, exit'
        exit
      end if
      m = m + alfa*g + beta*s
      s = alfa*g + beta*s
      ss= alfa*gg+ beta*ss
      sd= F
      if (present( mmov)) then
        mmov(:,iter) = m(:size(mmov,1)) ! report -----
      end if
      if (present( rmov)) then
        rmov(:,iter) = rd(:size(rmov,1))
      end if
      if (present( err )) then
        err( iter) = hpf(rd)
      end if
      if (present( verb)) then
        if (verb) then
          call hycdsolver_report(iter,m,g,rd)
        end if
      end if
      ofp=of
    end do
    write(0,*) 'INFO: -------------'
    write(0,*) 'INFO:  Diagnostic  '
    write(0,*) 'INFO: -------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  Total Number of forward operations =',nf_eval
    write(0,*) 'INFO:  Total Number of adjoint operations =',ng_eval
    write(0,*) 'INFO:'
    write(0,*) 'INFO: -------------'
    if (present( resd)) then
      resd = rd
    end if
    if (allocated(tmpd)) then
      deallocate(tmpd)
    end if
  end subroutine
end module
