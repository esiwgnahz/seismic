# 1 "<stdin>"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "<stdin>"
module hycdsolver_reg_mod
! 0 = W (F J m - d)
  use chain0_mod
  use hycdsolver_report_mod
  use hpf_mod ! 0 = A m
  implicit none
  logical, parameter, private :: T = .true., F = .false.
  contains
  subroutine hycdsolver_reg( m,d, Fop, Aop, Wop, Wmop, stepper, nAop,&
    & niter,eps , Jop,m0,rm0,err,resd,resm,mmov,rmov,verb)
    optional :: Jop,m0,rm0,err,resd,resm,mmov,rmov,verb
    interface
!-------------------------- begin definitions -----------
      integer function Fop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function Aop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function Wop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function Wmop(adj,add,m,d)
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
    real, dimension(:), intent(in) :: d, m0,rm0
    integer, intent(in) :: niter, nAop
    logical, intent(in) :: verb
    real, intent(in) :: eps
    real, dimension(:), intent(out) :: m,err, resd,resm
    real, dimension(:,:), intent(out) :: rmov,mmov
    real, dimension(size( m)) :: g, s, mtmp
    real, dimension(size( d) + nAop), target :: tmp, rs, gg, ss, tt,&
      & r
    real, dimension(:), pointer :: tmpd, gd, td, rd
    real, dimension(:), pointer :: tmpm, gm, tm, rm
    integer :: iter, stat, nf_eval=0&
      &,ng_eval=0,icount,ncount=10
    logical :: sd
    double precision :: alfa,alfa_k,beta,beta_k&
      &,of,ofp
    real :: f_new,f_old
    gd => gg(1:size(d))
    gm => gg(1+size(d):)
    td => tt(1:size(d))
    tm => tt(1+size(d):)
    rd => r(1:size(d))
    rm => r(1+size(d):)
    tmpd =>tmp(1:size(d))
    tmpm => tmp(1+size(d):)
    ss = 0. ! residual step
    s = 0. ! solution step
    alfa_k=0.
    beta_k=0.
    stat=Wop(F,F,-d,rd) ! rd = -W d
    rm = 0.
    if (present(rm0)) then
      rm=rm0 ! rm = R m0
    end if
    if (present( m0)) then
      m=m0 ! m = m0
      call chain0(Wop,Fop,F,T,m,rd,td)
      ! rd += W F m0
      call chain0(Wmop,Aop,F,T,m0,rm,tm)
      ! rm += W_m A m0
      nf_eval=nf_eval+1
    else
      m=0
    end if
    sd = T
!--------------------------- begin iterations -----------
    do iter = 1,niter
      tmpd=hpfp(rd)
      tmpm=eps*hpfp(rm) ! Rd, Rm
!Now computes \delta m
      call chain0(Wop,Fop,T,F,g,tmpd,td)
      ! g = (WF)'Rd
      call chain0(Wmop,Aop,T,T,g,tmpm,tm)
      ! g += A'Rm
      if (present(Jop)) then
        tm=g
        stat=Jop(F,F,tm,g )
      end if
! g = J g
      ng_eval=ng_eval+1
!Now computes \delta r
      call chain0(Wop,Fop,F,F,g,gd,td)
      ! Gd = (WF) g
      call chain0(Wmop,Aop,F,F,g,gm,tm)
      ! Gm = e W_m A g
      nf_eval=nf_eval+1
!Now computes coefficients
      stat = stepper(sd,alfa,beta,ss,rd,gg,rm,eps)
      ! m+=dm; R+=dR
      if (stat.eq.1) then
        exit ! got stuck descending
      end if
      mtmp = m + alfa*g + beta*s
      rd=-d
      stat=Fop(F,T,mtmp,rd) ! rd=Fm-d
      stat=Wop(F,F,rd,td)
      rd=td ! rd=W(Fm-d)
      call chain0(Wmop,Aop,F,F,mtmp,rm,tm)
      ! rm=W_m A m
      nf_eval=nf_eval+1
      of=hpf(rd)+eps*hpf(rm)
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
          call chain0(Wmop,Aop,F,F,mtmp,rm,tm)
          ! rm=W_m A m
          nf_eval=nf_eval+1
          icount=icount+1
          of=hpf(rd)+eps*hpf(rm)
        end do
      end if
      if (icount.ge.ncount) then
        write(0,*) 'INFO: Got stuck, keep last good m, exit'
        exit
      end if
      m = m + alfa*g + beta*s
      s = alfa*g + beta*s
      ss= alfa*gg+ beta*ss
      sd=F
      if (present( mmov)) then
        mmov(:,iter) = m(:size(mmov,1)) ! report -----
      end if
      if (present( rmov)) then
        rmov(:,iter) = r(:size(rmov,1))
      end if
      if (present( err )) then
        err( iter) = hpf(rd)+eps*hpf(rm)
      end if
      if (present( verb)) then
        if (verb) then
          call hycdsolver_report(iter,m,g,rd,rm,eps)
        end if
      end if
      ofp=of
    end do
    if (present(verb)) then
      if (verb) then
        write(0,*) 'INFO: -------------'
        write(0,*) 'INFO:  Diagnostic  '
        write(0,*) 'INFO: -------------'
        write(0,*) 'INFO:'
        write(0,*) 'INFO:  Total Number of forward operations ='&
          &,nf_eval
        write(0,*) 'INFO:  Total Number of adjoint operations ='&
          &,ng_eval
        write(0,*) 'INFO:'
        write(0,*) 'INFO: -------------'
      end if
    end if
    if (present( resd)) then
      resd = rd
    end if
    if (present( resm)) then
      resm = rm(:size(resm))
    end if
  end subroutine
end module
