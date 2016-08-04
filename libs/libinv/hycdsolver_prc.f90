# 1 "<stdin>"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "<stdin>"
module hycdsolver_prc_mod
! 0 = W (F S J p - d)
  use chain0_mod
  use hycdsolver_report_mod
  use hpf_mod ! 0 = I p
  implicit none
  logical, parameter, private :: T = .true., F = .false.
  contains
  subroutine hycdsolver_prc( m,d, Fop, Sop, Wop, Wmop, stepper, nSop,&
    & niter,eps , Jop,p0,rm0,err,resd,resm,mmov,rmov,verb)
    optional :: Jop,p0,rm0,err,resd,resm,mmov,rmov,verb
    interface
!-------------------------- begin definitions -----------
      integer function Fop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function
      integer function Sop(adj,add,m,d)
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
    real, dimension(:), intent(in) :: d, p0,rm0
    integer, intent(in) :: niter, nSop
    logical, intent(in) :: verb
    real, intent(in) :: eps
    real, dimension(:), intent(out) :: m,err, resd,resm
    real, dimension(:,:), intent(out) :: rmov,mmov
    real, dimension(size( m)) :: p , g, s, ptmp
    real, dimension(size( d) + nSop), target :: rr, gg, ss, tt, tmp
    real, dimension(:), pointer :: rd, gd, td, tmpd
    real, dimension(:), pointer :: rm, gm, tm, tmpm
    integer :: iter, stat,nf_eval=0&
      &,ng_eval=0,icount,ncount=10
    logical :: sd
    double precision :: alfa,beta,of,ofp
    rd => rr(1:size(d))
    rm => rr(1+size(d):)
    gd => gg(1:size(d))
    gm => gg(1+size(d):)
    td => tt(1:size(d))
    tm => tt(1+size(d):)
    tmpd => tmp(1:size(d))
    tmpm => tmp(1+size(d):)
    ss = 0. ! residual step
    s = 0. ! solution step
    stat=Wop(F,F,-d,rd)
    ! rd = -W d
    rm = 0.
    if (present(rm0)) then
      rm=rm0 ! rm = Rm0
    end if
    if (present( p0)) then
      p=p0 ! p = p0
      call chain0(Wop,Fop,Sop,F,T,p,rd,tm,td)
      ! rd += WFS p0
      stat=Wmop(F,T,p,rm)
      ! rm += e I p0
      nf_eval=nf_eval+1
    else
      p=0
    end if
    sd = T
!-------------------------- begin iterations ------------
    do iter = 1,niter
      tmpd=hpfp(rd)
      tmpm=eps*hpfp(rm)
      ! First derivative hyperbolic norm
!Now computes \delta p
      call chain0(Wop,Fop,Sop,T,F,g,tmpd,tm,td)
      ! g = (WFS)'Rd
      stat=Wmop(F,T,tmpm,g)
      ! g += e I'Rm
      if (present(Jop)) then
        tm=g
        stat=Jop(F,F,tm,g)
      end if
! g = J g
      ng_eval=ng_eval+1
!Now computes \delta r
      call chain0(Wop,Fop,Sop,F,F,g,gd,tm,td)
      ! Gd = (WFS) g
      stat=Wmop(F,F,g,gm)
      ! Gm = e W_mI g
      nf_eval=nf_eval+1
!Now computes coefficients
      stat = stepper(sd,alfa,beta,ss,rd,gg,rm,eps)
      ! m+=dm; R+=dR
      if (stat .eq.1) then
        exit
        ! got stuck descending
      end if
      ptmp= p + alfa*g + beta*s
      rd=-d
      call chain0(Fop,Sop,F,T,ptmp,rd,tm)
      ! rd=FSp-d
      stat=Wop(F,F,rd,td)
      rd=td ! rd=W(FSp-d)
      stat=Wmop(F,F,ptmp,rm)
      ! rm=eW_mIp
      nf_eval=nf_eval+1
      of=hpf(rd)+eps*hpf(rm)
! If we are not decreasing the value of the OF, we divide the steps by 2
      if (iter.gt.1) then
        icount=0
        do while ((of.gt.ofp).and.(icount.le.ncount))
          alfa=alfa/2
          beta=beta/2
          ptmp= p + alfa*g + beta*s
          rd=-d
          call chain0(Fop,Sop,F,T,ptmp,rd,tm)
          ! rd=FSp-d
          stat=Wop(F,F,rd,td)
          rd=td ! rd=W(FSp-d)
          stat=Wmop(F,F,ptmp,rm)
          ! rm=eW_mIp
          nf_eval=nf_eval+1
          icount=icount+1
          of=hpf(rd)+eps*hpf(rm)
        end do
      end if
      if (icount.ge.ncount) then
        write(0,*) 'INFO: Got stuck, keep last good m, exit'
        exit
      end if
      p = p + alfa*g + beta*s
      s = alfa*g + beta*s
      ss= alfa*gg+ beta*ss
      sd= F
      stat = Sop(F,F,p,m)
      ! m = S p
      if (present( mmov)) then
        mmov(:,iter) = m(:size(mmov,1)) ! report -----
      end if
      if (present( rmov)) then
        rmov(:,iter) = rr(:size(rmov,1))
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
    if (present( resm)) then
      resm = rm(:size(resm))
    end if
  end subroutine
end module
