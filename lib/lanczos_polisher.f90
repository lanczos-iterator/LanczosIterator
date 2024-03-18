subroutine lanczos_polisher(work,matvec,eval,evec,var,max_nstep,nstep_total,err,debug,report)
  use lanczos_m     !module supplying  Ritz vector procedures
  use tmatrix_m     !module supplying  T-matrix procedures
  use vector_m      !module supplying  Hilbert space basic linear algebra 
  implicit none
  integer,parameter :: dp = kind(1.0d0)
  external :: matvec
  integer, intent(in) :: max_nstep, report
  logical, intent(in) :: debug
  real (dp), intent(in) :: err
  real(dp),intent(inout) :: eval, var
  type (vector_t), intent(inout) :: evec
  type (vector_t), intent(inout), target :: work(3)
  integer, intent(out) :: nstep_total
  !-----------------------------------------------
  ! "polishes" an approximate eigenvector with a converged
  ! eigenvalue.   Does not reference previously-converged
  ! eigenvectors or eigenvalues.
  !----------------------------------------------
  integer :: lwork, i, nstep
  real(dp), allocatable :: cmat(:,:), w(:)
  real(dp) :: shift, eb, d, e, e2, reduction, prev_reduction, var_prev, muc, gap, vb
  real(dp), parameter:: zero = 0_dp, one = 1_dp
  type (vector_t), pointer::  vec1_p, vec2_p
  logical :: done
  type (vector_t) :: evec_not_used(1)
  real (dp) :: eval_not_used(1)

  gap = 0_dp
  nstep_total = 0
  ! polish evec
  print '("POLISH")'
  shift = eval
  if(var.gt.err/2) then
   lwork = 2*max_nstep
   allocate (cmat(2,max_nstep),w(lwork))
   nstep = 0
   var_prev = var
   call vector_scale_by_real(1/var_prev,work(1))
   call vector_copy(work(1),work(3))
   vec1_p => work(1)
   vec2_p => work(2)
   done = .false.
   eb = zero
   prev_reduction = zero
   do while (.not.done)
      if(nstep+1.ge.max_nstep) exit
      call lanczos_iterate(nstep,d,e,vec1_p,vec2_p,matvec,shift,0,&
           evec_not_used,eval_not_used,gap)
      cmat(1,nstep) = d - eb                 
      if(e.lt.err) then
         done = .true.
         e = zero
      else if (abs(cmat(1,nstep)).lt.err) then
         done = .true.
      else 
         cmat(2,nstep) = e/cmat(1,nstep)
         eb = e*cmat(2,nstep)
      endif
      e2 = e**2
      call lanczos_inverse_iteration(nstep,cmat,e2,reduction,vb,muc,w,lwork)
      if(debug.and.mod(nstep,report)==0) write(6,'("polish nstep:",i5,& 
           & " variance  reduction",1es18.9,l5)') nstep,reduction,done
      if(reduction*var_prev.lt.err/2) exit
      if (prev_reduction == zero) cycle
      if (prev_reduction < reduction) exit
      prev_reduction = reduction
   enddo
   call lanczos_invit_getvec(nstep,vb,muc,w,lwork)
   ! the improved t-matrix eigenvector is returned in w(1:nstep); rerun lanczos
   call vector_copy(work(3),vec1_p)
   call vector_scale_by_real(w(1),work(3))
   i = 0
   do while (i+1.lt.nstep)
      call lanczos_iterate(i,d,e,vec1_p,vec2_p,matvec,shift,&
           0,evec_not_used,eval_not_used,gap)
      call vector_axpy(cmplx(w(i+1),kind=dp),vec2_p,work(3))
   enddo
   
   call vector_normalize(work(3))
   call lanczos_maxpy(evec,work(1),.false.,matvec,shift,&
        0,evec_not_used,eval_not_used,zero)
   call lanczos_maxpy(work(3),work(2),.false.,matvec,shift,&
        0,evec_not_used,eval_not_used,zero)
   
   call polish_subtraction(evec,work(1),work(3),work(2))
   call residual(evec,work(1),matvec,eval,var)
   print '(" lanczos_polisher: variance ",1es10.3)', var
   nstep_total = nstep_total + nstep
   deallocate(cmat,w); 
end if
return
end subroutine lanczos_polisher
