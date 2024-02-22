! (c) F. D. M. Haldane, Princeton University, haldane@princeton.edu
! 2014-07-31
!-----------------------------------------------------------------
subroutine  lanczos_iterator(n,vector_is_real,nevec,evec,eval,var,err, &
     matvec,seed,maxstep,info,work,mask_val,mask_vec)
  use lanczos_m     !module supplying  Ritz vector procedures
  use tmatrix_m     !module supplying  T-matrix procedures
  use vector_m      !module supplying  Hilbert space basic linear algebra 
  implicit none 
  integer,parameter :: dp = kind(1.0d0)
  integer, intent(in)  :: n                                       ! Hilbert-space dimension
  logical, intent(in) :: vector_is_real                           ! selects real or complex hermitian
  integer, intent(in) ::  nevec                                   ! number of previously-found eigenvectors
  type(vector_t), intent(inout), dimension(nevec+1) :: evec       ! eigenvectors (type defined in module vector_m)
  real (kind=dp), intent(inout),dimension(nevec+1) :: eval, var   ! eigenvalues and eigenvector variances
  real (kind=dp), intent(inout) :: err                            ! resolution limit, calculated only when nevec = 0
  external matvec                                                 ! name of the  matrix-vector multiplication routine
  integer, intent(inout) :: seed(4)                               ! seed for randam number generator
  integer, intent(in)  :: maxstep                                 ! maximum allowed number of lanczos steps
  integer, intent(inout) :: info(2)                               ! controls information output about process
  type (vector_t), target, intent(inout) :: work(3)               ! three work vectors
  integer, intent(in), optional :: mask_val,mask_vec(n)           ! mask for initial state of lanczos process

  !--------------------------------------------
  ! find the (NEVEC+1)''th  eigenvalue EVAL(NEVEC) and corresponding 
  ! eigenvector TYPE (VECTOR_T) EVEC(NEVEC) of an NxN real or complex  Hermitian matrix
  ! specified implicitly by an external  matrix-vector multiplication  procedure MATVEC
  !
  ! Assumes that NEVEC previously-found eigenvalues and eigenvector sare available
  !
  !  EVEC(I)%N  = dimension of eigenvector (integer)
  !  EVEC(I)%IS_REAL  identifies it as real or complex (logical)
  !  EVEC(I)%R  is a double precision real pointer to the eigenvector if it is real
  !  EVEC(I)%C is a double precision complex  pointer to the eigenvector if it is complex
  !-----------------------------------------------
  !The external user-supplied matrix-vector multiplication
  !procedure MATVEC  must perform an action equivalent to:
  !   subroutine matvec(n,vec_in,vec_out,add)
  !    integer n
  !    double precision vec_in(n), vec_out(n)
  !    logical add
  !    double precision  matrix(n,n) !the matrix  
  !    if(.not.add) vec_out = 0d0
  !    vec_out  = vec_out + matmul(matrix,vec_in)
  !    endif
  !    return
  !   end subroutine matvec
  !-------------------------------------------------
  real (kind=dp), parameter :: zero = 0.d0, one = 1.d0
  real (kind=dp), save  ::property(3)
  type (vector_t), pointer :: vec1_p => null(), vec2_p => null()
  real (kind=dp) :: d,e,shift,ev,norm, var_prev,eb,e2,reduction,muc,vb,prev_reduction,gap
  integer :: run, nstep,i,report=1
  logical :: done
  real (kind=dp), allocatable, dimension(:,:):: cmat
  real (kind=dp), allocatable, dimension(:):: w
  integer :: lwork
  logical :: eval_converged = .false., no_var_prev = .true.
  logical:: debug = .false., give_info = .false.
  complex (kind=dp) :: sp
  
  real (kind=dp), allocatable, target, dimension(:)::  w1, w2
  real (kind=dp), allocatable, target, dimension(:,:):: tmat_1, tmat_2
  real (kind=dp), allocatable, target, dimension(:):: e0_1, e0_2

  real (kind=dp), pointer, dimension(:) :: w_p => null()
  real (kind=dp), pointer, dimension(:,:) :: tmat => null()
  real (kind=dp), pointer, dimension(:) :: e0 => null()

  integer, parameter :: initial_storage = 1024
  integer:: storage,  max_nstep, nstep0
  
  if(info(1) < 0) then
     debug = .true.
     report = -info(1)
     give_info = .true.
  endif
  if(nevec==0 .or. err == 0_dp) then
     call matrix_properties(n,vector_is_real,matvec,seed,property,give_info)   ! from lanczos_m
     err = property(3)    !this is the  maximum possible eigenvalue precision
     gap = 0_dp
  endif
  call vector_random_fill(evec(nevec+1),seed) ! place a random starting  vector in evec(nevec+1)

  ! optional projection with mask, if present
  if (present(mask_val)) then
     if( mask_val > 0 .and. mask_val <= n) then
        if (mask_val /= 0) then ! in this case set element i  with mask_vec(i) /= mask_val to zero
           if (count(mask_vec == mask_val) == 0) then
              write(6,'("LANCZOS_ITERATOR called with empty MASK_VEC = MASK_VAL =",i10)') mask_val
              stop
           endif
           call vector_mask(evec(nevec+1),n,mask_vec,mask_val) ! set element i to zero if mask_vec(i) = 0
        else   ! case where mask_val = 0
           if (count(mask_vec == 0) == n) then
              write(6,'("LANCZOS_ITERATOR called with empty MASK_VEC, MASK_VAL=:",i10)') mask_val 
              stop
           endif
           call vector_mask(evec(nevec+1),n,mask_vec)  
        endif
     endif
  endif
  eval_converged = .false.
  run = 0

  call new_tmat_storage(initial_storage)
  max_nstep = storage
  do while (.not.eval_converged)
     nstep = 0
     nstep0 = 0
     run =  run + 1
     call vector_normalize(evec(nevec+1))
     
     if (nevec > 0) then
        do i = 1, nevec   ! orthogonalize to previously-found eigenvectors
           call vector_orthogonalize(evec(i),evec(nevec+1))
        enddo
        call vector_normalize(evec(nevec+1))
     endif
     
     call vector_copy(evec(nevec+1),work(1))    ! initialize the lanczos iteration vectors 
     call vector_initialize(work(2),0)
     
     if (nevec > 0) then ! get an upper bound to gap to next eigenvalue
        call vector_matvec(work(1),work(2),matvec,.false.)
        call vector_scalar_product(work(1), work(2), sp) 
        gap = real(sp) - eval(nevec)
        call vector_initialize(work(2),0)
     endif
     done = .false.
     
     no_var_prev = .true.
     shift = zero
     var(nevec+1) = 0_dp
     
     vec1_p => work(1); vec2_p => work(2) !pointers to the two work vectors
     
     do while(.not.done)
        if (nstep + 1 > max_nstep) exit
        call lanczos_iterate(nstep,d,e,vec1_p,vec2_p,matvec,shift,nevec,evec,eval,gap)

        if (.not. associated(tmat)) then
           write(6,'("tmat error,storage =",3i10)') storage, nstep, max_nstep
           stop
        endif
  
        tmat(1,nstep) = d + shift
        tmat(2,nstep) = e
        if (nstep == storage)  then
           storage = 2*storage
           call new_tmat_storage(storage)
           max_nstep = storage
        endif
           
        if(e.lt.err) done = .true.
        !     carry out at least two steps unless  tmat(2,1)  < matverr.
        !     after  nstep=2 check for termination without eigenvector improvement.
        !     store the lowest eigenvalue of T(nstep) in e0(nstep)
        if(nstep.eq.1) then
           e0(1) = tmat(1,1)
           shift =  - e0(1)
           if(debug.and.mod(nstep,report)==0) write(6,'("run,nstep: ",2i5,f25.16)') &
                run,nstep,e0(nstep)
           cycle
        endif
        call tmatrix_eval(nstep,tmat,e0)
        
        shift = -e0(nstep)
        if(.not.eval_converged.and.e0(nstep).gt.e0(nstep-1)-err) eval_converged = .true.
        ev = e0(nstep)
        if(.not.eval_converged) then
           if(debug.and.mod(nstep,report)==0) write(6,'("run,nstep: ",2i5,f25.16)') &
                run,nstep,e0(nstep)
           if (maxstep == 0) cycle
           if (nstep < maxstep) cycle
           if (allocated(w1)) then
              if (size(w1) < nstep) deallocate(w1)
           endif
           if (.not.allocated(w1)) allocate (w1(nstep))
           call tmatrix_evec(nstep,tmat,ev,w1,var(nevec+1))
           w_p => w1
           exit
        endif
        if(nstep0 == 0) then
           nstep0 = nstep
           if (give_info) print '("eigenvalue has converged at nstep =",i6,f25.16 )',&
                nstep0, e0(nstep0)
           max_nstep = 2*nstep0
           if (allocated(w1)) then
              if (size(w1) < max_nstep) deallocate(w1)
           endif
           if (allocated(w2)) then
              if (size(w2) < max_nstep) deallocate(w2)
           endif
           if (.not. allocated(w1)) allocate (w1(max_nstep))
           if (.not. allocated(w2)) allocate (w2(max_nstep))
        endif
        if (modulo(nstep-nstep0,2) == 0) then
           w_p => w1
        else
           w_p => w2
        endif
        var_prev = var(nevec+1)
        call tmatrix_evec(nstep,tmat,ev,w_p,var(nevec+1))
        if(debug.and.mod(nstep,report)==0) &
             write(6,'("run,nstep: ",2i5,f25.16,1es12.3)') run,nstep,ev,var(nevec+1)
        if (nstep == nstep0) cycle
        if(var(nevec+1).gt.var_prev-err) then
           done = .true.
           nstep = nstep -1
           var(nevec+1) = var_prev
           if (modulo(nstep-nstep0,2) == 0) then
              w_p => w1
           else
              w_p => w2
           endif
           exit
        endif
     end do

     ! w(1:nstep) is the t-matrix eigenvector, e0(nstep) and var(nevec+1) are the predicted 
     ! eigenvalues and variance.
     ! now construct the Ritz vector
     
     call vector_copy(evec(nevec+1),vec1_p)
     call vector_initialize(vec2_p,0)
     call vector_scale_by_real(w_p(1),evec(nevec+1))
     
     i = 0
     shift = zero
     do while (i+1.lt.nstep)
        call lanczos_iterate(i,d,e,vec1_p,vec2_p,matvec,shift,nevec,evec,eval,gap)
        !     now i = 1,2,3...nstep-1
        !     vec1_p points to |i>  = |1>,|2>,...,|nstep-1>
        !     vec2_p points to |i+1>  = |2>,|3>,...,|nstep>  
        call vector_axpy(cmplx(w_p(i+1),kind=dp),vec2_p,evec(nevec+1))
        shift =  - e0(i)
     enddo
     ! if .not.eval_converged, we should repeat the lanczos process using
     ! evec(nevec+1) as the new starting vector
  enddo

  nullify(w_p)
  if (allocated(w1)) deallocate(w1)
  if (allocated(w2)) deallocate(w2)
  
  call residual(evec(nevec+1),work(1),matvec,eval(nevec+1),var(nevec+1))
  if(give_info) write(6,'(" basic lanczos :, nstep=",i12," eigenvalue =", &
       & f25.16," variance =",1es12.3)') nstep, eval(nevec+1),var(nevec+1)
  ! polish evec
  var_prev = var(nevec+1)
  info(1) = nstep

  nstep = 0
  polish:   if(var_prev.gt.err) then
   lwork = 2*max_nstep
   allocate (cmat(2,max_nstep),w(lwork))
   call vector_scale_by_real(1/var_prev,vec1_p)
   call vector_copy(vec1_p,work(3))
   done = .false.
   eb = zero
   shift = eval(nevec+1)
   do while (.not.done)
      if(nstep+1.ge.max_nstep) exit
      call lanczos_iterate(nstep,d,e,vec1_p,vec2_p,matvec,shift,nevec+1,evec,eval,gap)
      cmat(1,nstep) = d - eb                 
      if(e.lt.err) then
         done = .true.
         e = zero
      else if (cmat(1,nstep).lt.err) then
         done = .true.
      else 
         cmat(2,nstep) = e/cmat(1,nstep)
         eb = e*cmat(2,nstep)
      endif
      e2 = e**2
      call lanczos_inverse_iteration(nstep,cmat,e2,reduction,vb,muc,w,lwork)
      if(nstep.gt.1.and.(reduction/prev_reduction).gt.0.995d0) done = .true.
      if(reduction*var_prev.lt.err/2) done = .true.
      if(debug.and.mod(nstep,report)==0) write(6,'("polish nstep:",i5,& 
           & " variance  reduction",1es18.9,l5)') nstep,reduction,done
      if(reduction.gt.one) exit
      prev_reduction = reduction
   enddo
   if(done) then
      call lanczos_invit_getvec(nstep,vb,muc,w,lwork)
      ! the improved t-matrix eigenvector is returned in w(1:nstep); rerun lanczos
      call vector_copy(work(3),vec1_p)
      call vector_scale_by_real(w(1),work(3))
      i = 0
      do while (i+1.lt.nstep)
         call lanczos_iterate(i,d,e,vec1_p,vec2_p,matvec,shift,nevec+1,evec,eval,gap)
         call vector_axpy(cmplx(w(i+1),kind=dp),vec2_p,work(3))
      enddo

      call vector_normalize(work(3))
      call lanczos_maxpy(evec(nevec+1),work(1),.false.,matvec,shift,nevec+1,evec,eval,zero)
      call lanczos_maxpy(work(3),work(2),.false.,matvec,shift,nevec+1,evec,eval,zero)

      call polish_subtraction(evec(nevec+1),work(1),work(3),work(2))
      call residual(evec(nevec+1),work(1),matvec,eval(nevec+1),var(nevec+1))

   endif
   deallocate(cmat,w); 
end if polish

info(2) = nstep 

if(give_info) then
   write(6,'(" polish lanczos:, nstep=",i12," variance =",1es12.3)') nstep, var(nevec+1)
endif

return
contains
  subroutine new_tmat_storage(new)
    integer, intent(in) :: new
    logical :: alloc1, alloc2, error = .false.
    integer :: old

    nullify(tmat,e0)
    alloc1 = allocated(tmat_1)
    alloc2 = allocated(tmat_2)
    
    if (alloc1 .and. .not.alloc2) then
       old = size(e0_1)
       if (new == 0) then
          deallocate(tmat_1,e0_1)
       else if (new > old) then
          allocate(tmat_2(2,new), e0_2(new))
          tmat_2(:,1:old) = tmat_1(:,1:old); e0_2(1:old) = e0_1(1:old)
          tmat => tmat_2; e0 => e0_2
          deallocate(tmat_1, e0_1)
       else
          error = .true.
       end if
    else if (.not.alloc1 .and. alloc2) then
       old = size(e0_2)
       if (new == 0) then
          deallocate(tmat_2,e0_2)
       else if (new > old) then
          allocate(tmat_1(2,new), e0_1(new))
          tmat_1(:,1:old) = tmat_2(:,1:old); e0_1(1:old) = e0_2(1:old)
          tmat => tmat_1; e0 => e0_1
          deallocate(tmat_2, e0_2)
       else
          error = .true.
       end if
    else if (.not.alloc1 .and. .not.alloc2) then
       old = 0
       allocate(tmat_1(2,new), e0_1(new))
       tmat => tmat_1; e0 => e0_1
    else if (alloc1 .and. alloc2) then
       write(6,'("lanczos_interator storage error")')
       write(6,'("both tmat_1 and tmat_2 are allocated")')
       stop
    endif
    if (.not.error) then
       storage = new
       return
    endif
    write(6,'("lanczos_iterator storage error", 2i10)') old, new
    stop
  end subroutine new_tmat_storage
end subroutine lanczos_iterator

subroutine get_residual(evec, residual_vec, matvec, eigenvalue, variance, norm)
  use lanczos_m
  use vector_m
  implicit none
  external matvec
  integer, parameter :: dp = kind(1.0d0)
  type (vector_t), intent(inout) :: evec, residual_vec
  real (dp), intent(out) :: eigenvalue, variance, norm
  call residual(evec, residual_vec, matvec, eigenvalue, variance, norm)
  return
end subroutine get_residual
