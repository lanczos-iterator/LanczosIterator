subroutine rs_lanczos(n, nevec, eval, evec, variance, resolution,&
     matvec, maxstep, report, seed, range, nocheck)
  use vector_m
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in) :: n, nevec, maxstep, report, seed
  real(dp), intent(out) :: eval(nevec), variance(nevec)
  real(dp), intent(out), target :: evec(n,nevec)
  real(dp), intent(out)  :: resolution
  integer, intent(inout), optional :: range(2)
  logical, intent(in), optional :: nocheck
  external  matvec
!---------------------------------------------
!  finds lowest NEVEC eigenvalues and eigenvectors of a NxN real symmetric
!  matrix, supplied implicitly as an external subroutine MATVEC
!   N                  dimension of matrix (input)
!   NEVEC              number of eigenvectors requested (input)
!   EVAL(NEVEC)        eigenvalues (output)
!   EVEC(N,NEVEC)      eigenvectors (output)
!   VARIANCE(NEVEC)    variance of eigenvectors (output)
!   RESOLUTION         limit on eigenvalue resolution and variance (output)
!   MATVEC             name of matrix-vector multiplication procedure (input)
!   MAXSTEP            maximum number of lanczos steps allowed (input)
!   REPORT             if > 0, show lanczos progress after REPORT steps (input)
!   SEED               if > 0 on input, used as a random number seed (input).
!   RANGE(2)(optional) only compute eigenvectors [max(1,RANGE(1)):min(NEVEC,RANGE(2))]
!                      (if RANGE(1) > 1 data for eigenvalues [1:RANGE(1)-1] is required)
!                      on output: contains the  range of eigenvectors found
!   NOCHECK (optional) If RANGE is present and NOCHECK = .true., omit check of
!                      eigenvectors [1:RANGE(1)-1] (don't omit if .false.) 
!-------------------------------------------------
!
!   MATVEC(N,VEC_IN,VEC_OUT,ADD):
!
!    on  input:
!     vec_in(i)  = A(i), i = 1:n
!     vec_out(i) = B(i)  i = 1:n
!    on output:
!      vec_in, n, and add are unchanged
!      vec_out(i) = C(i) i = 1:n
!
!   (a) if add = .true.  
!     C(i)  = B(i) +  sum_j M(i,j)*A(j)
!
!   (b) if add = .false.
!     C(i)  =  sum_j M(i,j)*A(j)
!
!---------------------------------------------

    interface
       subroutine  lanczos_iterator(n,vector_is_real,nevec,evec,eval,var,err,&
            matvec,seed,maxstep,info,work,mask_val,mask_vec)
         use lanczos_m     !module supplying  Ritz vector procedures
         use tmatrix_m     !module supplying  T-matrix procedures
         use vector_m      !module supplying  Hilbert space basic linear algebra 
         implicit none 
         integer,parameter :: dp = kind(1d0)
         integer, intent(in)  :: n                                       ! Hilbert-space dimension
         logical, intent(in) :: vector_is_real                           ! selects real or complex hermitian
         integer, intent(in) ::  nevec                                   ! number of previously-found eigenvectors
         type (vector_t), intent(inout), dimension(nevec+1) :: evec      ! eigenvectors (type defined in module vector_m)
         real (kind=dp), intent(inout), dimension(nevec+1) :: eval, var  ! eigenvalues and eigenvector variances
         real (kind=dp), intent(inout) :: err                            ! resolution limit, calculated only when nevec = 0
         external matvec                                                 ! name of the  matrix-vector multiplication routine
         integer, intent(inout) :: seed(4)                               ! seed for randam number generator, if > 0
         integer, intent(in)  :: maxstep                                 ! maximum allowed number of lanczos steps
         integer, intent(inout) :: info(3)                               ! controls information output about process
         type (vector_t), intent(inout) :: work(3)                       ! three work vectors
         integer, intent(in), optional :: mask_val,mask_vec(n)           ! mask for initial state of lanczos process
       end subroutine lanczos_iterator
    end interface

    integer :: k, kmin, kmax, info(3), seed0
    integer, save , dimension(4) :: iseed  = (/222,333,444,555/)
    type (vector_t), allocatable :: vector_evec(:)
    type (vector_t) :: work(3)
    logical, parameter :: vector_is_real = .true.
    real (dp) :: eval_k, var_k, norm_k
    logical :: check_previous_eigenvectors

    kmin = 1
    kmax = nevec
    check_previous_eigenvectors = .false.
    if (present(range)) then
       check_previous_eigenvectors = .true.
       if (present(nocheck)) check_previous_eigenvectors = .not.nocheck
       kmin = max(1,range(1))
       kmax = min(nevec, range(2))
       range(1) = 1
       range(2) = 0
    endif

! if seed > 0, use it to create iseed(4) for LAPACK random number procedure dlarnv
    seed0 = seed
    if(seed0 > 0) then
       do  k = 1,4
          iseed(k)= mod(seed0,4096)
          seed0 = seed0/4096
       enddo
       if(mod(iseed(4),2)==1) iseed(4) = iseed(4) - 1
    endif                                                                         

    if(report>0) write(6,20) nevec, n, maxstep, report
20  format(' RS_LANCZOS: lowest', i5,'  eigenvectors of real symmetric matrix, dim=',i16,&
         &/,' max lanczos steps = ',i12,'; progress report step frequency:',i5)

    allocate (vector_evec(nevec))
    do k = 1, nevec
       vector_evec(k)%is_real = vector_is_real
       vector_evec(k)%n = n
       vector_evec(k)%r => evec(:,k)
       nullify(vector_evec(k)%c)
    enddo

    do k = 1,3
       call vector_create(n, vector_is_real, work(k))
    enddo

    if (check_previous_eigenvectors) then
       do k = 1, kmin - 1
          call get_residual(vector_evec(k), work(1), matvec, eval_k, var_k, norm_k)
          if (eval_k /= eval(k) .or. var_k /= variance(k)) then
             write(6,'("eigenvector ",i5," mismatch: norm = ",f25.16)') norm_k
             write(6,'("supplied eigenvalue, variance =",f25.16,1es12.3)') &
                  eval(k), variance(k)
             write(6,'("computed eigenvalue, variance =",f25.16,1es12.3)') &
                  eval_k, var_k
             stop
          else if (report > 0) then
             write(6,'("verified eigenvector ",i6," norm =",f25.16)') k, norm_k
          endif
       enddo
    endif

    resolution = 0_dp
    do k = kmin, kmax
       info(1) = 0
       if(report > 0) info(1) = -report
       if (report > 0) write(6,'(" finding eigenvector ",i5)') k
       call lanczos_iterator(n,vector_is_real,k-1,vector_evec,eval,variance,&
            resolution,matvec,iseed,maxstep,info,work)
       if (present(range)) range(2) = k
       if (report > 0) then
          write(6,'(" eigenvalue ",i4,":",f25.16,/" residual norm",1es12.3,&
               " (resolution limit:",1es12.3,"; ",&
               & i8," regular +",i5," polishing Lanczos steps)")') k,eval(k),&
               variance(k),resolution,info(1),info(2)

       endif
    enddo

    do k = 1,3
       call vector_erase(work(k))
    enddo

    do k = 1, nevec
       nullify(vector_evec(k)%r)
    enddo
    deallocate (vector_evec)
    return
  end subroutine rs_lanczos
