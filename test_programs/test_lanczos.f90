! simple sample ch_matvec for testing 
subroutine my_ch_matvec(n, vec_in, vec_out, add)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent (in) :: n
  complex (dp), intent(in) :: vec_in(n)
  complex (dp), intent(inout) :: vec_out(n)
  logical, intent(in) :: add
  integer :: i
  if (add) then
     forall (i = 1:n) vec_out(i) = vec_out(i) + vec_in(i) * i
  else
     forall (i = 1:n) vec_out(i) = vec_in(i) * i
  endif
end subroutine my_ch_matvec

! simple sample rs_matvec for testing 
subroutine  my_rs_matvec(n, vec_in, vec_out, add)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent (in) :: n
  real (dp), intent(in) :: vec_in(n)
  real (dp), intent(inout) :: vec_out(n)
  logical, intent(in) :: add
  integer :: i
  if (add) then
     forall (i = 1:n) vec_out(i) = vec_out(i) + vec_in(i) * i
  else
     forall (i = 1:n) vec_out(i) = vec_in(i) * i
  endif
end subroutine my_rs_matvec

program main
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer::  n, nevec, maxstep, report, seed
  real (dp) :: resolution
  real (dp), allocatable :: eval(:), variance(:)
  real (dp), allocatable :: evec_rs(:,:) 
  complex (dp), allocatable :: evec_ch(:,:)
  integer ::  i, use_blas, wrap_zdotc, matrix_type, n1
  integer :: range(2)  ! optional arguments to rs_lanczos, ch_lanzcos
  logical :: nocheck   ! optional arguments to rs_lanczos, ch_lanzcos
  external my_rs_matvec, my_ch_matvec

  ! the interface is needed if there are (as here) some calls to (rs,ch)_lanczos that
  ! omit  optional variables, and some that do not.
  interface
     subroutine rs_lanczos(n, nevec, eval, evec, variance, resolution,&
          matvec, maxstep, report, seed, range, nocheck)
       implicit none
       integer, parameter :: dp = kind(1.0d0)
       integer, intent(in) :: n, nevec, maxstep, report, seed
       real (kind=dp), intent(out) :: eval(nevec), variance(nevec)
       real (kind=dp), intent(out), target :: evec(n,nevec)
       real (kind=dp), intent(out)  :: resolution
       integer, intent(inout), optional :: range(2)
       logical, intent(in), optional :: nocheck
       external  matvec
     end subroutine rs_lanczos
  end interface

  interface
     subroutine ch_lanczos(n, nevec, eval, evec, variance, resolution,&
          matvec, maxstep, report, seed, range, nocheck)
       implicit none
       integer, parameter :: dp = kind(1.0d0)
       integer, intent(in) :: n, nevec, maxstep, report, seed
       real (kind=dp), intent(out) :: eval(nevec), variance(nevec)
       complex (kind=dp), intent(out), target :: evec(n,nevec)
       real (kind=dp), intent(out)  :: resolution
       integer, intent(inout), optional :: range(2)
       logical, intent(in), optional :: nocheck
       external  matvec
     end subroutine ch_lanczos
  end interface
  
  print '("Lanczos test program in FORTRAN")'
  print '("  n = dimension of matrix")'
  print '("  nevec = number of eigenvalues requested")'
  print '("  maxstep is number of Lanczos steps allowed")'
  print '("  maxstep is number of Lanczos steps allowed (0 for no limit)")'
  print '("  give report = r > 0 to see progress every r Lanczos steps")'
  print '("  (plus summary of Lanczos process for each eigenvector),1x")'
  print '("give n, nevec, maxstep, report")'
  read *, n, nevec, maxstep, report

  if (nevec > n) nevec = n
  seed = 12345678
  
  allocate (eval(nevec), variance(nevec))
  
  print '(" This program can use BLAS = Basic Linear Algebra Suprograms")'
  print '("enter 1 to use BLAS, 0 for no BLAS")'
  read *, use_blas
  select case (use_blas)
  case (0, 1)
     call set_use_blas(use_blas)
  case default
     print '("invalid choice")'
     stop
  end select
  
  print '(/"The test matrix is  diagonal with eigenvalues (0,1,2,3,....)")'
  print '("(The Lanczos algorithm is not aware that the matrix is diagonal)",/1x)'
  
  print '("give 1 for real symmetric, 2 for complex Hermitian")'
  read  * ,matrix_type

  select case (matrix_type)
  case (1)
     allocate (evec_rs(n, nevec)) 
     ! examples of the two computational modes: 
     ! compute the first n1 eigenvectors in all-in-one-go mode (range and
     ! nocheck not present)
     n1 = 1 + nevec/2
     seed = modulo(17*seed,7654321)
     call rs_lanczos(n, n1,  eval,  evec_rs,  variance, &
          resolution, my_rs_matvec, maxstep, report, seed)

     ! compute the remaining eigenvectors in one-by-one mode (range(1)=range(2))
     ! when (range(1),range(2)) is specified,the validity of the previously-found
     ! eigenvectors will by default (nocheck= .false.)  be checked by computing their
     ! eigenvalue and variance and comparing with the values in eval and variance.

     nocheck = .false.   ! check that the previous eigenvector data is correct     
     do i = n1 + 1, nevec
        range = i
        seed = modulo(17*seed,7654321)
        call rs_lanczos(n, i,  eval,  evec_rs,  variance, &
             resolution, my_rs_matvec, maxstep, report, seed, range, nocheck)
        nocheck = .true.  ! no point in re-checking
     enddo

  case (2)
     if (use_blas == 1) then
        print '("enter 0 to use BLAS without a ZDOTC wrapper, 1 for wrapper")'
        read *, wrap_zdotc
        select case  (wrap_zdotc) 
        case (0, 1)
           call set_wrap_zdotc(wrap_zdotc == 1)
        case default
           print '("invalid choice")'
           stop
        end select
     end if
     allocate (evec_ch(n, nevec))
     ! examples of the two computational modes: 
     ! compute the first n1 eigenvectors in all-in-one-go mode (range and
     ! nocheck not present)
     n1 = 1 + nevec/2
     seed = modulo(17*seed,7654321)
     call ch_lanczos(n, n1,  eval,  evec_ch,  variance, &
          resolution, my_ch_matvec, maxstep, report, seed)

     ! compute the remaining eigenvectors in one-by-one mode (range(1)=range(2))
     ! when (range(1),range(2)) is specified,the validity of the previously-found
     ! eigenvectors will by default (nocheck= .false.)  be checked by computing their
     ! eigenvalue and variance and comparing with the values in eval and variance.

     nocheck = .false.   ! check that the previous eigenvector data is correct
     do i = n1 + 1, nevec
        range = i
        seed = modulo(17*seed,7654321)
        call ch_lanczos(n, i,  eval,  evec_ch,  variance, &
             resolution, my_ch_matvec, maxstep, report, seed, range, nocheck)
        nocheck = .true.  ! no point in re-checking
     enddo
  case default
     print '("invalid choice")'
     stop
  end select

  do i = 1,nevec
     print '(i10, f25.16," variance =", es10.3)', i, eval(i), variance(i)
  enddo
  print '("  resolution: ",es10.3)', resolution
  print '("  eigenvalue range: [",i0,":",i0,"]")', range
  deallocate (eval, variance)
  if (allocated(evec_rs)) deallocate(evec_rs)
  if (allocated(evec_ch)) deallocate(evec_ch)
  stop
end program
