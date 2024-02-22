subroutine rs_lanczos_wrapper(n_c, nevec_c, eval_cp, evec_cp, variance_cp, &
     resolution_cp, rs_matvec_from_c, maxstep_c, report_c, seed_c, range_cp, nocheck_c) &
     bind(c,name='rs_lanczos_c')
  use, intrinsic :: iso_c_binding
  implicit none
  !--------------------------------------------------------------------
  ! This is a wrapper for FORTRAN rs_lanczos to allow calling from C
  !  *_cp are C pointers (double) for output of  results back to the calling C program
  !  *_c are C integers (int) providing input parameters for the FORTRAN program
  !-------------------------------------------------------------------
  integer (c_int), intent(in), value :: n_c, nevec_c, maxstep_c, report_c, seed_c
  logical (c_bool), intent(in), value :: nocheck_c
  type (c_ptr), value, intent(in) :: range_cp
  type (c_ptr), value, intent(in) :: eval_cp, evec_cp, variance_cp, resolution_cp
  type (c_funptr), value, intent(in) :: rs_matvec_from_c
  type (c_funptr), save :: rs_matvec_cp = c_null_funptr
  integer, parameter :: dp = kind(1.0d0)

  interface
     subroutine rs_lanczos(n, nevec, eval, evec, variance, resolution,&
          matvec, maxstep, report, seed, range, nocheck)
       implicit none
       integer, parameter :: dp = kind(1.0d0)
       integer, intent(in) :: n, nevec, maxstep, report, seed
       real(dp), intent(out) :: eval(nevec), variance(nevec)
       real(dp), intent(out), target :: evec(n,nevec)
       real(dp), intent(out)  :: resolution
       integer, intent(inout), optional :: range(2)
       logical, intent(in), optional :: nocheck
       external  matvec
     end subroutine rs_lanczos
  end interface
  
  integer :: n, nevec, maxstep, report, seed
  real (dp), pointer :: eval_p(:), variance_p(:), resolution_p(:)
  real (dp), pointer :: evec_p(:,:)
  integer, pointer :: range_p(:)
  logical :: nocheck
  
  rs_matvec_cp = rs_matvec_from_c
  n = n_c
  nevec = nevec_c
  
  call c_f_pointer (eval_cp, eval_p,(/ nevec /))
  call c_f_pointer (variance_cp, variance_p, (/ nevec /)) 
  call c_f_pointer (evec_cp, evec_p, (/ n, nevec /))
  call c_f_pointer (resolution_cp, resolution_p, (/ 1 /))
  maxstep = maxstep_c
  report = report_c
  seed = seed_c

  if (c_associated(range_cp)) then
     call c_f_pointer (range_cp, range_p, (/ 2 /))
     if (nocheck_c) then
        nocheck = .true.
     else
        nocheck = .false.
     endif
     call rs_lanczos(n, nevec, eval_p, evec_p, variance_p, resolution_p(1), &
          rs_matvec_wrapper, maxstep, report, seed, range_p, nocheck)
  else
     call rs_lanczos(n, nevec, eval_p, evec_p, variance_p, resolution_p(1), &
          rs_matvec_wrapper, maxstep, report, seed)
  endif

  return
contains
  subroutine rs_matvec_wrapper(n, vec_in, vec_out, add) 
    use iso_c_binding
    integer, parameter :: dp = kind(1.0d0)
    integer, intent(in) :: n
    real (dp), intent(in) :: vec_in(n)
    real (dp), intent(out) :: vec_out(n)
    logical, intent(in) :: add
    !-----------------------------------
    ! /* C code: */
    ! int *create_real_storage_c(int n) {
    !    return (double *) malloc(n * sizeof(double));
    ! }
    ! void destroy_real_storage_c(double *p) {
    !    free(p);
    !    p = NULL;
    ! }
    !------------------------------------------,
    
    interface
       function create_real_storage_c(n) result (p) bind(c)
         use, intrinsic :: iso_c_binding, only : c_ptr, c_int
         integer (c_int), intent(in), value :: n
         type (c_ptr) :: p
       end function create_real_storage_c
    end interface
    
    interface
       subroutine destroy_real_storage_c(p) bind(c)
         use, intrinsic :: iso_c_binding, only : c_ptr
         type (c_ptr), intent(in), value :: p
       end subroutine destroy_real_storage_c
    end interface
    
    !------------------------------------------,
    ! rs_matvec_c must be equivalent to
    !
    ! double matrix(int i, int j); /*returns matrix element */
    !
    ! void rs_matvec_c(int n, double *vec) {
    !   double  *mat_vec = (double *) calloc(n, sizeof(double))
    !   for (int i = 0; i < n; i++) {
    !      for (int j = 0; j < n; j++) {
    !          mat_vec[i] += matrix(i,j) * vec[j]
    !      }
    !   }
    !   memcpy(vec, mat_vec, n * sizeof(double))
    !   free (mat_vec)
    ! }
    !
    !-------------------------------------
    
    abstract interface
       subroutine rs_matvec_c(n, p) bind(c)
         use, intrinsic :: iso_c_binding, only : c_int, c_ptr
         integer (c_int), intent(in), value :: n
         type (c_ptr), intent(in), value :: p
       end subroutine rs_matvec_c
    end interface
    
    real (dp), pointer :: vec_pf(:)
    type (c_ptr) :: vec_pc
    procedure(rs_matvec_c), pointer :: rs_matvec_p
    
    vec_pc = create_real_storage_c(n)
    call c_f_pointer(vec_pc, vec_pf, (/ n /))  
    vec_pf = vec_in
    
    call c_f_procpointer (rs_matvec_cp, rs_matvec_p)
    call rs_matvec_p(n, vec_pc)
    
    if (add) then
       vec_out = vec_out + vec_pf
    else
       vec_out = vec_pf
    endif
    
    call destroy_real_storage_c(vec_pc)
    return
  end subroutine rs_matvec_wrapper
end subroutine rs_lanczos_wrapper
