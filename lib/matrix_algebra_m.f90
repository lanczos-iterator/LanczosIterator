module matrix_algebra_m

contains
  subroutine lowest_tmatrix_eigenvalue(n,d,e,e_min)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: d(n), e(n)
    double precision, intent(out) :: e_min
    
    !-----------------------------------------------------------------
    !   T(i,i)    = D((i) , i =  1 : N
    !   T(i,i+1)  = T(i+1,i)  = E(i)  ,  i = 1 : N-1
    !   E_MIN is lowest eigenvalue  of real symmetric tridiagonal matrix  
    !-----------------------------------------------------------------
    
    double precision, allocatable :: w(:), work(:)
    integer, allocatable :: iblock(:), isplit(:), iwork(:)
    double precision :: v
    integer :: m, nsplit, info
    
    interface  ! to LAPACK subprogram DSTEBZ
       subroutine dstebz(range, order, n, vl, vu, il, iu, abstol, d ,e ,m ,&
            & nsplit, w, iblock, isplit, work, iwork, info)
         character, intent (in) :: order, range
         integer, intent(in)  :: n, il, iu, info
         integer, intent(out)  ::  m, nsplit
         double precision, intent(in)  :: abstol, vl, vu
         integer, intent(out) ::  iblock(*), isplit(*), iwork(*)
         double precision, intent(in) :: d(*), e(*)
         double precision, intent(out) :: w(*), work(*)
       end subroutine dstebz
    end interface
    
    double precision :: dlamch
    external dlamch
    
    logical , save :: first_call= .true.
    double precision,  save :: abstol
    if (first_call) then
       abstol = 2*dlamch('S')         
       first_call = .false.
    endif
    
    allocate (w(n),iblock(n),isplit(n),work(4*n),iwork(3*n))
    call dstebz('I','E',n,v,v,1,1,abstol,d,e,m,nsplit,w,iblock,isplit,work,iwork,info)
    
    if(info /= 0 .or.  m /= 1) then
       write(6,'("LOWEST_TMATRIX_EIGENVALUE: DSTEBZ error, m=",i5," info =",i12)') m,info
       stop
    endif
    
    e_min = w(1)
    deallocate(w,iblock,isplit,work,iwork)
  
    return
  end subroutine lowest_tmatrix_eigenvalue
  
  subroutine tridiagonalize(n,ab,d,e,q)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: ab(3,n)
    double precision, intent(out) :: d(n), e(n), q(n,n)
    
    ! reduces a real symmetric band matrix with two off-diagonals to a 
    ! real symmetric tridiagonal matrix
    !   Q**T * A * Q = T     (Q is an orthogonal full matrix)
    !
    ! A(i,j) is input as 
    ! AB(1,) = A(i,i), AB(2,i) = A(i,i+1)), AB(3,i) = A(i,i+2)
    ! AB(3,n-1), AB(2,n), AB(3,n) are indeterminate
    ! ------------------------------------------------------------------ 
    ! D(i) = T(i,i),   i = 1:N
    ! E(n) = T(i,i+1), i = 1:N-1 
    !-----------------------------------------------------------------
    integer :: info
    double precision, allocatable :: w(:)
    interface  ! to LAPACK subprogarm DSBTRD
       subroutine dsbtrd(vect, uplo,n,kd,ab,ldab,d,e,q,ldq,work,info)
         character*1, intent(in) :: vect, uplo
         integer, intent(in) :: n, kd, ldab, ldq
         double precision, intent(inout) :: ab(ldab,*), q(ldq,*)
         double precision, intent(out)   :: d(*), e(*), work(*)
         integer, intent(out) :: info
       end subroutine dsbtrd
    end interface
  
    allocate (w(n))
  
    ! LAPACK subroutine DSBTRB transforms a band matrix to tridiagonal form
    !  using plane rotations  (orthogonal transformation)
    call dsbtrd('V','L',n,2,ab,3,d,e,q,n,w,info)
    
    if(info /= 0) then
       write(6,'("TRIDIAGONALIZE: DSBTRD, INFO=",i12)') info
       stop
    endif
  
    deallocate (w)
    
    return
  end subroutine tridiagonalize
end module matrix_algebra_m
