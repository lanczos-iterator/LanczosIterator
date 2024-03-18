! (c) F. D. M. Haldane, Princeton University, haldane@princeton.edu
! 2014-07-31   For use with lanczos_iterator subprogram,
!-------------------------------------------------------------------------
module tmatrix_m
  integer, parameter, private :: dp = kind(0d0)  ! double precision
  private u_inv,d_inv,l_inv,reverse_uinv,reverse_l_inv
  private cholesky, reverse_cholesky
  real(kind=dp), parameter, private ::  zero = 0.0d0, one = 1.0d0
contains
  
  
  subroutine  tmatrix_eval(n,tmat,e0)
    use matrix_algebra_m       
    implicit none
    integer, intent(in)::  n
    real (kind=dp), intent(in) :: tmat(2,n)
    real (kind=dp), intent(inout) :: e0(n)
    
    
    real (kind=dp), allocatable :: d(:),e(:),w(:),work(:)
    integer, allocatable :: iblock(:),isplit(:),iwork(:)
    integer :: m,nsplit,info,i
    real(kind=dp)::v
    
    allocate (d(n),e(n))
    d(1) = tmat(1,1)
    do i = 2,n
       d(i) = tmat(1,i)
       e(i-1) = tmat(2,i-1)
    enddo
    
    
    if(n == 1) then
       e0(1) = tmat(1,1)
       return
    elseif (n == 2) then
       e0(2) = (tmat(1,1) + tmat(1,2))/2 - dsqrt( ((tmat(1,1)-tmat(1,2))/2)**2 + tmat(2,1)**2)
       return
    endif
    
    call lowest_tmatrix_eigenvalue(n,d,e,e0(n))    ! in  module matrix_algebra_m

    deallocate(d,e)
    
    return
  end subroutine tmatrix_eval
  
  
  subroutine tmatrix_evec(n,tmat,ev,wf,var)
    implicit none
    integer, intent(in):: n
 
    real(kind=dp), intent(in):: tmat(2,n)
    real(kind=dp),intent(inout):: ev
    real(kind=dp), intent(out)::wf(n)
    !-----------------------------------------------
    ! obtain an approximation to the eigenvector corresponding
    ! to the lowest eigenvalue of the t-matrix, which is closely
    ! approximated by  ev
    !-------------------------------------------------      
    integer :: i
    real(kind=dp):: var,d_ev,var_in,rescale
    real(kind=dp), allocatable:: w(:)
    
    allocate (w(n))
    
    call tmatrix_cholesky_evec(n,tmat,ev,d_ev,var,w)
    ev = ev + d_ev
    
    ! one step of rayleigh quotient iteration
    call tmatrix_rqi(n,tmat,w,wf,ev,d_ev,var_in, var,rescale)
    ev = ev + d_ev
    
    deallocate(w)
    
    return
  end subroutine tmatrix_evec
  
  subroutine tmatrix2_evec(n,tmat,ev,var,wf)
    implicit none
    integer, intent(in):: n
    real(kind=dp), intent(in):: tmat(2,n)
    real(kind=dp), intent(inout):: ev
    real(kind=dp), intent(out):: var,wf(n)
    !--------------------------------------------
    !
    !  on input, tmat(2,n) contains the Lanczos t-matrix
    !  after n steps, and ev is a converged  approximation 
    !  to the lowest eigenvalue of T.
    !
    ! on output, wf(1),...,wf(n) is a "best" choice for the
    ! eigenvector
    !
    ! |EVEC> = sum_{i=1,n} wf(i)|I>
    !
    ! to minimize its variance.
    ! on output, EV and VAR are the eigenvalue and variance
    ! (residual nrm) based on an assumption that orthogonality
    !  <I|J> = delta_(I,J) remains exact.   In practice this
    ! is not so, and these quantitites must be recalculated using
    ! |EVEC>, insted of  WF(1),...,WF(N).
    !
    ! uses LAPACK subroutine DSBTRD to reduce a real symmetric band matrix to
    ! tridiagonal form.
    !----------------------------------------------
    integer ::i
    real(kind=dp):: norm,eval,eval1
    
    allocatable tmat2,q,w,cmat,u
    double precision tmat2(:,:),q(:,:),w(:),cmat(:,:),u(:)
    
    allocate (w(n))
    allocate (tmat2(2,n),q(n,n),cmat(2,n),u(n))
    
    ! set up tridiagonal representation  of band matrix P(1:n)((T-ev)**2)P(1:n)
    ! with (0,0,...,1) as first basis state |T2(1)>
    call reverse_tmatrix2_treduce(n,tmat,ev,tmat2,q)
    
    !  note, the orthogonal matrix q(i,j) is <n+1-i|T2(j)>
    
    ! Cholesky decomposition P_TTP_n = LDL^t
    
    call cholesky(n,tmat2,zero,cmat)
    !
    ! cmat(1,n) should be close to zero!
    ! 
    
    !  <T2(i)| P_n((T-ev)**2)P_n|T2(j)>  = LDL^t, with
    !
    !   D(i) = diag( cmat(1,i)), cmat(1,n) = zero (or close to zero)
    !   L(i,i) = 1, L(i+1,i) = cmat(2,i),i = 1,n-1,    L(i,j) = 0 ortherwise.
    !-----------------------------------
    !  solve sum_j L(j,i)*d(j)  =  0 (i<n),
    !        sum_j L(j,n)*d(j)  =  1
    
    u(1:n-1) = zero
    u(n) = one
    call u_inv(n,cmat,u)
    
    ! solve wf(n+1-i) = w(i)  = sum_j Q(i,j)*d(j),   sum_j L(j,i)*d(j) =   delta_(i,n).   
    
    w = zero
    do i = 1,n
       w = w + u(i)*q(:,i)
    enddo
    
    deallocate (tmat2,q,cmat,u)
    
    ! reverse the basis to obtain wf(1),..,wf(n)
    ! in original Lanczos t-matrix basis order.
    
    forall(i=1:n) wf(i) = w(n+1-i)
    
    if(wf(1) < zero) wf = -wf
    norm = sqrt(sum(wf**2))
    wf = wf/norm
    
    ! compute the predicted variance, assuming  exact arithmetic.
    
    ! act with tridiagonal matrix : w = tmat*wf
    call tmatvec(n,tmat,zero,wf,w)
    eval = sum(wf * w)
    w = w - eval1 * wf

    ! repeat, for full,accuracy.
    eval1 = sum(wf * w)
    w = w - eval1 * wf
    ! var is the predicted variance (residual norm)
    var = sqrt(sum(w**2))
    eval=  eval + eval1
    ! eval should be close to the input value ev.
    deallocate (w)
    ev = eval
    
    return
  end subroutine tmatrix2_evec
  
  subroutine reverse_tmatrix2_treduce(n,tmat,shift,tmat2,q)
    use matrix_algebra_m
    implicit none
    integer, intent(in)::  n
    real(kind=dp),intent(in):: tmat(2,n),shift
    real(kind=dp),intent(out):: tmat2(2,n),q(n,n)

    !-------------------------------------------------
    !  P(1:n)((T-shift)**2)P(1:n) = QT'Q^t
    !   T' is the trodiagonal matrix T'(i,i) = d(i), T'(i,i+1) = T'(i+1,i) = e(i).
    !
    !  Q(1:n,1) = (0,0,0,....,0,1)
    !
    !----------------------------------------------
    
    real(kind=dp),allocatable:: ab(:,:),d(:),e(:)
    integer:: i
    
    allocate (ab(3,n),d(n),e(n-1))
    
    d = tmat(1,:) - shift
   
    ! use a  reverse-order  basis for cholesky-decomposition stability of (T-ev)**2
    ! we will construct a tridiagonal representation of (T-ev)**2 projected into the Lanczos basis
    ! |1>, |2>,.., |n>, with initial vector |T2(1)>  = |n>
    !
    !   P_n((T-ev)**2)P_n|T2(j)> = d(j)|T2(j)> + e(j)|T2(j+1)> + e(j-1)|T2(j-1)>,   e(0) = 0.
    !  
    !  where P_n|j> = |j> for j <= n, P_n|j> = 0 for j > n.
    
    ab(1,n) = d(1)**2 + tmat(2,1)**2
    forall(i=1:n-1) ab(1,n-i) =   d(i+1)**2 + tmat(2,i+1)**2 + tmat(2,i)**2
    forall(i=1:n-1) ab(2,n-i)   = tmat(2,i)*(d(i+1) + d(i))
    forall(i=1:n-2) ab(3,n-i-1) = tmat(2,i)*tmat(2,i+1)
    

    !reduce the band matrix ab  to tridiagonal form (d,e)
    call tridiagonalize(n,ab,d,e,q)   ! in  module matrix_algebra_m
    
    tmat2(1,1:n)   = d
    tmat2(2,1:n-1) = e
    
    deallocate (ab,d,e)
    
    return
  end subroutine reverse_tmatrix2_treduce
  
  subroutine tmatrix_cholesky_evec(n,tmat,ev,d_ev,var,wf)
    implicit none
    integer, intent(in)::  n
    real(kind=dp), intent(in):: tmat(2,n),ev
    real(kind=dp), intent(out):: d_ev,wf(n),var
    !------------------------------------------
    ! if ev is close to the lowest eigenvalue of the t-matrix,
    ! wf(1:n) is the approximate eigenvector, ev + d_ev is
    ! its rayleigh quotient, and var is its residual norm.
    !----------------------------------------------
    real(kind=dp), allocatable:: cmat(:,:),w(:)
    real(kind=dp):: dd_ev,norm
    allocate (cmat(2,n),w(n))
    
    ! cholesky decomposition, in reverse order
    call  reverse_cholesky(n,tmat,ev,cmat) 
    
    ! cmat(1,n) should be close to zero, if ev is close to the lowest eigenvalue
    ! of tmat.
    
    !      wf(1) = one
    !      do i = 2,n
    !         wf(i) = -cmat(2,n+1-i)*wf(i-1)
    !      enddo
    
    wf(1) = one
    wf(2:n) = zero
    call reverse_u_inv(n,cmat,wf)
    
    norm = sqrt(sum(wf**2))
    wf = wf / norm
    if(wf(1) < zero) wf = -wf
    
    ! act with tridiagonal matrix : w = tmat*wf
    call tmatvec(n,tmat,ev,wf,w)
    
    ! d_ev is the predicted eigenvalue relative to ev_in
    d_ev = sum(wf*w)
    w = w - d_ev * wf
    ! var is the predicted variance (residual norm)
    if(var < abs(d_ev)) then
       dd_ev = sum(wf*w)
       w = w - dd_ev * wf
       var = sqrt(sum(w**2))
       d_ev = d_ev + dd_ev
    endif
    var = sqrt(var**2 + (tmat(2,n)*wf(n))**2)
    deallocate (w,cmat)
    
    return
  end subroutine tmatrix_cholesky_evec
  
  subroutine tmatrix_rqi(n,tmat,wf_in,wf_out,ev_in,d_ev,var_in, var_out,rescale)
    implicit none
    integer, intent(in)::n
    real(kind=dp), intent(in):: tmat(2,n),wf_in(n)
    real(kind=dp), intent(out):: ev_in, wf_out(n), var_in,var_out,rescale,d_ev
     !------------------------------------------
    ! rayleigh quotient iteration
    !   given wf_in (not necessarily normalized) find:
    !
    !   Its rayleigh quotient and residual norm !    ev_in =   (wf_in,T*wf_in)/(wf_in,wf_in)
    !    var_in = |(T-ev_in)*wf_in|/|wf_in|
    !
    !    Then solve  (T-ev_in)*y =   (wf_in/|wf_in|),
    !
    !     wf_out = rescale*y,  |wf_out| = 1    
    !
    !     rescale => 0 if wf_in becomes an  exact eigenvector
    !
    !    The change in Rayleigh quotient is 
    !    d_ev =  (wf_out,(T-ev_in)*wf_out)
    !
    !    var_out = |(T-ev_out)*wf_out|,   ev_out = ev_in + d_ev
    !------------------------------------------------
    real(kind=dp),allocatable :: cmat(:,:), w(:)
    real(kind=dp):: dd_ev,norm
    integer:: i
    
    allocate (cmat(2,n),w(n))
    
    norm = sqrt(sum(wf_in**2))
    rescale = 1/norm
    ! act with tridiagonal matrix : w = tmat*wf
    call tmatvec(n,tmat,zero,wf_in,w)
    ev_in = sum(wf_in * w)/(norm**2)
    call tmatvec(n,tmat,ev_in,wf_in,w)
    dd_ev = sum(wf_in * w)/(norm**2)
    w = w - dd_ev * wf_in

    ev_in = ev_in + dd_ev

    var_in = sqrt(sum(w**2))/norm
  
    ! form T-ev_in = LDU, U = L^t    
    call  reverse_cholesky(n,tmat,ev_in,cmat)
    
    !      if(abs(cmat(1,n)).gt.1.0d-4) then
    !         write(6,102)
    ! 102     format(' TMATRIX_RQI error')
    !         stop
    !      endif
    
    ! cmat(1,n) should be close to zero, if ev_in is close to the lowest eigenvalue
    ! of tmat.   Reverse order and the Lanczos structure guarantees that the approximate zero occurs
    ! only in the last position D(n) of the (T-ev) = LDT^t decomposition.
    
    wf_out = wf_in
    if (cmat(1,n) == zero) then   ! we can't proceed, too close to an eigenvector!
       rescale = zero
       var_out = var_in
       d_ev = zero
       return
    endif

    call l_inv(n,cmat,wf_out)
    rescale = rescale*cmat(1,n)
    wf_out(1:n-1) = wf_out(1:n-1)*cmat(1,n)
    forall(i=1:n-1) wf_out(i) = wf_out(i)/cmat(1,i)
    call reverse_u_inv(n,cmat,wf_out)
    norm = sqrt(sum(wf_out**2))
    rescale = norm*rescale
    wf_out = wf_out/norm

    ! Rayleigh quotient iteration is finished when rescale is small enough
    
    call tmatvec(n,tmat,ev_in,wf_out,w)
    d_ev = sum(wf_out * w)
    w = w  - d_ev * wf_out
    var_out = sqrt(sum(w**2))
    if(var_out < abs(d_ev)) then
       dd_ev = sum(wf_out*w)
       w = w - dd_ev * wf_out
       d_ev = d_ev + dd_ev
       var_out = sqrt(sum(w**2))
    endif
    var_out = sqrt(var_out**2 + (tmat(2,n)*wf_out(n))**2)
    deallocate (w,cmat)
    
    return
  end subroutine tmatrix_rqi
  
  integer  function tmatrix_slice(n,tmat,t)
    implicit none
    integer, intent(in):: n
    real(kind=dp), intent(in):: tmat(2,n),t
    !--------------------------------
    ! slice = number of eigenvalues of the real symmetric
    ! n x n  tridiagonal matrix tmat(2,n)  which  are less
    ! than t.
    !
    !  T(i,i) = TMAT(1,i)  i = 1,..,n
    !  T(i,i+1) = T(i+1,i) = TMAT(2,i) > 0, i = 1,..,n-1
    !  T(i,j) = 0  |i-j| > 1.
    !
    !-------------------------------------
    real(kind=dp):: y,x,shift,emax,x_min, small, rnd, eps
    logical ::restart
    integer :: m
    logical, save :: first_call = .true.
    real(kind=dp), save :: sfmin
    real(kind=dp), parameter :: one = 1_dp, zero = 0_dp
    integer :: expon 

    if(first_call) then      ! get overflow threshold 1/sfmin (based on LAPACK reference dlamch)
       rnd = one
       if(rnd == one) then   ! assume rounding not chopping
          eps = epsilon(zero) * 0.5
       else
          eps = epsilon(zero)
       endif
       sfmin = tiny(zero)
       expon = 1 - maxexponent(zero) - minexponent(zero)
       ! LAPACK dlamch reference uses small = one/huge(zero), which causes  underflow if expon < 0
       ! one/huge(zero) = tiny(zero) * (radix(zero)**expon) / (one - epsilon(zero)/digits(zero))
       if (expon < 0) then
          small = zero
       else  if (expon > 0) then
          small = sfmin * radix(zero)**(abs(expon))
       else
          small = sfmin
       endif
       if (small >= sfmin) then        !make sfmin just a bit bigger 
          sfmin = small * (one + eps)  !than small, to avoid
       endif                           !overflow due to rounding
       first_call = .false.
    endif
    
    if(n < 0) then
       write(6,'("TMATRIX_SLICE: called with n = ",i12)') n
       stop
    endif
    
    select case (n)
    case(1)
       if(t > tmat(1,1)) then
          tmatrix_slice = 1 
       else
          tmatrix_slice = 0
       endif
    case(2)
       x = (tmat(1,1)-t)*(tmat(1,2)-1)  -tmat(2,1)**2
       if(x < zero) then
          tmatrix_slice = 1
          return
       else if(x > zero) then
          if(t < tmat(1,1)) then 
             tmatrix_slice = 0
          else
             tmatrix_slice = 2
          endif
          return
       else
          if(t < min(tmat(1,1),tmat(1,2))) then
             tmatrix_slice = 0
          else
             tmatrix_slice = 1
          endif
       endif
       return
    end select
    
    !      tmatrix_slice = -1
    !      return
    
    ! use cholesky LDU factorization for n > 2
    emax = zero
    tmatrix_slice = 0
    m = 1 
    restart = .true.
    shift = t
    do while (m <= n)
       
       !  use emax  as the scale factor that substitutes for   t(2,n).
       if(tmat(2,m) > emax) emax = tmat(2,m)
       
       if (restart) then
          if(m < n) then
             x =   (tmat(1,m)  - shift)/tmat(2,m)
          else
             x =   (tmat(1,m)  - shift)/emax
          endif
          restart = .false.
       else
          x = y/x
       endif
       !     x = x(m), y = y(m)
       !     if x(m) is close to zero, we may overflow when forming x(m+1) or y(m+2)
       !     if eiher m+1 or m+2 are less than n.   In this case jump from m to m+2,
       !     restart, and increase tmatrix_slice by 1.
       
       if(m+1 < n) then
          y = (tmat(1,m+1) - t)/tmat(2,m+1)
          y = (y*x) -(tmat(2,m)/tmat(2,m+1))
       else
          y = (tmat(1,m+1) - t)/emax
          y = (y*x) -(tmat(2,m)/emax)
       endif
       
       !     will x(m+1) = y(m+1)/x(m) cause overflow in the next step? 
       x_min = 2*abs(y)*sfmin
       !     will y(m+2) = (tmat(1,m+2)-t)/tmat(2,m+2))*x(m+1) cause overflow in the step after that?
       if(m+2 <= n) then
          if(m+2 < n) then
             if (abs(tmat(1,m+2)-t) > tmat(2,m+2))  x_min = (x_min*abs(tmat(1,m+2)-t))/tmat(2,m+2)
          else
             if (abs(tmat(1,m+2)-t) > emax)   x_min = (x_min*abs(tmat(1,m+2)-t))/emax
          endif
       endif
       if(abs(x) < x_min) then
          tmatrix_slice  = tmatrix_slice + 1
          restart = .true.
          if(m+2 <= n)  shift = (t*y+tmat(2,m+1)*x)/y
          m = m + 2
       else
          m = m + 1
          if(x < zero) tmatrix_slice = tmatrix_slice + 1
       endif
    enddo
      
    return
  end function tmatrix_slice
  
  subroutine lanczos_inverse_iteration(n,cmat,e2,reduction, vb,muc, work,lwork)
    implicit none
    integer, intent(in):: n,lwork
    real(kind=dp), intent(in):: cmat(2,n),e2
    real(kind=dp), intent(out):: muc,vb,reduction
    real(kind=dp), intent(inout):: work(lwork)
    !-----------------------------------------------------------
    ! solve T(n)P(1)T(n)|v> = mu(T(n)T(n) + e2*P(n))|v>
    !
    !   mu = 1-muc  ,  0< mu <= 1.
    !
    !  call lanczos_invit_getvec after this to get
    !   <i|v> returned in work(1),..work(n),
    !  normalized so <1|T|v> = 1.
    !
    !------------------------------------------
    !  here T is an nxn positive definite triadiagonal matrix,
    !  T|1> = d(1)|1> + e(1)|2>
    !  T|j+1> = e(j)|j> + d(j+1)|j+1> + e(j+1)|J+1>
    !
    !   T(n) is this matrix projected into the subspace spanned by 
    !  basis vectors |1>, |2>,...|n>   and e2  = e(n)**2.
    !  P(1) = |1><1|,  P(n) = |n><n|.
    !
    !  there is a unique solution with mu > 0, and has mu < 1 for non vanishing e2.
    !  mu = 1-muc.
    !
    !--------------------------------------------------
    ! T(n) is presented as the Cholesky decompositon
    !
    !  T(n) = L(n)D(n)L(n)^T
    !
    !   D(n) = diag(cmat(1,1),cmat(1,2),...,cmat(1,n)), 
    !   for all i = 1,2.,,.n, !mat(1,i) > 0
    !   positivity of D(n) guarantees T(n) is positive definite.
    !
    !   L(i,i) = 1, i = 1,..,n
    !   L(i+1,i) = cmat(2,i), i =1,...,n-1
    !   L(i,j) = 0 if (i-j) .ne. 0 or 1.
    !---------------------------------------
    !
    !   T|a> = |1>,    TT|b> = |n>
    !    |v> = |a> + vb|b>
    !    mu  = 1 +  vb*<1|T|b>
    !    vb = -e2(<N|a> + vb<N|b>)   = 0 if e2 = 0
    !    <1|T|v> = mu
    !
    !    T = LDL^t = LDU     
    !    L|x> = |N> -> |x> = |N
    !    U|x> = |1> -> |x> = |1>
    
    !   T|a> = |1>     
    !   L|x> = |1>
    !   D|y> = |x>
    !   U|a> = |y>
    !-------------------------------------
    !     variance reduction factor = dsqrt((1-p1)/p1)
    !     where  T|v> has projection  p1  on |1>
    !     
    !     p1 =  <v|TP(1)T|v>        = mu
    !     -------------------
    !     <v|TT + e2*P(n)|v>
    !
    !  with normalization <1|T|v> = 1
    !
    !                           sqrt(1-mu)
    !     reduction factor = ------------------------
    !                        || |0> - x\mu |v> ||
    !                                                            
    !     reduction factor = dsqrt(muc/mu)
    !----------------------------------------
    real(kind=dp):: mu,bt1,an,bn
    
    if(2*n > lwork) then
       write(6,10) lwork,2*n
10     format('LANCZOS_INVERSE_ITERATION: LWORK=',i12,' must be at least',i12)
       stop
    endif
    
    work(1) = one  
    work(2:n) = zero
    
    call l_inv(n,cmat,work)
    call d_inv(n,cmat,work)
    call u_inv (n,cmat,work)     
    
    if(e2 == zero) then
       !     solution is complete, and is already in work(1),....,work(n)
       reduction = zero
       vb = zero
       mu = one
       muc = zero
       return
    endif
    
    !     If e2 .ne. 0  also need solution of         
    !     TT|b> = |N>   LDULDL|b> = |N>
    !     ULDU|b> =   (1/D(N))|N>
    !     
    !   LDU|x> = |N>                                                                                                                                       
    !    DU|x> = |N>         
    !     U|x> = (1/D(N)|N>
    !     L|y> = |x>
    !     D|z> = |y>
    !     U|b> = |z>
    
    !     T|b>  =|x>   T|x>= |N>
    
    work(n+1:2*n-1) = zero
    work(2*n) = 1/cmat(1,n)
    call u_inv(n,cmat,work(n+1))
    
    !     work(n+1) currently contains <1T|b>
    !     save <1|T|b> = bt1
    bt1 = work(n+1)
    
    call l_inv(n,cmat,work(n+1))
    call d_inv(n,cmat,work(n+1))
    call u_inv(n,cmat,work(n+1))
    !  work(n+1:n+n) now contains |b> 
    !     solution is |v> = |a> - vb*|b>
    !     vb*(1+e2<N|b>) = (e2*<N|a>)
    !     mu = 1 + <1|T|b>(B/A) = <1|T|v>
    an = work(n)
    bn = work(2*n)
    vb = e2*an
    vb = vb/(one + e2*bn)
    muc = vb*bt1
    mu = one - muc
    reduction = sqrt(muc)
    return
  end subroutine lanczos_inverse_iteration
  
  subroutine lanczos_invit_getvec(n,vb,muc,work,lwork)
    implicit none
    integer, intent(in):: n, lwork
    real(kind=dp), intent(in):: vb,muc
    real(kind=dp), intent(inout):: work(lwork)
    !------------------------------------
    !     create the  inverse iteration solution in work.
    !     call after lanczos_inverse_iteration (with unchanged
    !     n,vb ,muc,work) to complete the solution produced there, 
    !     if the eigenvector is wanted.
    !     eigenvector is created in work(1),,,,work(n).
    !-----------------------------------------------
    real (kind=dp) ::  mu
    mu = one - muc
    work(1:n) = work(1:n) - vb * work(n+1:2*n)
    !     rescale so <1|T|v> = 1
    work(1:n) = work(1:n) / mu
    return
  end subroutine lanczos_invit_getvec
  
  subroutine tmatvec(n,tmat,shift,x,y)
    implicit none
    integer, intent(in):: n
    real(kind=dp),intent(in):: shift,tmat(2,n),x(n)
    real(kind=dp),intent(out):: y(n)
    !-------------------------
    !   y = (T-shift)*x
    !  where T is a  nxn real symmetric tridiaginal matrix
    ! with T(i,i) = tmat(1,i), T(i,i+1) = tmat(2,i)
    !----------------------------------------------
    integer i
    forall(i = 1:n)   y(i) = (tmat(1,i)-shift)*x(i)
    forall(i = 1:n-1) y(i) = y(i) + tmat(2,i)*x(i+1)
    forall(i = 2:n)   y(i) = y(i) + tmat(2,i-1)*x(i-1)
      
    return
  end subroutine tmatvec
  
  subroutine cholesky(n,tmat,ev,cmat)
    implicit none
    integer, intent(in):: n
    real(kind=dp),intent(in):: tmat(2,n),ev
    real(kind=dp),intent(out)::cmat(2,n)
    !----------------------------------
    ! cholesky decomposition of tmat - ev
    !-------------------------------    
    integer:: i

    cmat(1,1) = tmat(1,1)-ev
    do i = 1,n-1
       cmat(2,i) = tmat(2,i)/cmat(1,i)
       cmat(1,i+1) = (tmat(1,i+1)-ev) - cmat(2,i)*tmat(2,i)
    enddo
    
    return
  end subroutine cholesky
  
  subroutine u_inv(n,cmat,x)
    implicit none
    integer, intent(in):: n
    real(kind=dp), intent(in):: cmat(2,n)
    real(kind=dp), intent(out):: x(n)
    !--------------------------------------
    !     part of solver for T =LDL^t =  LDU
    !     solves Ux(out) = x(in)
    !     for i = 1...,n
    !     U(i,i) =  1
    !     U(i,i+1) = cmat(2,i)
    !-----------------------------
    integer:: i
    do i = n-1,1,-1
       x(i) = x(i) - cmat(2,i)*x(i+1)
    enddo
    return
  end subroutine u_inv
  
  subroutine l_inv(n,cmat,x)
    implicit none
    integer, intent(in)::n
    real(kind=dp), intent(in):: cmat(2,n)
    real(kind=dp), intent(out):: x(n)
    !--------------------------------------
    !     part of solver for T =LDL^t =  LDU
    !     solves Lx(out) = x(in)
    !     for i = 1...,n
    !     L(i,i) =  1
    !     L(i+1,i) = cmat(2,i)
    !--------------------------------
    integer:: i
    do i = 1,n-1
       x(i+1) = x(i+1) - cmat(2,i)*x(i)
    enddo
    return
  end subroutine l_inv
  
  subroutine d_inv(n,cmat,x)
    implicit none
    integer, intent(in)::n
    real(kind=dp), intent(in):: cmat(2,n)
    real(kind=dp), intent(out):: x(n)
    !--------------------------------------
    !     part of solver for T =LDL*=  LDU
    !     replaces x(in) by x(out), where  Dx(out) = x(in)
    !     for i = 1...,n
    !     D(i,i) = cmat(1,i)
    !--------------------------------
    integer:: i
    forall(i = 1:n) x(i) = x(i)/cmat(1,i)
    return
  end subroutine d_inv
  
  subroutine reverse_cholesky(n,tmat,ev,cmat)
    implicit none
    integer, intent(in):: n
    real(kind=dp),intent(in):: tmat(2,n),ev
    real(kind=dp),intent(out):: cmat(2,n)
    !-------------------
    ! cholesky decomposition, in reverse order
    !-----------------------
    integer:: i

    cmat(1,1) = tmat(1,n)-ev
    do i = 1,n-1
       cmat(2,i) = tmat(2,n-i)/cmat(1,i)
       cmat(1,i+1) = (tmat(1,n-i)-ev) - cmat(2,i)*tmat(2,n-i)
    enddo
    
    return
  end subroutine reverse_cholesky
  
  subroutine reverse_u_inv(n,cmat,x)
    implicit none
    integer, intent(in):: n
    real(kind=dp),intent(in):: cmat(2,n)
    real(kind=dp),intent(out):: x(n)
    !-------------------------------------
    integer:: i
    do i = 2,n
       x(i) = -cmat(2,n+1-i)*x(i-1)
    enddo
    return
  end subroutine reverse_u_inv
  
  subroutine reverse_l_inv(n,cmat,x)
    implicit none
    integer, intent(in):: n
    real(kind=dp),intent(in):: cmat(2,n)
    real(kind=dp),intent(out):: x(n)
    !-------------------------------------
    integer:: i
    do i = n-1, 1, -1
       x(i) = x(i) - cmat(2,n-i)*x(i+1)
    enddo
    return
  end subroutine reverse_l_inv  
  
end module tmatrix_m
