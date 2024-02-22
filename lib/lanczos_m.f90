
! (c) F. D. M. Haldane, Princeton University, haldane@princeton.edu
! 2014-07-31   For use with lanczos_iterator subprogram,
!-------------------------------------------------------------------------
module lanczos_m
  use vector_m     ! provides all linear algebra done in Hilbert space underlying matrix. 
 
  integer, parameter, private :: dp = kind(0d0)  ! double precision
  real  (kind=dp), parameter, private :: zero = 0_dp

contains
  
  subroutine lanczos_iterate(nstep, d, e, vec1_p, vec2_p, matvec, shift, nevec, evec, eval, gap) 
    implicit none
    integer  :: nstep,  nevec
    type (vector_t), pointer :: vec1_p, vec2_p  
    real (kind=dp) ::  d, e, shift, eval(nevec), gap
    external matvec
    type (vector_t) :: evec(nevec)

    intent(in) :: shift, nevec, evec, eval, gap
    intent(out) :: d, e
    intent(inout) :: nstep
    !--------------------------------------
    ! LANCZOS iteration step,  
    !
    ! the Matrix-vector multiplication is implemented by
    !
    !    call maxpy(n,x,y,add....))
    !    <i|y> = <i|A|x> + <i|y>  i = 1,..,n if (add)
    !    <i|y> = <i|A|x>          i = 1,..,n if (.not.add)
    !
    !  A is a real Hermitian (symmetric) matrix
    !-----------------------------------------------------------
    ! on initial call with NSTEP=0, VEC(1),...VEC(N) 
    ! must contain a normalized starting vector  |1>
    ! LVEC must be at least 2*N.
    ! 
    ! on output,  NSTEP has increased by 1
    ! vec1_p points to  |NSTEP>, 
    ! vec2_p points to  |NSTEP+1>, unless E  = 0
    ! D = TMAT(1,NSTEP)
    ! E = TMAT(2,NSTEP)
    !    
    !  A|1> = TMAT(1,1)|1> + TMAT(2,1)|2>
    !  A|N> = TMAT(2,N-1)|N-1> + TMAT(1,N)|N> + TMAT(2,N)|N+1>,   N > 1.
    !
    ! <N|N> = 1, <N|N+1> = 0,  
    !
    !---------------------------------------      
    real (kind=dp), save ::  eprev
    integer, save :: nprev
    complex (kind=dp) :: diag,diag1  
    logical :: add
    type (vector_t), pointer :: temp_p   
    
    nullify(temp_p)
    
    if(nstep.eq.0) then
       nstep = 1
       add = .false.
    else
       if(nstep.ne.nprev) then
          write(6,30) nstep,nprev
30        format('LANCZOS_ITERATE: on entry,  NSTEP = ',i6, &
               & ' must match previous value =',i6) 
          stop
       endif
       nstep = nstep+1
       if(eprev.eq.zero) then
          write(6,20) nstep
20        format('LANCZOS_ITERATE: error at NSTEP=',i8,': EPREV= 0')
          stop
       endif
       
       ! exchange the targets of pointers vec1_p and vec2_p
       temp_p => vec1_p
       nullify(vec1_p)
       vec1_p => vec2_p
       vec2_p => temp_p
       nullify(temp_p)
       add = .true.
       call vector_scale_by_real(-eprev,vec2_p)
    endif
    
    call lanczos_maxpy(vec1_p,vec2_p,add, matvec,shift,nevec,evec,eval,gap)
    
   call vector_orthogonalize(vec1_p,vec2_p,diag)
    call vector_norm(vec2_p,e)
    if(e.lt.abs(diag)) then
       call vector_orthogonalize(vec1_p,vec2_p,diag1)
       call vector_norm(vec2_p,e)      
       diag = diag + diag1
    endif
    
    d = real(diag)
    if(e.ne.zero)  call vector_scale_by_real(1/e,vec2_p)
    nprev = nstep
    eprev = e
    return
  end subroutine lanczos_iterate
  
   
  subroutine lanczos_maxpy(x,y,add,matvec,shift,nevec,evec,eval,gap)
    implicit none
    integer nevec
    type (vector_t) :: x, y, evec(nevec)
    logical add
    external matvec
    real (kind=dp) :: shift,eval(nevec),gap
    !---------------------------------------
    ! for use with lanczos_iterate.
    ! if |y> = A|x>  is the basic matrix vector multiplication
    ! provided by MATVEC(n,x,y,add)
    ! maxpy supplements it so
    ! |y> = A|x> - shift|x> + sum_ij |evec(i)>M(i,j)<evec(j)|x>
    ! where |evex(i)>, i = 1,nevec are given, and M(i,J) is real hermitian
    
    !--------------------------------------
    complex  (kind=dp) :: a
    integer i
    
    if(x%is_real) then
       call matvec(x%n,x%r,y%r,add)
    else
       call matvec(x%n,x%c,y%c,add)
    endif

    if(shift.ne.zero) call vector_axpy(-cmplx(shift,kind=dp),x,y)
    
    do i = 1,nevec
       if(gap.eq.zero.and.i.eq.nevec) exit 
       call vector_scalar_product(evec(i),x,a)
       a = a*(gap + eval(nevec)-eval(i))
       call vector_axpy(a,evec(i),y)
    enddo
      
    return
  end subroutine lanczos_maxpy

  subroutine residual(evec,residual_vec,matvec,eval,var,norm) 
    implicit none
    type (vector_t), intent (inout)  :: evec
    type (vector_t), intent (inout)  :: residual_vec
    real (kind=dp), intent (out) ::  eval,var
    real (kind=dp), intent (out), optional :: norm
    external matvec
    complex (kind=dp) :: fac

    if (present(norm)) then !check evec normalization, leave it unchanged
       call vector_norm(evec, norm)
    else
       call vector_normalize(evec)
    endif
    call vector_matvec(evec,residual_vec,matvec,.false.)
    call vector_orthogonalize(evec,residual_vec,fac)
    eval = real(fac)
    ! one more time, for accurate orthogonalization
    call vector_orthogonalize(evec,residual_vec,fac)
    eval = eval + real(fac)
    call vector_norm(residual_vec,var)
    return
  end subroutine residual
  
  subroutine polish_subtraction(evec,h_evec,vec,h_vec)
    implicit none
    type(vector_t), intent(inout):: evec
    type(vector_t), intent(in):: vec,h_evec,h_vec 
    ! (H-eval)|evec>  = |h_evec>
    ! (H-eval)|vec>   = |h_vec>      
    ! find fac  so (H-eval)(|evec> - fac*|vec> has smallest residual
    !       
    real (kind=dp) ::norm
    complex (kind=dp) ::fac
    call vector_norm(h_vec,norm)
    call vector_scalar_product(h_vec,h_evec,fac)
    ! norm of h_evec should be close to abs(fac)/norm
    fac = fac/(norm**2)
    call vector_axpy(-fac,vec,evec)
    
    return
  end subroutine polish_subtraction

  subroutine matrix_properties(n,vector_is_real,matvec,seed,property,info)
    implicit none
    integer, intent(in)::  n
    logical, intent(in) :: vector_is_real
    integer, intent(inout) :: seed(4)
    logical, intent(in) :: info
    real (kind=dp), intent(out) :: property(3)
    external matvec
    !--------------------------------------------
    ! measure three characteristic properties of matvec using
    ! random state generated by seed((4)  
    ! (valid  seeds have seed(i) in range [0,4095], with iseed(4) odd)
    !
    !  property(1) = average diagonal element  d = <A>
    !  property(2) = average variance  e =  sqrt(<(A-<A>)^2>)
    !  property(3) = random error err:
    !
    !       A.x =  (A.x)_{exact}  + err |x| r_x,
    !       
    !    r_x is a random unit vector
    !------------------------------------------------   
      
    real (kind=dp) :: d,e,err,norm,eps,ten=10_dp,norm1
    type (vector_t) :: x,y,z,w
    complex (kind=dp) :: dg,myx1,myx2
    integer i

    if (info) then
       if (vector_is_real) then
          write(6,'(" Real Symmetric eigenproblem")')
       else
          write(6,'(" Complex Hermitian eigenproblem")')
       endif
    end if
    
    ! machine precision
    eps = epsilon(d)
    
    call vector_create(n,vector_is_real,x)
    call vector_create(n,vector_is_real,y)
    call vector_create(n,vector_is_real,z)
    call vector_create(n,vector_is_real,w)

    call vector_random_fill(x,seed)         ! x = random(1)
    call vector_normalize(x)                ! |x| = 1
    call vector_matvec(x,z,matvec,.false.)  ! z = M*x
    call vector_norm(z,norm1)               ! norm1 = |z|, z = z/|z|

    ! check that y=z+x where y and z result from calls to matvec(n,x,y,.true.) 
    ! and matvec*n,x,z,.false.)
 
    call vector_copy(z,y)                  ! y = z
    call vector_change_sign(y)             ! y = -y
    call vector_matvec(x,y,matvec,.true.)  ! y = y + M*x

    call vector_norm(y,norm)                !norm = |y|
    if(norm.gt.sqrt(ten*n)*eps*norm1) then
       write(6,10) norm1,norm
10     format('MATRIX_PROPERTIES(zlanczos_m) ERROR:',&
            & /'matvec(n,x,y,.true.) and matvec(n,x,z,.false.) ', &
            & /'produce inconsistent y,z:',&
              & /' for |x| = 1, |y| and |z| = |y-y| were',2d12.3)
       stop
    endif

    call vector_copy(z,y)                   ! y = z 
    call vector_orthogonalize(x,y,dg)       ! y = y - x(x,y)

    d = real(dg)
    call vector_norm(y,e)                   ! e = |y|

    call vector_random_fill(y,seed)         ! y = random(2)
    call vector_normalize(y)                ! |y| = 1

    ! hermiticity test:
    call vector_matvec(y,w,matvec,.false.)  ! w = M*y
    call vector_norm(w,norm)                ! norm = |w|
    call vector_scalar_product(y,z,myx1)    ! myx1 =  (y,z)
    call vector_scalar_product(w,x,myx2)    ! my2  = (w,x)

    if(abs(myx1-myx2).gt.sqrt(norm*norm1*ten*n)*eps) then
       write(6,20) myx1,myx2,myx1-myx2
20     format('MATRIX_PROPERTIES (zlanczos_m) WARNING:',&
            &/'<x|H|y>  =  ',2f25.16, &
            &/'<y|H|x>^* = ',2f25.16, &
            &/'difference = ',2d12.3, &
            &/' user-supplied MATVEC may be non-Hermitian') 
    endif
    
    call vector_xpy(y,x)                  ! x =  x + y
    call vector_matvec(y,z,matvec,.true.) ! z = z + M*y
    call vector_change_sign(z)            ! z = -z
    call vector_matvec(x,z,matvec,.true.) ! z = z + M*x
    call vector_norm(z,err)               ! err = |z|
    property(1) = d
    property(2) = e
    property(3) = err
    
    call vector_erase(x) ;    call vector_erase(y) 
    call vector_erase(z) ;    call vector_erase(w) 

    if (info) write(6,'(" matrix_properties: resolution limit =",1es15.3)') err
   
    return
  end subroutine matrix_properties
end module lanczos_m

