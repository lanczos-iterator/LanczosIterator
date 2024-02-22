! (c) F. D. M. Haldane, Princeton University, haldane@princeton.edu
! 2014-07-31   For use with lanczos_iterator subprogram,
!-------------------------------------------------------------------------
module vector_m
  public

  logical :: use_blas = .true. , wrap_zdotc= .false.

  integer, parameter :: dp = kind(0d0)  ! kind number of double precision real numbers 
  private :: dp

  type vector_t
     integer :: n
     logical :: is_real
     real(kind=dp), pointer :: r(:) => null()
     complex(kind=dp), pointer :: c(:) => null()
  end type vector_t

  
  interface   ! to BLAS  (Basic Linear Algebra Subroutines, provided with LAPACK)
     function ddot(n,dx,incx,dy,incy) result(res)
       integer, intent(in) :: n, incx, incy
       double precision, intent(in) :: dx(*),dy(*)
       double precision :: res
     end function ddot
     function  zdotc(n,zx,incx,zy,incy) result(res)
       integer, intent(in) :: n, incx, incy
       double complex, intent(in) :: zx(*),zy(*)
       double complex :: res
     end function zdotc
     subroutine zdotc_wrapper(res,n,zx,incx,zy,incy)
       integer, intent(in) :: n, incx, incy
       double complex, intent(in) :: zx(*),zy(*)
       double complex, intent(out) :: res
     end subroutine zdotc_wrapper
     function dnrm2(n,dx,incx) result(res)
       integer, intent(in) :: n, incx
       double precision, intent(in) :: dx(*)
       double precision :: res
     end function dnrm2
     function dznrm2(n,zx,incx) result(res)
       integer, intent(in) :: n, incx
       double complex, intent(in) :: zx(*)
       double precision :: res
     end function dznrm2
     subroutine daxpy(n,da,dx,incx,dy,incy)
       integer, intent(in) :: n, incx, incy
       double precision, intent(in) :: da,dx(*)
       double precision, intent(inout) :: dy(*)
     end subroutine daxpy
     subroutine zaxpy(n,za,zx,incx,zy,incy)
       integer, intent(in) :: n, incx, incy
       double complex, intent(in) :: zx(*), za
       double complex, intent(inout):: zy(*)
     end subroutine zaxpy
     subroutine dscal(n,da,dx,incx)
       integer, intent(in) :: n, incx
       double precision, intent(in) :: da
       double precision, intent(inout) :: dx(*)
     end subroutine dscal
     subroutine zdscal(n,da,zx,incx)
       integer, intent(in) :: n, incx
       double precision, intent(in) :: da
       double complex, intent(inout):: zx(*)
     end subroutine zdscal
  end interface
       
  private :: ddot, zdotc, dnrm2, dznrm2, daxpy, zaxpy, dscal, zdscal

contains

  subroutine vector_create(n,vector_is_real,x)             ! allocate memory for the vector
    implicit none
    integer, intent(in) :: n
    logical, intent(in) :: vector_is_real
    type(vector_t), intent(inout) :: x
    
    x%n = n
    x%is_real = vector_is_real
    if (x%is_real) then
       allocate (x%r(n))
       nullify(x%c)
    else
       allocate (x%c(n))
       nullify(x%r)
    endif

    return
  end subroutine vector_create
  
  subroutine vector_erase(x)                                 ! deallocate the vector array
    implicit none
    type(vector_t), intent(inout) :: x
    
    if (x%is_real) then
       if(associated(x%r)) deallocate(x%r)
    else
       if(associated(x%c)) deallocate(x%c)
    endif

    return
  end subroutine vector_erase
  
  subroutine vector_random_fill(x,seed)                           ! x = random(seed)
    implicit none
    type (vector_t), intent(inout) :: x
    integer, intent(inout) :: seed(4)
    
    interface ! to LAPACK random number generator   LARNV
       subroutine dlarnv(idist,iseed,n,x)
         integer, intent(in) :: idist, n
         integer, intent(inout) :: iseed(4)
         double precision, intent(out) :: x(*)
       end subroutine dlarnv
       subroutine zlarnv(idist,iseed,n,x)
         integer, intent(in) :: idist, n
         integer, intent(inout) :: iseed(4)
         double complex, intent(out) :: x(*)
       end subroutine zlarnv
    end interface
    
    integer :: i
 
    !  ensure that seed(4) is valid for use by Lapack's random number generator
    do i = 1,4
       if(seed(i).lt.0) seed(i) = -seed(i)
       if (seed(i) >= 4096) seed(i) = mod(seed(i),4096)
    enddo
    if(mod(seed(4),2).eq.0) seed(4) = mod(seed(4)+1,4096)
    
    if (x%is_real) then
       call dlarnv(2,seed,x%n,x%r)
    else
       call zlarnv(2,seed,x%n,x%c)
    endif
    return
  end subroutine vector_random_fill
  
  subroutine vector_mask(x,n,mask,mask_val)                   ! set some vector elements to zero
    implicit none                                             ! (as determined by mask(n) and mask_val)
    integer, intent(in) :: n
    integer, intent(in) :: mask(n)
    integer, optional :: mask_val
    type (vector_t) :: x

    integer :: i

    if (x%n /= n) then
       write(6,'("VECTOR_MASK: incompatible mask:",2i16)') n, x%n
       stop
    end if
    
    if (present(mask_val)) then  !set elements to zero if mask(i) /= mask_val
       do i = 1,n
          if (mask(i) == mask_val) cycle  
          if (x%is_real) then
             x%r(i) = real(0,kind=dp)
          else
             x%c(i) = cmplx(0,kind=dp)
          endif
       enddo
    else                         !set elements to zero if mask(i) == 0
       do i = 1,n
          if(mask(i) == 0) then
             if (x%is_real) then
                x%r(i) = real(0,kind=dp)
             else
                x%c(i) = cmplx(0,kind=dp)
             endif
          endif
       enddo
    endif

    return
  end subroutine vector_mask
  
  subroutine vector_compatibility_test(x,y,name)         ! .true. only if x and y are vectors of same type
    implicit none
    type (vector_t), intent(in) :: x,y
    character (len=50), intent(in) :: name
    
    if (x%n /= y%n .or. x%is_real .neqv. y%is_real) then
       write(6,'(a,": x and y are incompatible",2l5,2i10)')  trim(name),&
            & x%is_real, y%is_real, x%n, y%n
       stop
    endif

    return
  end subroutine vector_compatibility_test
  
  subroutine vector_matvec(x,y,matvec,add)              ! if (add) y = y + M.x ,  else y = M.x
    implicit none
    external matvec
    type (vector_t), intent(in) , target :: x
    type (vector_t), intent (inout) ,target :: y
    logical, intent(in) :: add
    
    character (len=50) :: name = "VECTOR_MATVEC"
    call vector_compatibility_test(x,y,name)

    if(x%is_real) then
       call matvec(x%n,x%r,y%r,add)
    else
       call matvec(x%n,x%c,y%c,add)
    endif
    
    return
  end subroutine vector_matvec
  
  subroutine vector_orthogonalize(x,y,sp)                ! y = y - x(x,y),   sp = (x,y)
    implicit none
    type (vector_t), intent(in) :: x
    type (vector_t),  intent(inout):: y
    complex (kind=dp),intent(out), optional :: sp
    
    complex (kind=dp) :: sp0
    real (kind=dp) :: sp0_r
    character (len=50) :: name = "VECTOR_ORTHOGONALIZE"

    call vector_compatibility_test(x,y,name)
    
    if (x%is_real) then
       if (use_blas) then
          sp0_r = ddot(x%n,x%r,1,y%r,1)
          call daxpy(x%n,-sp0_r,x%r,1,y%r,1)
       else
          sp0_r  = sum(x%r * y%r)
          y%r = -sp0_r*x%r + y%r
       endif
       sp0 = cmplx(sp0_r,kind=dp)
    else
       if (use_blas) then
          if(wrap_zdotc) then
             call zdotc_wrapper(sp0,x%n,x%c,1,y%c,1)
          else
             sp0 = zdotc(x%n,x%c,1,y%c,1)
          endif
          call zaxpy(x%n,-sp0,x%c,1,y%c,1)
       else
          sp0 =  sum(conjg(x%c)*y%c)
          y%c = -sp0*x%c + y%c
       endif
    endif

    if (present(sp)) sp = sp0

    return
  end subroutine vector_orthogonalize

  subroutine vector_normalize(x,norm)                 ! x = x/|x|,  norm = |x|
    implicit none
    
    type (vector_t), intent(inout) :: x
    real(kind=dp), intent(out), optional :: norm

    real(kind=dp) :: norm1

    if (x%is_real) then
       if (use_blas) then
          norm1 = dnrm2(x%n,x%r,1)
          call dscal(x%n,1/norm1,x%r,1)
       else
          norm1 = sqrt(sum((x%r)**2))
          x%r = x%r/norm1
       endif
    else
       if(use_blas) then
          norm1 = dznrm2(x%n,x%c,1)
          call zdscal(x%n,1/norm1,x%c,1)
       else
          norm1 = sqrt(sum(real(x%c)**2 + aimag(x%c)**2))
          x%c = x%c/norm1
       endif
    end if

    if(present(norm)) norm = norm1
    return
  end subroutine vector_normalize

  subroutine vector_norm(x,norm)                   ! norm = |x|
    implicit none
    type (vector_t), intent(in) :: x
    real (kind=dp), intent(out) :: norm

    if (x%is_real) then
       if (use_blas) then
          norm = dnrm2(x%n,x%r,1)
       else
          norm = sqrt(sum(x%r**2))
       endif
    else
       if(use_blas) then
          norm = dznrm2(x%n,x%c,1)
       else
          norm = sqrt(sum(real(x%c)**2 + aimag(x%c)**2))
       endif
    endif

    return
  end subroutine vector_norm

  subroutine vector_scale_by_real(a,x)              ! x = a*x
    implicit none
    real (kind=dp), intent(in) :: a
    type (vector_t), intent(inout) :: x

    if(x%is_real) then
       if(use_blas) then
          call dscal(x%n,a,x%r,1)
       else
          x%r = a * x%r
       endif
    else
       if(use_blas) then
          call zdscal(x%n,a,x%c,1)
       else
          x%c = a * x%c
       endif
    end if
    return
  end subroutine vector_scale_by_real

  subroutine vector_scalar_product(x,y,sp)          ! sp = (x,y)
    implicit none
    type (vector_t), intent(in) :: x, y
    complex (kind=dp), intent(out) :: sp

    character (len=50) :: name = "VECTOR_SCALAR_PRODUCT"

    call vector_compatibility_test(x,y,name)

    if(x%is_real) then
       if (use_blas) then
          sp = ddot(x%n,x%r,1,y%r,1)
       else
          sp = sum(x%r * y%r )
       endif
    else
       if (use_blas) then
          if(wrap_zdotc) then
             call zdotc_wrapper(sp,x%n,x%c,1,y%c,1)
          else
             sp = zdotc(x%n,x%c,1,y%c,1)
          endif
       else
          sp = sum(conjg(x%c)*y%c)
       endif
    endif
    
    return
  end subroutine vector_scalar_product

  subroutine vector_change_sign(x)                  ! x = -x
    implicit none
    type (vector_t) :: x

    if  (x%is_real) then
       x%r = -x%r
    else
       x%c = -x%c
    endif
    return
  end subroutine vector_change_sign

  subroutine vector_xpy(x,y)                        ! y = x + y
    implicit none
    type (vector_t), intent(in) :: x
    type (vector_t), intent(inout) :: y
   
    character (len=50) :: name = "VECTOR_XPY"
    call vector_compatibility_test(x,y,name)

    if(x%is_real) then
       y%r = x%r + y%r
    else
       y%c = x%c + y%c
    endif
    
    return
  end subroutine vector_xpy
  
  subroutine vector_axpy(a,x,y)                     ! y = a*x + y
    implicit none
    complex (kind=dp), intent(in) :: a
    type (vector_t), intent (in) :: x
    type(vector_t), intent(inout) :: y
  
    character (len=50) :: name = "VECTOR_AXPY"
    call vector_compatibility_test(x,y,name)   


    if (x%is_real) then
       if(use_blas) then
          call daxpy(x%n,real(a),x%r,1,y%r,1)
       else
          y%r = a*x%r + y%r
       endif
    else
       if(use_blas) then
          call zaxpy(x%n,a,x%c,1,y%c,1)
       else
          y%c = a*x%c + y%c
       endif
    endif
    
    return
  end subroutine vector_axpy
  
  subroutine vector_initialize(x,i)               ! x = 0 (i=0), or  x = e(i) (i > 0)
    implicit none
    type (vector_t), intent(inout) :: x  
    integer, intent(in) :: i

    if (x%is_real) then
       x%r = real(0,kind=dp)
       if(i > 0 .and. i <= x%n)  x%r(i) = real(1,kind=dp)
    else
       x%c = cmplx(0,kind=dp)
       if(i > 0 .and. i <= x%n)  x%c(i) = cmplx(1,kind=dp)
    endif

    return
  end subroutine vector_initialize

  logical function vector_element_is_null(x,i)    ! test if x(i) == 0
    implicit none
    type (vector_t), intent(in) :: x
    integer, intent(in) :: i
    
    vector_element_is_null = .true.
    if (i > 0 .and. i <= x%n) then
       if (x%is_real) then
          if(x%r(i) /=  real(0,kind=dp))  vector_element_is_null = .false.
       else
          if(x%c(i) /= cmplx(0,kind=dp))  vector_element_is_null = .false.
       endif
    endif
   return
  end function vector_element_is_null

  subroutine vector_copy(x,y)                     !  y = x
    implicit none
    type (vector_t), intent(inout) :: x,y

    character (len=50) :: name = "VECTOR_COPY"
    call vector_compatibility_test(x,y,name)   
    
    if (x%is_real) then
       y%r = x%r
    else
       y%c = x%c
    endif

    return
  end subroutine vector_copy

end module vector_m

subroutine set_use_blas(use_blas_in)
  use vector_m
  implicit none
  logical, intent(in) :: use_blas_in

! allow control of usage of BLAS routines in vector_m

  use_blas = use_blas_in

  return
end subroutine set_use_blas

subroutine set_wrap_zdotc(wrap_zdotc_in)
  use vector_m
  implicit none
  logical, intent(in) :: wrap_zdotc_in

! activate zdotc_wrapper

  wrap_zdotc = wrap_zdotc_in

  return
end subroutine set_wrap_zdotc
