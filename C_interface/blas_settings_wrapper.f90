subroutine set_wrap_zdotc_wrapper(wrap_zdotc) bind(c, name='set_wrap_zdotc')
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  integer (c_int), intent(in), value :: wrap_zdotc
  interface
     subroutine set_wrap_zdotc(wrap_zdotc)
       logical, intent(in) :: wrap_zdotc
     end subroutine set_wrap_zdotc
  end interface
  
  logical :: setting
  setting = (wrap_zdotc /= 0)
  call set_wrap_zdotc (setting)

  return
end subroutine set_wrap_zdotc_wrapper

subroutine set_use_blas_wrapper(use_blas) bind(c, name='set_use_blas')
  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  integer (c_int), intent(in), value :: use_blas
  interface
     subroutine set_use_blas(use_blas)
       logical, intent(in) :: use_blas
     end subroutine set_use_blas
  end interface
  logical :: setting

  setting = (use_blas /= 0)
  call set_use_blas (setting)

  return
end subroutine set_use_blas_wrapper
