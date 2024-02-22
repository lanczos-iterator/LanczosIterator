  implicit none
  integer, parameter :: dp = kind(1.0d0)
  complex(kind=dp), allocatable :: c1(:),c2(:)
  real(kind=dp), allocatable  :: a(:), b(:)
  integer :: n, test
  complex(dp) :: z
  
  interface
     function zdotc(n,z1,inc1,z2,inc2) result(z)
       implicit none
       complex(kind(1.0d0)) :: z
       integer, intent(in) :: n
       complex(kind(1.0d0)), intent(in) :: z1(n)
       integer, intent(in) :: inc1
       complex(kind(1.0d0)), intent(in) :: z2(n)
       integer, intent(in) :: inc2
     end function zdotc
  end interface

  interface
     subroutine zdotc_wrapper(z,nz,az,n1z,bz,n2z)
       implicit none
       complex(kind(1.0d0)), intent(out) :: z
       integer, intent(in) :: nz, n1z, n2z
       complex(kind(1.0d0)), intent(in) :: az(nz), bz(nz)
     end subroutine zdotc_wrapper
  end interface
  
  n = 2000

  allocate (a(n), b(n))

  call  random_number(a)
  a = 2*a -1_dp
  call  random_number(b)
  b = 2*b -1_dp
  c1 = cmplx(a,b,kind=dp)

  call  random_number(a)
  a = 2*a -1_dp
  call  random_number(b)
  b = 2*b -1_dp
  c2 = cmplx(a,b,kind=dp)


  
  print '("zdotc in your LAPACK/BLAS may need a wrapper")'
  print '("enter 0 to test with no wrapper, 1 for a wrapper")'
  read *, test
  if (test == 0) then
     print '("*** Testing with unwrapped zdotc:",/1x,&
       /" If this fails (incorrect answer or crash) you will need to use the"&
       /" zdotc wrapper when using the BLAS library linked to this test program",/1x)'
  else
     print '("*** Testing with wrapped zdotc:",/1x,&
       /" If this gives the correct result, you can use the zdotc wrapper"&
       /" when using the BLAS library linked to this test program",/1x)'
  endif


  print '("correct result is:")'
  


  print '(2e25.14)', dot_product(c1,c2)
  print '(/"zdotc reports: ")'
  

  if (test == 0) then
     z = zdotc(n,c1,1,c2,1)
  else
     call zdotc_wrapper(z,n,c1,1,c2,1)
  endif
  write(6,'(2e25.14)') z

  stop
end program
