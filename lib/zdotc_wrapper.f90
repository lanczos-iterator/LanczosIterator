subroutine zdotc_wrapper(z1_dot_z2,n,z1,inc1,z2,inc2) 
  implicit none
  integer, intent(in) :: n, inc1, inc2
  complex (kind(1.0d0)), intent(in) :: z1(n), z2(n)
  complex (kind(1.0d0)), intent(out) :: z1_dot_z2
  interface
     subroutine zdotc(z1_dot_z2,n,z1,inc1,z2,inc2)
       implicit none
       complex (kind(1.0d0)), intent(out) :: z1_dot_z2
       integer, intent(in) :: n
       complex (kind(1.0d0)), intent(in) :: z1(n)
       integer, intent(in) :: inc1
       complex (kind(1.0d0)), intent(in) :: z2(n)
       integer, intent(in) :: inc2
     end subroutine zdotc
  end interface
  call zdotc(z1_dot_z2,n,z1,inc1,z2,inc2)
  return
end subroutine zdotc_wrapper

