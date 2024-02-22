subroutine matrix_subspace_detector(n, vector_is_real, matvec, n_subspace, min_dimension, max_dimension, split,max_allowed)
  use vector_m
  use lanczos_m
  implicit none
  integer, intent(in) ::  n
  logical, intent(in) :: vector_is_real
  integer, intent(inout) :: n_subspace
  integer, intent(out) :: split(n), min_dimension, max_dimension
  integer, intent(in), optional :: max_allowed

  external :: matvec

  !--------------------------------------------------------------------------------------------
  ! provides information about any trivial invariant (disconnected)  subspaces in the MATVEC basis. 
  ! detailed output is placed in the output structure subspace_structure.
  !
  ! on input: 
  !   n_subspace = -k < 0  outputs (prints) more detailed data on subspaces 1 : min(k, n_subspace) 
  !   max_allowed (optional)  : terminates subspace search if more than max_allowed subspaces found
  ! 
  ! on output:
  !   n_subspace          : number of disconnected subspac   (set to 0 is max_allowed is exceeded)
  !   min_dimension       : smallest subspace dimension
  !   max_dimension       : largest subspace dimension
  !   split(n)            : index in [1: n_subspace] of each state in the basis
  !    
  !------------------------------------------------------------------------------------------
  type (vector_t), target :: vec1,vec2
  type (vector_t), pointer :: vec_in, vec_out, vec_tmp
  integer  :: subspace_radius, subspace_dimension
  integer, pointer :: temp1_p(:), temp2_p(:), temp3_p(:)
  integer :: top, pos(2), count, oldcount, p, i, k, max_nsubspace, report
  logical :: found_top

  report = 0
  if(n_subspace < 0) report = - n_subspace
  
  if (n <= 0) then
     write(6,'(" invalid input to MATRIX_SUBSPACE: n = ",i20)') n
     stop
  endif
    
  
  call vector_create(n,vector_is_real,vec1) ; call vector_create(n,vector_is_real,vec2)
  vec_in => vec1 ; vec_out => vec2
  split = 0
  top = 1
  n_subspace = 1
  subspace_radius= 0
  subspace_dimension  = 1
  call vector_initialize(vec_in,top)
  split(top) = n_subspace
  count =  1
  
  do while(count.lt.n)
     oldcount = count
     top = 0   
     found_top = .false.
     call vector_matvec(vec_in,vec_out,matvec,.false.)
     do i = 1,n
        if(split(i) == 0) then   
           if(vector_element_is_null(vec_out,i))  then
              if (.not. found_top) then
                 found_top= .true.
                 top = i
              endif
           else
              count = count + 1
              split(i) = n_subspace
              subspace_dimension = subspace_dimension + 1
           endif
        endif
     enddo
     ! if found_top == .false. there is nothing more to do. 
     ! if count == oldcount, the subspace is complete.  
     if(count > oldcount) then
        subspace_radius = subspace_radius + 1
     else if(count == oldcount) then
        if ( report > 0 .and. n_subspace <= report ) write(6,&
             &'(" n_subspace =",i5," dim=",i10," total count =",i12," radius =",i3)') n_subspace, &
             & subspace_dimension, count,subspace_radius
        if (n_subspace == report) write(6,'("remaining subspaces have total dimension: ",i12)') n - count
        if (.not.found_top) exit  ! hilbert space is exhausted
        if (n_subspace  == 1) then
           min_dimension = subspace_dimension
           max_dimension = subspace_dimension
        endif
        if (subspace_dimension > max_dimension) max_dimension = subspace_dimension
        if (subspace_dimension < min_dimension) min_dimension = subspace_dimension
        n_subspace = n_subspace + 1
        if (present(max_allowed) .and. n_subspace > max_allowed) then
           write(6,'("MATRIX_SUBSPACE_DETECTOR: more than max_allowed =",i6," subspaces detected", &
           & /"MATRIX_SUBSPACE_DETECTOR: terminating subspace search")') max_allowed
           n_subspace = 0
           return
        endif
        subspace_radius = 0
        subspace_dimension = 1
        call vector_initialize(vec_out,top)
        split(top) = n_subspace
        count = count + 1
     endif
     vec_tmp => vec_in   ! swap vec_out and vec_in
     vec_in => vec_out
     vec_out => vec_tmp
  enddo
    
  call vector_erase(vec1); call vector_erase(vec2)
      
  return
end subroutine matrix_subspace_detector

  
  
