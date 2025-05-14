program driver

  !!Essentially both modules had variables with the same name
  !!so they had to be locally aliased for the code to compile and run
  !!Obviously only one set of variables is needed so the other is useless
  
  use forward_adj
  use forward_tgt, only: forward_problem_d

  use global_variables_adj 
  use global_variables_tgt, only: Md, Vd

  implicit none 
  
  real(8), parameter :: eps = 1.d-7
  real(8), dimension(0:nx) :: accuracy, M_fd
  real(8) :: V_orig
  integer :: i
  
  
  !! Forward run
  call forward_problem()
  V_orig = V
  
  !! Adjoint run
  M(:) = 0.0
  Vb = 1.0
  call forward_problem_b()

  open (unit = 1, file = "results.txt", action="write",status="replace")  
  write(1,*) "         #                Reverse                           FD",&
             "                          Tangent                     Relative accuracy"
  write(1,*) "______________________________________________________________",&
             "_______________________________________________________________________"
  
  !! Finite differences and Tangent Linear Model
  do i=0,nx
        
    !! TLM
    M(:) = 0.0
    Md(:) = 0.0
    Md(i) = 1.0
    V = 0.0
    Vd = 0.0
    call forward_problem_d()

    !! FD
    M(:) = 0.0
    M(i) = eps
    V = 0.0
    call forward_problem()
    M_fd(i) =  (V - V_orig)/eps

    if ( M_fd(i).NE. 0. ) then
        accuracy(i) = 1.d0 - Vd/M_fd(i)
    else
        accuracy(i) = 0.
    end if
    write(1,*) i, "    ", Mb(i), "    ", M_fd(i),"    ", Vd,"    ", accuracy(i)

  end do
  
  close(1)

end program driver

