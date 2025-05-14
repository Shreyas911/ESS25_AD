module forward

use global_variables
implicit none

contains

subroutine forward_problem()
    
  implicit none
  
  integer :: t,i

  do i=0,nx

    xarr(i) = i*dx
    M(i) = M(i) + M0 - xarr(i)*M1
    b(i) = 1.0 + xarr(i)*bx

    zs_old(i) = b(i)
    zs(i) = b(i)

    H_old(i) = 0.0
    H(i) = 0.0

  end do

  !$AD BINOMIAL-CKP nt+1 20 1
  do t=1,nt

    do i=0,nx-1
      D(i) = C * ((H_old(i) + H_old(i+1)) / 2)**(n+2) * ((zs_old(i+1) - zs_old(i)) / dx)**(n-1)
      phi(i) = - D(i) * (zs_old(i+1) - zs_old(i)) / dx
    end do

    zs(1) = b(1)
    H(1) = 0.0

    zs(nx) = b(nx)
    H(nx) = 0.0

    do i=1,nx-1
      zs(i) = zs_old(i) + M(i)*dt - dt/dx * (phi(i) - phi(i-1))

      if (zs(i) < b(i)) then
        zs(i) = b(i)
      end if
  
      H(i) = zs(i) - b(i)
    end do

    do i=0,nx
      zs_old(i) = zs(i)
      H_old(i) = H(i)
    end do
  
  end do
  
  open (unit = 2, file = "results_forward_run.txt", action="write",status="replace")
  write(2,*) "         #                H                h                  b"
  write(2,*) "_______________________________________________________________________________"

  do i=0,nx
    V = V + H(i)*dx
    write(2,*) i, "    ", H(i), "    ", zs(i), "    ", b(i)
  end do
  
  close(2)
 
end subroutine forward_problem

end module forward
