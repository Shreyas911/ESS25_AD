module global_variables

implicit none
save

real(8), parameter :: rho = 920.0
real(8), parameter :: g = 9.2
real(8), parameter :: n = 3
real(8), parameter :: A = 1e-16
real(8), parameter :: C = 2*A/(n+2)*(rho*g)**n*(1.e3)**n
real(8), parameter :: bx = -0.0001
real(8), parameter :: M0 = 0.004
real(8), parameter :: M1 = 0.0002

real(8), parameter :: xend = 30
real(8), parameter :: dx = 1
integer, parameter :: nx = int(xend/dx)

real(8), parameter :: tend = 5000
real(8), parameter :: dt = 1/12.0
integer, parameter :: nt = int(tend/dt)

real(8), dimension(0:nx) :: zs_old, zs
real(8), dimension(0:nx) :: H_old, H
real(8), dimension(0:nx) :: xarr
real(8), dimension(0:nx) :: M 
real(8), dimension(0:nx) :: b 
real(8), dimension(0:nx-1) :: D, phi

real(8) :: V

end module global_variables
