      PROGRAM FWD_DRIVER

      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      DOUBLE PRECISION EPS
      PARAMETER (N=100, MAX_ITER=100000, EPS=1.e-3)

      double precision J, J_ORIG
      double precision XXS(N)

C -   initialize control parameter to zero
      DO I = 1, N
        XXS(I) = 0.0D0
      END DO

      J = 0.0d0

      call budyko_sellers( XXS, J )

      J_ORIG = J
      print *, 'value of J ', J
C     open a file to save J
      open(unit=11, file='J_forward.txt')
      write(unit=11,fmt='(F24.17,A)') J
      close(unit=11)

      DO I = 1, N
        XXS(I) = 0.0D0
      END DO

C     open a file to save gradients (dJ/dx)
      open(unit=111,file='dJdX_from_finite_differences.txt')
      J = 0.0d0

      DO I = 1, N
C         perturb one element of xxs
          XXS(I) = EPS
      
          call budyko_sellers( XXS, J )

C         reset perturbation to zero
          XXS(I) = 0.
          print *, 'values of gradient for I = ', I, ': ',
     &             (J-J_ORIG)/EPS
          write(unit=111,fmt='(F24.17,A)') (J-J_ORIG)/EPS
      END DO

      close(unit=111)

      END
