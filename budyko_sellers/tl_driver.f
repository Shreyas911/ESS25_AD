      PROGRAM TL_DRIVER

      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision J
      double precision J_TL
      double precision XXS(N)
      double precision XXS_TL(N)

      DO I = 1, N
        XXS(I) = 0.0D0
        XXS_TL(I) = 0.0D0
      END DO

C     open a file to save gradients (dJ/dx)
      open(unit=11,file='dJdx_from_tangent_linear.txt')      

      J = 0.0d0

      DO I = 1, N
C         perturb one element of xxs 
          XXS_TL(I) = 1.
      
          call budyko_sellers_tl( XXS, XXS_TL, J, J_TL )

C         reset perturbation to zero          
          XXS_TL(I) = 0.

          print *, 'values of gradient for I = ', I, ': ', J_TL 
          write(unit=11,fmt=*) J_TL
      END DO

      close(unit=11)
      
      END
