      PROGRAM drive_f_tl
      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision j
      double precision j_tl
      double precision xxs(n)
      double precision xxs_tl(n)

      DO I = 1, N
        XXS(I) = 0.0D0
        xxs_tl(I) = 0.0D0
      END DO

C     open a file to save gradients (dJ/dx)
      open(unit=11,file='dJdx_from_tangent_linear.txt')      

      j = 0.0d0

      DO I = 1, N
C         perturb one element of xxs 
          xxs_tl(I) = 1.
      
          call f_tl( xxs, xxs_tl, j, j_tl )

C         reset perturbation to zero          
          xxs_tl(I) = 0.

          print *, 'values of gradient ', j_tl 
          write(unit=11,fmt=*) j_tl         
      END DO

      close(unit=11)
      
      END
