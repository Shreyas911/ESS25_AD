      PROGRAM AD_DRIVER

      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision J
      double precision J_AD
      double precision XXS(N)
      double precision XXS_AD(N)

      DO I = 1, N
C       Initialize the control vector to zero
        XXS(I) = 0.0D0
C       Initialize the gradient vector to zero 
        XXS_AD(I) = 0.0D0
      END DO

      J = 0.

C     Compute gradient of J with respect to all element of XXS
C     simultaneously

C     initialize with J_AD = 1
      J_AD = 1.
      call budyko_sellers_ad( XXS, XXS_AD, J, J_AD )
                    
      DO I = 1, N
          print *, 'gradient of J w.r.t. XXS ', I, XXS_AD(i) 
      END DO

C     Save the gradient vector to disk 
      open(unit=10, file='dJdX_from_adjoint.txt')
      do I = 1, N
         write(unit=10,fmt='(F24.17,A)') XXS_AD(I)
      end do
      close(unit=10)
 
      END
