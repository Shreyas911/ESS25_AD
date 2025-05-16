      PROGRAM FWD_DRIVER

      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision J
      double precision xxs(n)

C -   initialize control parameter to zero
      DO I = 1, N
        XXS(I) = 0.0D0
      END DO

      J = 0.0d0

      call budyko_sellers( xxs, J )

      print *, 'value of J ', J
C     open a file to save J
      open(unit=11, file='J_forward.txt')      
      write(unit=11,fmt=*) J
      close(unit=11)
      
      END
