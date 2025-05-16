      PROGRAM drive_f
      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision j
      double precision xxs(n)

C -   initialize control parameter to zero
      DO I = 1, N
        XXS(I) = 0.0D0
      END DO

      j = 0.0d0

      call f( xxs, j)

      print *, 'value of J ', j
C     open a file to save J
      open(unit=11, file='J_forward.txt')      
      write(unit=11,fmt=*) j
      close(unit=11)
      
      END
