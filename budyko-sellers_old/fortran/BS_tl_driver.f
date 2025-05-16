      PROGRAM drive_f_tl
      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision j
      double precision j_tl
      double precision xxs(n)
      double precision xxs_tl(n)

      double precision JS(n)
      double precision JTLS(n)
      
      DO I = 1, N
        JS(I) = 0.0D0
        JTLS(I) = 0.0D0
        
        XXS(I) = 0.0D0
        xxs_tl(I) = 0.0D0
      END DO

      j = 0.0d0

      DO I = 1, N
          xxs_tl(I) = 1.
      
          call f_tl( xxs, xxs_tl, j, j_tl )
          
          xxs_tl(I) = 0.
          JS(I) = j
          JTLS(I) = j_tl
          
          print *, 'values of f and df ', j, j_tl 
      END DO

      open(unit=10,file='JS_TL.txt')
      open(unit=11,file='JTLS_TL.txt')      
      
      do i = 1, n
         write(unit=10,fmt=*) JS(i)
         write(unit=11,fmt=*) JTLS(i)         
      end do
      
      close(unit=10)
      close(unit=11)
      
      END