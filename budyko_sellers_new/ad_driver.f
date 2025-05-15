      PROGRAM drive_f_ad
      IMPLICIT NONE
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)

      double precision j
      double precision j_ad
      double precision xxs(n)
      double precision xxs_ad(n)

      DO I = 1, N
        
        XXS(I) = 0.0D0
        xxs_ad(I) = 0.0D0
      END DO

      j = 0.
      j_ad = 1.

      call f_ad( xxs, xxs_ad, j, j_ad )
                    
      DO I = 1, N
          print *, 'gradient of J w.r.t. xxs ', I, xxs_ad(i) 
      END DO

      print *, 'j_ad ', j_ad
      

      open(unit=10,file='xxs_ad.txt')
      do i = 1, n
         write(unit=10,fmt=*) xxs_ad(i)
      end do
      close(unit=10)
      
      
      
      END