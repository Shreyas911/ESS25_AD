      PROGRAM AD_DRIVER

      IMPLICIT NONE
      
      INTEGER N, I, ITER_GD, MAX_ITER_GD
      DOUBLE PRECISION ETA
      PARAMETER (ETA=0.0001, N=100, MAX_ITER_GD=300)

      double precision J
      double precision J_AD
      double precision XXS
      double precision XXS_AD

C       Initialize the control vector to zero
        XXS = 0.0D0
C       Initialize the gradient vector to zero 
        XXS_AD = 0.0D0

      J = 0.
      ITER_GD = 0

      XXS = -0.50
      call budyko_sellers(XXS, J)
      print *, "GradDes I ", ITER_GD, " Cost J: ", J, " XXS: ", XXS

      DO ITER_GD = 1, MAX_ITER_GD
C       initialize with J_AD = 1
        J_AD = 1.
        XXS_AD = 0.0D0
        call budyko_sellers_ad( XXS, XXS_AD, J, J_AD )
        
          XXS = XXS - ETA*XXS_AD

        call budyko_sellers(XXS, J)

        print *, "GradDes I ", ITER_GD, " Cost J: ", J, " XXS: ", XXS,
     &           " XXSAD: ", XXS_AD
      END DO                    

      END
