     
C ---------------------------------
      SUBROUTINE BUDYKO_SELLERS(XXS, J)
      
      INTEGER N, I, ITER, MAX_ITER
      PARAMETER (N=100, MAX_ITER=100000)
      
      REAL*8 X(N), XEDGES(N+1), DX
      REAL*8 LAT(N), T(N), TNEW(N), ALPHA(N)
      REAL*8 FIN(N), FOUT(N), FDIFF(N), DT, DIFF, TOL
      REAL*8 S0, Q, SIGMA, EPSILON
      REAL*8 SX(N), DTDX, DTDX_M, DTDX_P
      REAL*8 XM, XP, A1, A2, T1, T2
      
      REAL*8 XXS
      REAL*8 J

CADJ INIT tapex = 'TAF_tape'

C --- Constants
      S0 = 1366.D0
      Q = S0 / 4.D0
      SIGMA = 5.67D-8
      EPSILON = 0.63D0
      DIFF = 0.6D0 + XXS
      TOL = 1.D-9

C --- Grid setup (x = sin(lat))
      DO I = 1, N+1
         XEDGES(I) = -1.D0 + 2.D0*(I-1)/FLOAT(N)
      END DO
      DO I = 1, N
         X(I) = 0.5D0 * (XEDGES(I) + XEDGES(I+1))
      END DO
      DX = X(2) - X(1)

C --- Latitude in degrees
      DO I = 1, N
         LAT(I) = ASIN(X(I)) * 180.D0 / 3.14159265D0
      END DO

C --- Insolation
      DO I = 1, N
         SX(I) = Q*(1.D0 - 0.482D0 * X(I)**2)*1.35D0-55.D0
      END DO

C --- Initial temperature guess
      DO I = 1, N
         T(I) = 288.D0 + 60.D0 * (1.D0 - X(I)**2) - 20.D0
      END DO

C --- Time step (CFL-like)
      DT = DX**2 / (2.D0 * DIFF) * 0.5D0

C --- Albedo parameters
      T1 = 223.15D0
      T2 = 293.15D0
      A1 = 0.9D0
      A2 = 0.2D0

C --- Iterative solver to find T(x) that
C --- balances Energy In = Energy Out

      DO ITER = 1, MAX_ITER
CADJ STORE t = tapex

C --- Albedo (linear)
         DO I = 1, N
            IF (T(I) .LE. T1) THEN
               ALPHA(I) = A1
            ELSE IF (T(I) .GE. T2) THEN
               ALPHA(I) = A2
            ELSE
               ALPHA(I) = A1 + (A2 - A1)*(T(I) - T1)/(T2 - T1)
            END IF
         END DO

C --- Radiation terms
         DO I = 1, N
            FIN(I) = SX(I) * (1.D0 - ALPHA(I))
            FOUT(I) = EPSILON * SIGMA * T(I)**4
         END DO

C --- Diffusion
         DO I = 1, N
            IF (I .EQ. 1) THEN
               DTDX_M = 0.D0
            ELSE
               DTDX_M = (T(I) - T(I-1)) / DX
            END IF
            IF (I .EQ. N) THEN
               DTDX_P = 0.D0
            ELSE
               DTDX_P = (T(I+1) - T(I)) / DX
            END IF
            XM = XEDGES(I)
            XP = XEDGES(I+1)
            FDIFF(I) = DIFF * ((1.D0 - XP**2) * DTDX_P -
     &                        (1.D0 - XM**2) * DTDX_M) / DX
         END DO

C --- Update temperature
         DO I = 1, N
            TNEW(I) = T(I) + DT * (FIN(I) - FOUT(I) + FDIFF(I))
         END DO
         DO I = 1, N
            T(I) = TNEW(I)
         END DO
      
      END DO


C --- J is the Equator to Pole temperature difference
C --- Takes the average of the two bands around the equator
C --- 1/2 * ( T(50) + T(51) ) 
C --- and the average the average of the band around the north pole
C --- and south pole 
C --- 1/2 * ( T(1)  + T(N) )

      J = 0.5*(0.5*( T(50) + T(51)) - 0.5* ( T(1) + T(N) ) - 45.0)**2
C ---     Avg. T around Eq.    -   Avg T around poles
      
      END
