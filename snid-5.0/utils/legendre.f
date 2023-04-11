* All sorts of routines to fit and evaluate Legendre polynomials, and to
* convert coefficients to those of plain polynomials
*
* Copyright (C) 1980 John L. Tonry
*
*** NB Bug in polyc corrected 7/1/97.
*
      SUBROUTINE FITLPOLY(NPT,X,Y,SCALE,NCOEFF,COEFF)
*     X and Y are data arrays with NPT points to fit.
*     SCALE receives XMIN and XMAX, and returns A and B.
*     where Z = A * (X - B)
*     A polynomial of degree NCOEFF-1 is fit and the coefficients are
*     returned in COEFF.
      PARAMETER (MAXFIT=8)
      REAL*8 COV(MAXFIT*MAXFIT), VEC(MAXFIT)
      REAL*8 COEFF(1)
      REAL*8 A, B, Z, LPI, DET
      REAL*4 X(NPT), Y(NPT), SCALE(*)
      
      A = 2.0D0 / (SCALE(2) - SCALE(1))
      B = 0.5D0 * (SCALE(2) + SCALE(1))
      
      DO 5 J = 1,NCOEFF
         VEC(J) = 0
 5    CONTINUE
      DO 6 I = 1,NCOEFF*NCOEFF
         COV(I) = 0
 6    CONTINUE
      DO 10 N = 1,NPT
         Z = A * (DBLE(X(N)) - B)
         
         DO 21 J = 1,NCOEFF
            COEFF(J) = LPI(J,Z)
            VEC(J) = VEC(J) + COEFF(J) * Y(N)
            DO 20 I = 1,J
               COV((J-1)*NCOEFF+I) = COV((J-1)*NCOEFF+I) + 
     $              COEFF(I) * COEFF(J)
 20         CONTINUE
 21      CONTINUE
 10   CONTINUE
      
      DO 31 J = 2,NCOEFF
         DO 30 I = 1,J-1
            COV((I-1)*NCOEFF+J) = COV((J-1)*NCOEFF+I)
 30      CONTINUE
 31   CONTINUE
      CALL INVERT(NCOEFF,COV,COEFF,DET)
      IF(DET.EQ.0) THEN
         WRITE(6,*) 'FITLPOLY: singular covariance matrix'
         RETURN
      END IF
      
      DO 41 J = 1,NCOEFF
         COEFF(J) = 0
         DO 40 I = 1,NCOEFF
            COEFF(J) = COEFF(J) + COV((J-1)*NCOEFF+I) * VEC(I)
 40      CONTINUE
 41   CONTINUE
      SCALE(1) = A
      SCALE(2) = B
      
      RETURN
      END
      
      SUBROUTINE FITWLPOLY(NPT,X,Y,W,SCALE,NCOEFF,COEFF)
* X and Y are data arrays with NPT points to fit.
* W is a weight array, so that Sum W * (Y-Yfit)**2 is minimized
* SCALE receives XMIN and XMAX from caller, and returns A and B,
* where Z = A * (X - B) is the variable given to the polynomials.
* A polynomial of degree NCOEFF-1 is fit and the coefficients are
* returned in COEFF.
      PARAMETER (MAXFIT=8)
      REAL*8 COV(MAXFIT*MAXFIT), VEC(MAXFIT)
      REAL*8 COEFF(1)
      REAL*8 A, B, Z, LPI, DET
      REAL*4 X(NPT), Y(NPT), SCALE(2), W(NPT)
      
      A = 2.0D0 / (SCALE(2) - SCALE(1))
      B = 0.5D0 * (SCALE(2) + SCALE(1))
      
      DO 5 J = 1,NCOEFF
         VEC(J) = 0
 5    CONTINUE
      DO 6 I = 1,NCOEFF*NCOEFF
         COV(I) = 0
 6    CONTINUE
      DO 10 N = 1,NPT
         Z = A * (DBLE(X(N)) - B)
         
         DO 21 J = 1,NCOEFF
            COEFF(J) = LPI(J,Z)
            VEC(J) = VEC(J) + COEFF(J) * Y(N) * W(N)
            DO 20 I = 1,J
               COV((J-1)*NCOEFF+I) = COV((J-1)*NCOEFF+I) + 
     $              COEFF(I) * COEFF(J) * W(N)
 20         CONTINUE
 21      CONTINUE
 10   CONTINUE
      
      DO 31 J = 2,NCOEFF
         DO 30 I = 1,J-1
            COV((I-1)*NCOEFF+J) = COV((J-1)*NCOEFF+I)
 30      CONTINUE
 31   CONTINUE
      CALL INVERT(NCOEFF,COV,COEFF,DET)
      IF(DET.EQ.0) THEN
         WRITE(6,*) 'FITWLPOLY: singular covariance matrix'
         RETURN
      END IF
      
      DO 41 J = 1,NCOEFF
         COEFF(J) = 0
         DO 40 I = 1,NCOEFF
            COEFF(J) = COEFF(J) + COV((J-1)*NCOEFF+I) * VEC(I)
 40      CONTINUE
 41   CONTINUE
      SCALE(1) = A
      SCALE(2) = B
      
      RETURN
      END
      
      REAL FUNCTION LPOLY(X,SCALE,NCOEFF,COEFF)
      
      REAL*8 COEFF(1)
      REAL*8 TEMP, Z, LPI
      REAL*4 SCALE(*)
      
      Z = SCALE(1) * (X - SCALE(2))
      TEMP = 0
      DO 10 I = 1,NCOEFF
         TEMP = TEMP + COEFF(I) * LPI(I,Z)
 10   CONTINUE
      LPOLY = TEMP
      RETURN
      END

      SUBROUTINE FIT2DLPOLY(NPT,X,Y,Z,SCALE,NCX,NCY,COEFF)
* X, Y and Z are data arrays with NPT points to fit.
* SCALE receives XMIN, XMAX, YMIN and YMAX and returns AX, BX, AY, BY
* where X' = AX * (X - BX)  and Y' = AY * (Y - BY)
* A two-dimensional polynomial of degree NCX-1 in X and NCY-1 in Y
* is fit and the coefficients are returned in COEFF.
      PARAMETER (MAXFIT=8)
      REAL*8 COV(MAXFIT*MAXFIT*MAXFIT*MAXFIT), VEC(MAXFIT*MAXFIT)
      REAL*8 COEFF(1)
      REAL*8 AX, BX, AY, BY, XSC, YSC, LPI, DET
      REAL*4 X(1), Y(1), Z(1), SCALE(*)

      AX = 2.0D0 / (SCALE(2) - SCALE(1))
      BX = 0.5D0 * (SCALE(2) + SCALE(1))
      AY = 2.0D0 / (SCALE(4) - SCALE(3))
      BY = 0.5D0 * (SCALE(4) + SCALE(3))

      NFIT = NCX*NCY

      DO 5 J = 1,NFIT
5     VEC(J) = 0
      DO 6 I = 1,NFIT*NFIT
6     COV(I) = 0

      DO 10 N = 1,NPT
      XSC = AX * (X(N) - BX)
      YSC = AY * (Y(N) - BY)

      DO 20 J = 1,NFIT
      JX = (J-1)/NCY + 1
      JY = MOD(J-1,NCY) + 1
      COEFF(J) = LPI(JX,XSC) * LPI(JY,YSC)
      VEC(J) = VEC(J) + COEFF(J) * Z(N)
      DO 20 I = 1,J
20    COV(I+(J-1)*NFIT) = COV(I+(J-1)*NFIT) + COEFF(I) * COEFF(J)

10    CONTINUE

      DO 30 J = 2,NFIT
      DO 30 I = 1,J-1
30    COV(J+(I-1)*NFIT) = COV(I+(J-1)*NFIT)

      CALL INVERT(NFIT,COV,COEFF,DET)
      IF(DET.EQ.0) THEN
          WRITE(6,*) 'FIT2DLPOLY: singular covariance matrix'
          RETURN
      END IF

      DO 40 J = 1,NFIT
      COEFF(J) = 0
      DO 40 I = 1,NFIT
40    COEFF(J) = COEFF(J) + COV(I+(J-1)*NFIT) * VEC(I)

      SCALE(1) = AX
      SCALE(2) = BX
      SCALE(3) = AY
      SCALE(4) = BY

      RETURN
      END

      REAL FUNCTION LPOLY2D(X,Y,SCALE,NCX,NCY,COEFF)
* The desired fit has NCX coefficients in X and NCY in Y.
      REAL*8 COEFF(1), XSC, YSC, LPI, TEMP
      REAL*4 SCALE(*)

      NFIT = NCX*NCY

      XSC = SCALE(1) * (X - SCALE(2))
      YSC = SCALE(3) * (Y - SCALE(4))

      TEMP = 0

      DO 20 J = 1,NFIT
      JX = (J-1)/NCY + 1
      JY = MOD(J-1,NCY) + 1
20    TEMP = TEMP + COEFF(J) * LPI(JX,XSC) * LPI(JY,YSC)

      LPOLY2D = TEMP

      RETURN
      END

      REAL FUNCTION POLY2D(X,Y,NCX,NCY,COEFF)
* The desired fit has NCX coefficients in X and NCY in Y.
      REAL*8 COEFF(1), T1, T2, DX, DY

      DX = DBLE(X)
      DY = DBLE(Y)

      T2 = 0

      DO 900 J = NCX*NCY,NCY,-NCY
          T1 = 0
          DO 901 I = J,J-NCY+1,-1
              T1 = T1 * DY + COEFF(I)
901       CONTINUE
          T2 = T2 * DX + T1
900   CONTINUE

      POLY2D = T2

      RETURN
      END

      FUNCTION LPI(I,Z)
* Legendre polynomial of order i-1
      REAL*8 LPI
      INCLUDE 'legendre.inc'
      GOTO(1,2,3,4,5,6,7,8) I
1     LPI = LP0(Z)
      RETURN
2     LPI = LP1(Z)
      RETURN
3     LPI = LP2(Z)
      RETURN
4     LPI = LP3(Z)
      RETURN
5     LPI = LP4(Z)
      RETURN
6     LPI = LP5(Z)
      RETURN
7     LPI = LP6(Z)
      RETURN
8     LPI = LP7(Z)
      RETURN
      END

      FUNCTION EVALEG(ORDER,X,SCALE,C)
      REAL*8 C(*)
      REAL*4 SCALE(2)
      INTEGER*4 ORDER
* Evaluates sum ( C(i) * LPi(a*(x-b)) )
* Legendre polynomials 0 - 7
      INCLUDE 'legendre.inc'

      EVALEG = 0.0D0

      Z = SCALE(1) * (X - SCALE(2))
      NCOEFF = ORDER + 1
      IF(NCOEFF.LT.1.OR.NCOEFF.GT.8) THEN
          WRITE(6,*) ' EVALEG: wrong number of coeffs',NCOEFF
          RETURN
      END IF

      GOTO(1,2,3,4,5,6,7,8) NCOEFF
1     EVALEG = C(1)*LP0(Z)
      RETURN
2     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z)
      RETURN
3     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z) + C(3)*LP2(Z)
      RETURN
4     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z) + C(3)*LP2(Z) +
     1    C(4)*LP3(Z)
      RETURN
5     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z) + C(3)*LP2(Z) +
     1    C(4)*LP3(Z) + C(5)*LP4(Z)
      RETURN
6     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z) + C(3)*LP2(Z) +
     1    C(4)*LP3(Z) + C(5)*LP4(Z) + C(6)*LP5(Z)
      RETURN
7     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z) + C(3)*LP2(Z) +
     1    C(4)*LP3(Z) + C(5)*LP4(Z) + C(6)*LP5(Z) +
     1    C(7)*LP6(Z)
      RETURN
8     EVALEG = C(1)*LP0(Z) + C(2)*LP1(Z) + C(3)*LP2(Z) +
     1    C(4)*LP3(Z) + C(5)*LP4(Z) + C(6)*LP5(Z) +
     1    C(7)*LP6(Z) + C(8)*LP7(Z)
      RETURN
      END

      SUBROUTINE POLYC(NC,SCALE,COEFF,PCOEFF)
* Routine to convert coefficients of a 1-d Legendre polynomial to
* coefficients of a plain polynomial
* The NC coefficients are stored with the low order coefficient first
      REAL*8 COEFF(1), PCOEFF(1), SC(2), TEMP
      REAL*4 SCALE(2)
      INCLUDE 'legendre.inc'

      SC(1) = DBLE(SCALE(1))
      SC(2) = DBLE(SCALE(2))
      DO 900 K = 0,NC-1
          TEMP = 0
          DO 901 I = K,NC-1
              TEMP = TEMP + COEFF(I+1) * LCOEFF(K+1,I+1)
901       CONTINUE
          PCOEFF(K+1) = TEMP
900   CONTINUE
      DO 902 K = 0,NC-1
         IF(SC(2).NE.0) THEN
            TEMP = 0
            DO 903 I = K,NC-1
               TEMP = TEMP + PCOEFF(I+1) * SC(1)**I *
     1              (-SC(2))**(I-K) * DFLOAT(IBC(K,I))
 903        CONTINUE
            PCOEFF(K+1) = TEMP
         ELSE
            PCOEFF(K+1) = PCOEFF(K+1) * SC(1)**K
         END IF
902   CONTINUE

      RETURN
      END

      SUBROUTINE POLYCC(NX,NY,SCALE,COEFF,PCOEFF)
* Routine to convert coefficients of a 2-d Legendre polynomial to
* coefficients of a plain polynomial
* The coefficients are stored as NX polynomials of order NY-1.
* Low order coefficient comes first, and the low order coefficient
* of the x polynomial is obtained by evaluating the first y polynomial.
      REAL*8 COEFF(1), PCOEFF(1), XSC(2), YSC(2), TEMP
      REAL*4 SCALE(4)
      INCLUDE 'legendre.inc'

      XSC(1) = DBLE(SCALE(1))
      XSC(2) = DBLE(SCALE(2))
      YSC(1) = DBLE(SCALE(3))
      YSC(2) = DBLE(SCALE(4))

      DO 10 J = 0,NX-1
      DO 11 K = 0,NY-1
          TEMP = 0
          DO 12 I = K,NY-1
              TEMP = TEMP + COEFF(J*NY+I+1) * LCOEFF(K+1,I+1)
12        CONTINUE
          PCOEFF(J*NY+K+1) = TEMP
11    CONTINUE
      DO 13 K = 0,NY-1
          TEMP = 0
          DO 900 I = K,NY-1
              IF(YSC(2).NE.0) THEN
                  TEMP = TEMP + PCOEFF(J*NY+I+1) * YSC(1)**I *
     1            (-YSC(2))**(I-K) * DFLOAT(IBC(K,I))
              ELSE
                  TEMP = TEMP + PCOEFF(J*NY+I+1) * YSC(1)**I
              END IF
900           CONTINUE
          PCOEFF(J*NY+K+1) = TEMP
13    CONTINUE
10    CONTINUE

      DO 20 J = 0,NY-1
      DO 21 K = 0,NX-1
          TEMP = 0
          DO 22 I = K,NX-1
              TEMP = TEMP + PCOEFF(I*NY+J+1) * LCOEFF(K+1,I+1)
22        CONTINUE
          PCOEFF(K*NY+J+1) = TEMP
21    CONTINUE
      DO 23 K = 0,NX-1
          TEMP = 0
          DO 901 I = K,NX-1
              IF(XSC(2).NE.0) THEN
                  TEMP = TEMP + PCOEFF(I*NY+J+1) * XSC(1)**I *
     1              (-XSC(2))**(I-K) * DFLOAT(IBC(K,I))
              ELSE
                  TEMP = TEMP + PCOEFF(I*NY+J+1) * XSC(1)**I
              END IF
901       CONTINUE
          PCOEFF(K*NY+J+1) = TEMP
23    CONTINUE
20    CONTINUE

      RETURN
      END

      FUNCTION IBC(M,N)
* Returns the (m n) binomial coefficient
      IBC = 1
      DO 900 I = M+1,N
          IBC = IBC * I
900   CONTINUE
      DO 901 I = 2,N-M
          IBC = IBC / I
901   CONTINUE
      RETURN
      END
