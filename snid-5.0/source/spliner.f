************************************************************************
**                                                                    **
**                 Spline routines lifted from Vista                  **
**                                                                    **
**                                                                    **
**                  Copyright (C) 1999 John L. Tonry                  **
**                                                                    **
************************************************************************

c Included in this file:
c
c subroutine spline
c subroutine splineval



************************************************************************
* subroutine spline -- Generate a natural cubic spline
************************************************************************

      SUBROUTINE SPLINE(X,Y,NP,ERR)

C       This routine will generate a natural cubic spline on the
C       set of NP 'knot points' specified by positions {X} and values
C       {Y}.  If the knot boundries are exceeded, the initial
C       and final slope values are used to extrapolate.

C       Input:          X       An ordered list of the X coordinates
C                               of the knot points.
C                       Y       The function values at the knot points.
C                       NP      The number of knot points.

C       Output:         ERR     Set TRUE if an error accurs.

C       Author: Tod R. Lauer                    6/14/83

      PARAMETER (MAXKNOT=150)   ! Maximum number of knotpoints
      DIMENSION X(NP), Y(NP), XPOS(MAXKNOT)
      DIMENSION U(MAXKNOT), Y2(MAXKNOT), A(MAXKNOT), B(MAXKNOT)
      DIMENSION C(MAXKNOT), F(MAXKNOT), H(MAXKNOT), R(MAXKNOT)
      DIMENSION YT(MAXKNOT)
      REAL L(MAXKNOT)
      LOGICAL ERR

      COMMON /SPLINEV/ Y2, H, XPOS

C     Check parameters

      ERR     =.FALSE.
      IF (NP .GT. MAXKNOT) THEN
         PRINT *,'The number of spline points must be <=',MAXKNOT
         ERR     =.TRUE.
         RETURN
      END IF

C     Generate spline equations

      N       =NP-1
      NPT     =NP
      DO 2757 I=1, NP
         XPOS(I) =X(I)
 2757 CONTINUE

      DO 2758 I=1, N
         H(I)    =X(I+1)-X(I)
         IF (H(I) .LE. 0.0) THEN
            PRINT *,'Spline X value is out of order...'
            ERR=    .TRUE.
            do 10 k = 1,np
               write(6,'(i5,2f10.3)') k, x(k), y(k)
 10         continue



            RETURN
         END IF

         F(I)    =Y(I+1)-Y(I)
         R(I)    =F(I)/H(I)
 2758 CONTINUE

      NA      =N-1
      DO 2759 I=1, NA
         B(I)=6.0*(R(I+1)-R(I))
         C(I)=H(I+1)
         A(I)=2.0*(H(I+1)+H(I))
 2759 CONTINUE

C     Solve for second derivatives

      Y2(1)   =0.0              ! Natural cubic spline has second
      Y2(NP)  =0.0              ! derivatives of zero at the end points
      U(1)    =A(1)
      DO 2760 I=2, NA
         L(I)    =C(I-1)/U(I-1)
         U(I)    =A(I)-L(I)*C(I-1)
 2760 CONTINUE

      YT(1)   =B(1)
      DO 2761 I=2, NA
         YT(I)   =B(I)-L(I)*YT(I-1)
 2761 CONTINUE

      Y2(NA+1)=YT(NA)/U(NA)
      NA1     =NA-1
      DO 2762 I=1, NA1
         FAC             =YT(NA-I)-C(NA-I)*Y2(NA+2-I)
         Y2(NA+1-I)      =FAC/U(NA-I)
 2762 CONTINUE

      DO 2763 I=1, NP
         Y2(I)   =Y2(I)/6.0
 2763 CONTINUE

      RETURN
      END



************************************************************************
* subroutine splineval -- Find the value of the spline at XP
************************************************************************

      SUBROUTINE SPLINEVAL(XP,YP,X,Y,NPT)

      PARAMETER (MAXKNOT=150)   ! Maximum number of knotpoints
      DIMENSION X(NPT), Y(NPT), XPOS(MAXKNOT)
      DIMENSION Y2(MAXKNOT), H(MAXKNOT)

      COMMON /SPLINEV/ Y2, H, XPOS


C     Find the value of the spline at XP

C     Locate XP

      LOC     =0
      DO 2764 I=1, NPT
         IF (XP .GT. XPOS(I)) LOC=LOC+1
 2764 CONTINUE

      IF (LOC .GT. 0 .AND. LOC .LT. NPT) THEN
         SP1     =Y2(LOC)/H(LOC)
         SP2     =Y2(LOC+1)/H(LOC)
         SP3     =Y(LOC+1)/H(LOC)-Y2(LOC+1)*H(LOC)
         SP4     =Y(LOC)/H(LOC)-Y2(LOC)*H(LOC)
         DXP     =X(LOC+1)-XP
         DXM     =XP-XPOS(LOC)
         YP      =SP1*DXP**3+SP2*DXM**3
         YP      =YP+SP4*DXP+SP3*DXM

      ELSE IF (LOC .EQ. 0) THEN

C     If values are needed beyond the knot points, extrapolate with
C     the slope values on either end.

         DERIV   =-3.0*Y2(1)/H(1)*(XPOS(2)-XPOS(1))**2.0
         DERIV   =DERIV+(Y(2)-Y(1))/H(1)+(Y2(1)-Y2(2))*H(1)
         DEL     =XP-XPOS(1)
         YP      =DEL*DERIV+Y(1)

      ELSE
         DERIV   =(Y(NPT)-Y(NPT-1))/H(NPT-1)+(Y2(NPT-1)-Y2(NPT))
     1        *H(NPT-1)
         DERIV   =DERIV+3.0*Y2(NPT)/H(NPT-1)
     1        *(XPOS(NPT)-XPOS(NPT-1))**2.0
         DEL     =XP-XPOS(NPT)
         YP      =DEL*DERIV+Y(NPT)

      END IF

      RETURN

      ENTRY SPLINEDEL(XP,YP)

C     Find the derivative of the spline at XP

C     Locate XP

      LOC     =0
      DO 2765 I=1, NPT
         IF (XP .GT. XPOS(I)) LOC=LOC+1
 2765 CONTINUE

C      IF (LOC .GT. 0 .AND. LOC .LT. NP) THEN
      IF (LOC .GT. 0 .AND. LOC .LT. NPT) THEN
         SP1     =Y2(LOC)/H(LOC)
         SP2     =Y2(LOC+1)/H(LOC)
         SP3     =Y(LOC+1)/H(LOC)-Y2(LOC+1)*H(LOC)
         SP4     =Y(LOC)/H(LOC)-Y2(LOC)*H(LOC)
         DXP     =XPOS(LOC+1)-XP
         DXM     =XP-XPOS(LOC)
         YP      =-3.0*SP1*DXP**2+3.0*SP2*DXM**2-SP4+SP3

      ELSE IF (LOC .EQ. 0) THEN

C     If values are needed beyond the knot points, extrapolate with
C     the slope values on either end.

         DERIV  =-3.0*Y2(1)/H(1)*(XPOS(2)-XPOS(1))**2.0
         DERIV  =DERIV+(Y(2)-Y(1))/H(1)+(Y2(1)-Y2(2))*H(1)
         YP     =DERIV

      ELSE
         DERIV  =(Y(NPT)-Y(NPT-1))/H(NPT-1)+(Y2(NPT-1)-Y2(NPT))*H(NPT-1)
         DERIV  =DERIV+3.0*Y2(NPT)/H(NPT-1)*(XPOS(NPT)-XPOS(NPT-1))**2.0
         YP     =DERIV

      END IF

      RETURN
      END
