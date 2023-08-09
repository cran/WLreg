
SUBROUTINE WRLOGISTICA(N,P,Z,WLSTATUS,BETA0,TOL,MAXITER,BETA,UBETA,VU,VBETA,IMATRIX,WTOTAL,LTOTAL,ERR,ITER)
   !USE NUMERICAL_LIBRARIES
   !USE LINEAR_OPERATORS
   IMPLICIT NONE
   INTEGER,INTENT(IN)::N,P,MAXITER,WLSTATUS(N*(N-1)/2)
   REAL(8),INTENT(IN)::Z(N,P),TOL
   INTEGER,INTENT(IN)::BETA0(P)
   REAL(8),INTENT(OUT)::BETA(P,1),UBETA(P,1),VU(P,P),VBETA(P,P),IMATRIX(P,P),WTOTAL(P,1),LTOTAL(P,1),ERR
   INTEGER,INTENT(OUT)::ITER
   
   INTEGER::I,J,INDEXIJ
   REAL(8)::B0(P,1),B1(P,1),XB(P,1),XTEMP,PIJ,QIJ,AIB(N,P),PPMX(P,P),P2PMX(P,P),PNMX(P,N),TXB(1,P),INVM(P,P),P1MX(P,1)
   INTEGER::WIJ,LIJ
   
   ERR=1.0;ITER=0
   B0(:,1)=BETA0;B1=B0
   DO WHILE (ERR>TOL .AND. ITER<MAXITER)
      ITER=ITER+1
      UBETA=0.0;VBETA=0.0;IMATRIX=0.0;WTOTAL=0.0;LTOTAL=0.0
      DO I=2,N,1
         DO J=1,I-1,1
            INDEXIJ=(I-1)*(I-2)/2+J
            IF (WLSTATUS(INDEXIJ) /= 0) THEN
               WIJ=0;LIJ=0
               IF (WLSTATUS(INDEXIJ)>0) THEN
                  WIJ=1
               ELSE 
                  LIJ=1
               END IF
               XB(:,1)=Z(I,:)-Z(J,:)
               TXB=TRANSPOSE(XB) 
               XTEMP=DEXP(SUM(XB(:,1)*B0(:,1)))
               PIJ=XTEMP/(1.0+XTEMP);QIJ=1.0-PIJ
               UBETA(:,1)=UBETA(:,1)+XB(:,1)*(REAL(WIJ,8)-PIJ)
               PPMX=MATMUL(XB,TXB)
               IMATRIX=IMATRIX+PPMX*PIJ*QIJ 
               WTOTAL(:,1)=WTOTAL(:,1)-XB(:,1)*REAL(LIJ,8)
               LTOTAL(:,1)=LTOTAL(:,1)-XB(:,1)*REAL(WIJ,8)
            END IF
         END DO
      END DO
      PPMX=2.0*IMATRIX/REAL(N,8)/REAL(N,8)
      UBETA=2.0*UBETA/REAL(N,8)/REAL(N,8)
      CALL INV1(PPMX,P,INVM)
      P1MX=MATMUL(INVM,UBETA)
      B1=B0+P1MX
      BETA=B0
      ERR=MAXVAL(DABS(B1(:,1)-B0(:,1)))
      B0=B1
   END DO

   !Variance
      AIB=0.0;B0=BETA
      DO I=1,N,1
         DO J=1,N,1
            IF (J<I) THEN
               INDEXIJ=(I-1)*(I-2)/2+J
               IF (WLSTATUS(INDEXIJ) /= 0) THEN
                  WIJ=0;LIJ=0
                  IF (WLSTATUS(INDEXIJ)>0) THEN
                     WIJ=1
                  ELSE 
                     LIJ=1
                  END IF
                  XB(:,1)=Z(I,:)-Z(J,:)
                  XTEMP=DEXP(SUM(XB(:,1)*B0(:,1)))
                  PIJ=XTEMP/(1.0+XTEMP)
                  AIB(I,:)=AIB(I,:)+XB(:,1)*(REAL(WIJ,8)-PIJ)
               END IF
            ELSE IF (J>I) THEN
               INDEXIJ=(J-1)*(J-2)/2+I
               IF (WLSTATUS(INDEXIJ) /= 0) THEN
                  WIJ=0;LIJ=0
                  IF (WLSTATUS(INDEXIJ)>0) THEN
                     WIJ=1
                  ELSE 
                     LIJ=1
                  END IF
                  XB(:,1)=Z(J,:)-Z(I,:)
                  XTEMP=DEXP(SUM(XB(:,1)*B0(:,1)))
                  PIJ=XTEMP/(1.0+XTEMP)
                  AIB(I,:)=AIB(I,:)+XB(:,1)*(REAL(WIJ,8)-PIJ)
               END IF
            END IF
         END DO
      END DO

   AIB=AIB/REAL(N,8)
   PNMX=TRANSPOSE(AIB)
   PPMX=MATMUL(PNMX,AIB)
   VU=4.0*PPMX/REAL(N,8)
   PPMX=MATMUL(INVM,VU)
   P2PMX=TRANSPOSE(INVM)
   VBETA=MATMUL(PPMX,P2PMX)
   RETURN
END SUBROUTINE WRLOGISTICA

! -------------------------------------------------------------------- 
SUBROUTINE INV1 (a,n,ai)       ! Invert matrix by Gauss method 
! --------------------------------------------------------------------
IMPLICIT NONE 
INTEGER :: n 
REAL(8) :: a(n,n)
REAL(8),INTENT(OUT) :: ai(n,n)
! - - - Local Variables - - - 
REAL(8) :: b(n,n), c, d, temp(n) 
INTEGER :: i, j, k, m, imax(1), ipvt(n) 
! - - - - - - - - - - - - - - 
b = a 
ipvt = (/ (i, i = 1, n) /) 
DO k = 1,n 
   imax = MAXLOC(ABS(b(k:n,k))) 
   m = k-1+imax(1) 
   IF (m /= k) THEN 
      ipvt( (/m,k/) ) = ipvt( (/k,m/) ) 
      b((/m,k/),:) = b((/k,m/),:) 
   END IF 
   d = 1/b(k,k) 
   temp = b(:,k) 
   DO j = 1, n 
      c = b(k,j)*d 
      b(:,j) = b(:,j)-temp*c 
      b(k,j) = c 
   END DO 
   b(:,k) = temp*(-d) 
   b(k,k) = d 
END DO 
a(:,ipvt) = b 
ai=a
END SUBROUTINE INV1

