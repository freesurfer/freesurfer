cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Header: /space/repo/1/dev/dev/talairach_avi/eigen.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c $Log: eigen.f,v $
c Revision 1.1  2007/05/04 22:33:59  nicks
c new talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.3  2007/04/21  21:18:11  avi
c g77 compliant
c
c Revision 1.2  2000/12/13  18:37:59  avi
c copyright
c
      subroutine eigentst
      external transpos, matmul, eigen
      real*4 q(3,3),w(3,3),t(3,3),r(3,3)
      do 1 i=1,3
      do 1 j=1,3
    1 t(i,j)=rand(0)
      do 2 i=1,3
      do 2 j=1,3
      q(i,j)=0.
      do 2 k=1,3
    2 q(i,j)=q(i,j)+t(i,k)*t(j,k)
      do 4 i=1,3
      do 4 j=1,3
    4 r(i,j)=q(i,j)
      call eigen(q,w,3)
      do 3 i=1,3
    3 write(*,101)(r(i,k),k=1,3),(q(i,k),k=1,3),(w(i,k),k=1,3)
  101 format(3(3f8.4,5x))
      call transpos(w,t,3)
      call matmul(w,q,r,3)
      call matmul(r,t,q,3)
      write(*,"()")
      do 5 i=1,3
    5 write(*,101)(q(i,k),k=1,3)
      return
      end
      SUBROUTINE EIGEN(A,P,N)
      REAL*4 A(N,N),P(N,N)
      LOGICAL*4 LIND
      DATA RANGE/1.E-6/
      DO 40 L=1,N
   40 P(L,L)=1.
      DO 41 L=1,N-1
      DO 41 M=L+1,N
      P(L,M)=0.
   41 P(M,L)=0.
      THR=0.
      DO 30 L=1,N-1
      DO 30 M=L+1,N
   30 THR=THR+A(L,M)**2
      THR=SQRT(2.*THR)
   45 THR=THR/FLOAT(N)
   46 LIND=.FALSE.
      DO 50 L=1,N-1
      DO 55 M=L+1,N
      IF(ABS(A(L,M)).LE.THR)GOTO 55
      LIND=.TRUE.
      ALAMDA=-A(L,M)
      AMU=.5*(A(L,L)-A(M,M))
      OMEGA=SIGN(1.,AMU)*ALAMDA/SQRT(ALAMDA**2+AMU**2)
      ST=OMEGA/SQRT(2.*(1.+SQRT(1.-OMEGA**2)))
      ST2=ST**2
      CT2=1.-ST2
      CT=SQRT(CT2)
      STCT=ST*CT
      DO 125 I=1,N
      IF(I.EQ.L.OR.I.EQ.M)GOTO 115
      X=CT*A(I,L)-ST*A(I,M)
      A(I,M)=ST*A(I,L)+CT*A(I,M)
      A(M,I)=A(I,M)
      A(I,L)=X
      A(L,I)=X
  115 X=CT*P(I,L)-ST*P(I,M)
      P(I,M)=ST*P(I,L)+CT*P(I,M)
      P(I,L)=X
  125 CONTINUE
      X=2.*A(L,M)*STCT
      Y=A(L,L)*CT2+A(M,M)*ST2-X
      X=A(L,L)*ST2+A(M,M)*CT2+X
      A(L,M)=(A(L,L)-A(M,M))*STCT+A(L,M)*(CT2-ST2)
      A(M,L)=A(L,M)
      A(L,L)=Y
      A(M,M)=X
   55 CONTINUE
   50 CONTINUE
      IF(LIND)GOTO 46
      DIAG=ABS(A(N,N))
      OFFD=0.
      DO 48 L=1,N-1
      DIAG=DIAG+ABS(A(L,L))
      DO 48 M=L+1,N
   48 OFFD=OFFD+ABS(A(L,M))
      IF(OFFD/(.5*FLOAT(N-1)).GT.RANGE*DIAG)GOTO 45
      DO 185 L=1,N-1
      DO 185 M=L+1,N
      IF(A(L,L).GE.A(M,M))GOTO 185
      X=A(L,L)
      A(L,L)=A(M,M)
      A(M,M)=X
      DO 180 I=1,N
      X=P(I,L)
      P(I,L)=P(I,M)
  180 P(I,M)=X
  185 CONTINUE
      RETURN
      END
