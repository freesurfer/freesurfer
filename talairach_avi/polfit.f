c Revision 1.1  2007/05/04 22:34:03  nicks
c new talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.1  1996/04/19  17:12:51  ty7777
c Initial revision
c
c
      SUBROUTINE POLFIT(X,Y,NPTS,A,ARRAY,NTERM,CHISQR)
C     BEVINGTON p. 141
      character*256 rcsid /'$Header: /space/repo/1/dev/dev/talairach_avi/polfit.f,v 1.1 2007/05/04 22:34:03 nicks Exp $'/
      REAL*4 X(NPTS),Y(NPTS)
      REAL*4 A(NTERM),ARRAY(NTERM,NTERM)
      REAL*4 SUMX(19),SUMY(10)
      IF(NTERM.GT.10)STOP 'POLFIT NTERM TOO LARGE'
      NMAX=2*NTERM-1
      DO 13 N=1,NMAX
   13 SUMX(N)=0.
      DO 15 J=1,NTERM
   15 SUMY(J)=0.
      CHISQR=0.
      DO 50 I=1,NPTS
      XTERM=1.
      DO 44 N=1,NMAX
      SUMX(N)=SUMX(N)+XTERM
   44 XTERM=XTERM*X(I)
      YTERM=Y(I)
      DO 48 N=1,NTERM
      SUMY(N)=SUMY(N)+YTERM
   48 YTERM=YTERM*X(I)
      CHISQR=CHISQR+Y(I)**2
   50 CONTINUE
      DO 54 J=1,NTERM
      DO 54 K=1,NTERM
   54 ARRAY(J,K)=SUMX(J+K-1)
      if(.true.)then
      CALL MATINV(ARRAY,NTERM,DET)
      DO 61 L=1,NTERM
      A(L)=0.
      do 61 J=1,NTERM
   61 A(L)=A(L)+ARRAY(L,J)*SUMY(J)
      else
      DELTA=DETERM(ARRAY,NTERM)
      DO 70 L=1,NTERM
      DO 66 J=1,NTERM
      DO 65 K=1,NTERM
   65 ARRAY(J,K)=SUMX(J+K-1)
   66 ARRAY(J,L)=SUMY(J)
   70 A(L)=DETERM(ARRAY,NTERM)/DELTA
      endif
      DO 75 J=1,NTERM
      CHISQR=CHISQR-2.*A(J)*SUMY(J)
      DO 75 K=1,NTERM
   75 CHISQR=CHISQR+A(J)*A(K)*SUMX(J+K-1)
      CHISQR=CHISQR/FLOAT(NPTS-NTERM)
      RETURN
      END
      FUNCTION DETERM(ARRAY,NORDER)
      REAL*4 ARRAY(NORDER,NORDER)
      DETERM=1.
      DO 50 K=1,NORDER
      IF(ARRAY(K,K).NE.0.)GOTO 41
      DO 23 J=K,NORDER
      IF(ARRAY(J,K).NE.0.)GOTO 31
   23 CONTINUE
      DETERM=0.
      RETURN
   31 DO 34 I=K,NORDER
      SAVE=ARRAY(I,J)
      ARRAY(I,J)=ARRAY(I,K)
   34 ARRAY(I,K)=SAVE
      DETERM=-DETERM
   41 DETERM=DETERM*ARRAY(K,K)
      IF(K.GE.NORDER)GOTO 50
      DO 46 I=K+1,NORDER
      DO 46 J=K+1,NORDER
   46 ARRAY(I,J)=ARRAY(I,J)-ARRAY(I,K)*ARRAY(K,J)/ARRAY(K,K)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE POLTST
      PARAMETER (NPTS=51)
      REAL*4 X(NPTS),Y(NPTS)
      PARAMETER (NORDER=10)
      REAL*4 A(NORDER),ARRAY(NORDER,NORDER),C(NORDER)
      write(*,"(' ENTER NTERM')")
      read(*,*)NTERM
      IF(NTERM.GT.NORDER)STOP 'POLTST NTERM TOO LARGE'
      write(*,"(' ENTER NTERM COEFFICIENTS')")
      read(*,*)(C(N),N=1,NTERM)
      DO 50 I=1,NPTS
      T=-1.+2.*FLOAT(I-1)/FLOAT(NPTS-1)
      X(I)=T
      Y(I)=C(1)
      DO 2 N=2,NTERM
      Y(I)=Y(I)+C(N)*T
    2 T=T*X(I)
   50 CONTINUE
      CALL POLFIT(X,Y,NPTS,A,ARRAY,NTERM,CHISQR)
      DO 10 N=1,NTERM
   10 write(*,101)N,C(N),A(N)
  101 FORMAT(I3,2F10.4)
      write(*,102)CHISQR
  102 FORMAT(' CHISQR',F10.4)
      STOP
      END
