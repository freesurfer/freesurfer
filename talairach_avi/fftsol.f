c$Header: /space/repo/1/dev/dev/talairach_avi/fftsol.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c$Log: fftsol.f,v $
cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  2007/03/05  00:26:13  avi
c dimension nfac(11) -> dimension nfac(13)
c nfac array index limit (maxnfac) check installed
c
c Revision 1.1  2007/03/04  23:20:21  avi
c Initial revision
c
      subroutine fft(a,b,nseg,n,nspn,isn)
      real*4 a(nseg*n*nspn),b(nseg*n*nspn)
      ntot=nseg*n*nspn
      nspan=n*nspn
      call fft2(a,b,ntot,n,nspan,isn)
      if (isn.gt.0)then
        on=1./float(n)
        do i=1,ntot
          a(i)=a(i)*on
          b(i)=b(i)*on
        enddo
      endif
      return
      end
      subroutine fft2(a,b,ntot,n,nspan,isn)
c  multivariate complex fourier transform, computed in place
c    using mixed-radix fast fourier transform algorithm.
c  by r. c. singleton, stanford research institute, sept. 1968
c  arrays a and b originally hold the real and imaginary
c    components of the data, and return the real and
c    imaginary components of the resulting fourier coefficients.
c  multivariate data is indexed according to the fortran
c    array element successor function, without limit
c    on the number of implied multiple subscripts.
c    the subroutine is called once for each variate.
c    the calls for a multivariate transform may be in any order.
c  ntot is the total number of complex data values.
c  n is the dimension of the current variable.
c  nspan/n is the spacing of consecutive data values
c    while indexing the current variable.
c  the sign of isn determines the sign of the complex
c    exponential, and the magnitude of isn is normally one.
c  a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
c    is computed by
c      call fft(a,b,n1*n2*n3,n1,n1,1)
c      call fft(a,b,n1*n2*n3,n2,n1*n2,1)
c      call fft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
c  for a single-variate transform,
c    ntot = n = nspan = (number of complex data values), e.g.
c      call fft(a,b,n,n,n,1)
c  the data can alternatively be stored in a single complex array c
c    in standard fortran fashion, i.e. alternating real and imaginary
c    parts. then with most fortran compilers, the complex array c can
c    be equivalenced to a real array a, the magnitude of isn changed
c    to two to give correct indexing increment, and a(1) and a(2) used
c    to pass the initial addresses for the sequences of real and
c    imaginary values, e.g.
c       complex c(ntot)
c       real    a(2*ntot)
c       equivalence (c(1),a(1))
c       call fft(a(1),a(2),ntot,n,nspan,2)
c  arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)
c    are used for temporary storage.  if the available storage
c    is insufficient, the program is terminated by a stop.
c    maxf must be .ge. the maximum prime factor of n.
c    maxp must be .gt. the number of prime factors of n.
c    in addition, if the square-free portion k of n has two or
c    more prime factors, then maxp must be .ge. k-1.
      implicit double precision (a-h,o-z)
      real a(1),b(1)
c  array storage in nfac for a maximum of 15 prime factors of n.
c  if n has more than one square-free factor, the product of the
c    square-free factors must be .le. 210
      dimension nfac(13),np(209)
c  array storage for maximum prime factor of 23
      dimension at(23),bt(23),ck(23),sk(23)
      equivalence (i,ii)
c  the following two constants should agree with the array dimensions.
      maxp=209
C
C Date: Wed, 9 Aug 1995 09:38:49 -0400
C From: ldm@apollo.numis.nwu.edu
      maxf=23
      maxnfac=13
C
      if(n .lt. 2) return
      inc=isn
      c72= 0.30901699437494742411
      s72= 0.95105651629515357211
      s120=0.86602540378443864675
      rad= 6.28318530717958647688
      if(isn .ge. 0) go to 10
      s72=-s72
      s120=-s120
      rad=-rad
      inc=-inc
   10 nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*float(jc)*0.5
      i=0
      jf=0
c  determine the factors of n
      m=0
      k=n
      go to 20
   15 m=m+1
      nfac(m)=4
      k=k/16
   20 if(k-(k/16)*16 .eq. 0) go to 15
      j=3
      jj=9
      go to 30
   25 m=m+1
      nfac(m)=j
      k=k/jj
   30 if(mod(k,jj) .eq. 0) go to 25
      j=j+2
      jj=j**2
      if(jj .le. k) go to 30
      if(k .gt. 4) go to 40
      kt=m
      nfac(m+1)=k
      if(k .ne. 1) m=m+1
      go to 80
   40 if(k-(k/4)*4 .ne. 0) go to 50
      m=m+1
      nfac(m)=2
      k=k/4
   50 kt=m
      j=2
   60 if(mod(k,j) .ne. 0) go to 70
      m=m+1
      nfac(m)=j
      k=k/j
   70 j=((j+1)/2)*2+1
      if(j .le. k) go to 60
   80 if(kt .eq. 0) go to 99
      j=kt
   90 m=m+1
      nfac(m)=nfac(j)
      j=j-1
      if(j .ne. 0) go to 90
   99 if (.false.) write(*,"('m='i2,' nfac ',13i2)")m,nfac
      if (m .gt. maxnfac) goto 997
c  compute fourier transform
  100 sd=radf/float(kspan)
      cd=2.0*dsin(sd)**2
      sd=dsin(sd+sd)
      kk=1
      i=i+1
      if(nfac(i) .ne. 2) go to 400
c  transform for factor of 2 (including rotation factor)
      kspan=kspan/2
      k1=kspan+2
  210 k2=kk+kspan
      ak=a(k2)
      bk=b(k2)
      a(k2)=a(kk)-ak
      b(k2)=b(kk)-bk
      a(kk)=a(kk)+ak
      b(kk)=b(kk)+bk
      kk=k2+kspan
      if(kk .le. nn) go to 210
      kk=kk-nn
      if(kk .le. jc) go to 210
      if(kk .gt. kspan) go to 800
  220 c1=1.0-cd
      s1=sd
  230 k2=kk+kspan
      ak=a(kk)-a(k2)
      bk=b(kk)-b(k2)
      a(kk)=a(kk)+a(k2)
      b(kk)=b(kk)+b(k2)
      a(k2)=c1*ak-s1*bk
      b(k2)=s1*ak+c1*bk
      kk=k2+kspan
      if(kk .lt. nt) go to 230
      k2=kk-nt
      c1=-c1
      kk=k1-k2
      if(kk .gt. k2) go to 230
      ak=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
      c1=2.0-(ak**2+s1**2)
      s1=c1*s1
      c1=c1*ak
      kk=kk+jc
      if(kk .lt. k2) go to 230
      k1=k1+inc+inc
      kk=(k1-kspan)/2+jc
      if(kk .le. jc+jc) go to 220
      go to 100
c  transform for factor of 3 (optional code)
  320 k1=kk+kspan
      k2=k1+kspan
      ak=a(kk)
      bk=b(kk)
      aj=a(k1)+a(k2)
      bj=b(k1)+b(k2)
      a(kk)=ak+aj
      b(kk)=bk+bj
      ak=-0.5*aj+ak
      bk=-0.5*bj+bk
      aj=(a(k1)-a(k2))*s120
      bj=(b(k1)-b(k2))*s120
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      kk=k2+kspan
      if(kk .lt. nn) go to 320
      kk=kk-nn
      if(kk .le. kspan) go to 320
      go to 700
c  transform for factor of 4
  400 if(nfac(i) .ne. 4) go to 600
      kspnn=kspan
      kspan=kspan/4
  410 c1=1.0
      s1=0
  420 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      akp=a(kk)+a(k2)
      akm=a(kk)-a(k2)
      ajp=a(k1)+a(k3)
      ajm=a(k1)-a(k3)
      a(kk)=akp+ajp
      ajp=akp-ajp
      bkp=b(kk)+b(k2)
      bkm=b(kk)-b(k2)
      bjp=b(k1)+b(k3)
      bjm=b(k1)-b(k3)
      b(kk)=bkp+bjp
      bjp=bkp-bjp
      if(isn .lt. 0) go to 450
      akp=akm-bjm
      akm=akm+bjm
      bkp=bkm+ajm
      bkm=bkm-ajm
      if(s1 .eq. 0) go to 460
  430 a(k1)=akp*c1-bkp*s1
      b(k1)=akp*s1+bkp*c1
      a(k2)=ajp*c2-bjp*s2
      b(k2)=ajp*s2+bjp*c2
      a(k3)=akm*c3-bkm*s3
      b(k3)=akm*s3+bkm*c3
      kk=k3+kspan
      if(kk .le. nt) go to 420
  440 c2=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
      c1=2.0-(c2**2+s1**2)
      s1=c1*s1
      c1=c1*c2
      c2=c1**2-s1**2
      s2=2.0*c1*s1
      c3=c2*c1-s2*s1
      s3=c2*s1+s2*c1
      kk=kk-nt+jc
      if(kk .le. kspan) go to 420
      kk=kk-kspan+inc
      if(kk .le. jc) go to 410
      if(kspan .eq. jc) go to 800
      go to 100
  450 akp=akm+bjm
      akm=akm-bjm
      bkp=bkm-ajm
      bkm=bkm+ajm
      if(s1 .ne. 0) go to 430
  460 a(k1)=akp
      b(k1)=bkp
      a(k2)=ajp
      b(k2)=bjp
      a(k3)=akm
      b(k3)=bkm
      kk=k3+kspan
      if(kk .le. nt) go to 420
      go to 440
c  transform for factor of 5 (optional code)
  510 c2=c72**2-s72**2
      s2=2.0*c72*s72
  520 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      k4=k3+kspan
      akp=a(k1)+a(k4)
      akm=a(k1)-a(k4)
      bkp=b(k1)+b(k4)
      bkm=b(k1)-b(k4)
      ajp=a(k2)+a(k3)
      ajm=a(k2)-a(k3)
      bjp=b(k2)+b(k3)
      bjm=b(k2)-b(k3)
      aa=a(kk)
      bb=b(kk)
      a(kk)=aa+akp+ajp
      b(kk)=bb+bkp+bjp
      ak=akp*c72+ajp*c2+aa
      bk=bkp*c72+bjp*c2+bb
      aj=akm*s72+ajm*s2
      bj=bkm*s72+bjm*s2
      a(k1)=ak-bj
      a(k4)=ak+bj
      b(k1)=bk+aj
      b(k4)=bk-aj
      ak=akp*c2+ajp*c72+aa
      bk=bkp*c2+bjp*c72+bb
      aj=akm*s2-ajm*s72
      bj=bkm*s2-bjm*s72
      a(k2)=ak-bj
      a(k3)=ak+bj
      b(k2)=bk+aj
      b(k3)=bk-aj
      kk=k4+kspan
      if(kk .lt. nn) go to 520
      kk=kk-nn
      if(kk .le. kspan) go to 520
      go to 700
c  transform for odd factors
  600 k=nfac(i)
      kspnn=kspan
      kspan=kspan/k
      if(k .eq. 3) go to 320
      if(k .eq. 5) go to 510
      if(k .eq. jf) go to 640
      jf=k
      s1=rad/float(k)
      c1=dcos(s1)
      s1=dsin(s1)
      if(jf .gt. maxf) go to 998
      ck(jf)=1.0
      sk(jf)=0.0
      j=1
  630 ck(j)=ck(k)*c1+sk(k)*s1
      sk(j)=ck(k)*s1-sk(k)*c1
      k=k-1
      ck(k)=ck(j)
      sk(k)=-sk(j)
      j=j+1
      if(j .lt. k) go to 630
  640 k1=kk
      k2=kk+kspnn
      aa=a(kk)
      bb=b(kk)
      ak=aa
      bk=bb
      j=1
      k1=k1+kspan
  650 k2=k2-kspan
      j=j+1
      at(j)=a(k1)+a(k2)
      ak=at(j)+ak
      bt(j)=b(k1)+b(k2)
      bk=bt(j)+bk
      j=j+1
      at(j)=a(k1)-a(k2)
      bt(j)=b(k1)-b(k2)
      k1=k1+kspan
      if(k1 .lt. k2) go to 650
      a(kk)=ak
      b(kk)=bk
      k1=kk
      k2=kk+kspnn
      j=1
  660 k1=k1+kspan
      k2=k2-kspan
      jj=j
      ak=aa
      bk=bb
      aj=0.0
      bj=0.0
      k=1
  670 k=k+1
      ak=at(k)*ck(jj)+ak
      bk=bt(k)*ck(jj)+bk
      k=k+1
      aj=at(k)*sk(jj)+aj
      bj=bt(k)*sk(jj)+bj
      jj=jj+j
      if(jj .gt. jf) jj=jj-jf
      if(k .lt. jf) go to 670
      k=jf-j
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      j=j+1
      if(j .lt. k) go to 660
      kk=kk+kspnn
      if(kk .le. nn) go to 640
      kk=kk-nn
      if(kk .le. kspan) go to 640
c  multiply by rotation factor (except for factors of 2 and 4)
  700 if(i .eq. m) go to 800
      kk=jc+1
  710 c2=1.0-cd
      s1=sd
  720 c1=c2
      s2=s1
      kk=kk+kspan
  730 ak=a(kk)
      a(kk)=c2*ak-s2*b(kk)
      b(kk)=s2*ak+c2*b(kk)
      kk=kk+kspnn
      if(kk .le. nt) go to 730
      ak=s1*s2
      s2=s1*c2+c1*s2
      c2=c1*c2-ak
      kk=kk-nt+kspan
      if(kk .le. kspnn) go to 730
      c2=c1-(cd*c1+sd*s1)
      s1=s1+(sd*c1-cd*s1)
      c1=2.0-(c2**2+s1**2)
      s1=c1*s1
      c2=c1*c2
      kk=kk-kspnn+jc
      if(kk .le. kspan) go to 720
      kk=kk-kspan+jc+inc
      if(kk .le. jc+jc) go to 710
      go to 100
c  permute the results to normal order---done in two stages
c  permutation for square factors of n
  800 np(1)=ks
      if(kt .eq. 0) go to 890
      k=kt+kt+1
      if(m .lt. k) k=k-1
      j=1
      np(k+1)=jc
  810 np(j+1)=np(j)/nfac(j)
      np(k)=np(k+1)*nfac(j)
      j=j+1
      k=k-1
      if(j .lt. k) go to 810
      k3=np(k+1)
      kspan=np(2)
      kk=jc+1
      k2=kspan+1
      j=1
      if(n .ne. ntot) go to 850
c  permutation for single-variate transform (optional code)
  820 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 820
  830 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 830
      j=1
  840 if(kk .lt. k2) go to 820
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 840
      if(kk .lt. ks) go to 830
      jc=k3
      go to 890
c  permutation for multivariate transform
  850 k=kk+jc
  860 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=k2+inc
      if(kk .lt. k) go to 860
      kk=kk+ks-jc
      k2=k2+ks-jc
      if(kk .lt. nt) go to 850
      k2=k2-nt+kspan
      kk=kk-nt+jc
      if(k2 .lt. ks) go to 850
  870 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 870
      j=1
  880 if(kk .lt. k2) go to 850
      kk=kk+jc
      k2=kspan+k2
      if(k2 .lt. ks) go to 880
      if(kk .lt. ks) go to 870
      jc=k3
  890 if(2*kt+1 .ge. m) return
      kspnn=np(kt+1)
c  permutation for square-free factors of n
      j=m-kt
      nfac(j+1)=1
  900 nfac(j)=nfac(j)*nfac(j+1)
      j=j-1
      if(j .ne. kt) go to 900
      kt=kt+1
      nn=nfac(kt)-1
      if(nn .gt. maxp) go to 998
      jj=0
      j=0
      go to 906
  902 jj=jj-k2
      k2=kk
      k=k+1
      kk=nfac(k)
  904 jj=kk+jj
      if(jj .ge. k2) go to 902
      np(j)=jj
  906 k2=nfac(kt)
      k=kt+1
      kk=nfac(k)
      j=j+1
      if(j .le. nn) go to 904
c  determine the permutation cycles of length greater than 1
      j=0
      go to 914
  910 k=kk
      kk=np(k)
      np(k)=-kk
      if(kk .ne. j) go to 910
      k3=kk
  914 j=j+1
      kk=np(j)
      if(kk .lt. 0) go to 914
      if(kk .ne. j) go to 910
      np(j)=-j
      if(j .ne. nn) go to 914
      maxf=inc*maxf
c  reorder a and b, following the permutation cycles
      go to 950
  924 j=j-1
      if(np(j) .lt. 0) go to 924
      jj=jc
  926 kspan=jj
      if(jj .gt. maxf) kspan=maxf
      jj=jj-kspan
      k=np(j)
      kk=jc*k+ii+jj
      k1=kk+kspan
      k2=0
  928 k2=k2+1
      at(k2)=a(k1)
      bt(k2)=b(k1)
      k1=k1-inc
      if(k1 .ne. kk) go to 928
  932 k1=kk+kspan
      k2=k1-jc*(k+np(k))
      k=-np(k)
  936 a(k1)=a(k2)
      b(k1)=b(k2)
      k1=k1-inc
      k2=k2-inc
      if(k1 .ne. kk) go to 936
      kk=k2
      if(k .ne. j) go to 932
      k1=kk+kspan
      k2=0
  940 k2=k2+1
      a(k1)=at(k2)
      b(k1)=bt(k2)
      k1=k1-inc
      if(k1 .ne. kk) go to 940
      if(jj .ne. 0) go to 926
      if(j .ne. 1) go to 924
  950 j=k3+1
      nt=nt-kspnn
      ii=nt-inc+1
      if(nt .ge. 0) go to 924
      return
  997 write(*,"('fft2 factor count limit (',i2,') exceeded')")maxnfac
c  error finish, insufficient array storage
  998 isn=0
      print 999
      stop
  999 format(44h0array bounds exceeded within subroutine fft)
      end
C-------------------------------------------------------------------
C SUBROUTINE:  REALS
C USED WITH FFT TO COMPUTE FOURIER TRANSFORM OR INVERSE
C FOR REAL DATA
C-------------------------------------------------------------------
C
      SUBROUTINE REALS(A, B, N, ISN)
C
C IF ISN=-1, THIS SUBROUTINE COMPLETES THE FOURIER TRANSFORM
C      OF 2*N REAL DATA VALUES, WHERE THE ORIGINAL DATA VALUES ARE
C      STORED ALTERNATELY IN ARRAYS A AND B, AND ARE FIRST
C      TRANSFORMED BY A COMPLEX FOURIER TRANSFORM OF DIMENSION N.
C      THE COSINE COEFFICIENTS ARE IN A(1),A(2),...A(N),A(N+1)
C      AND THE SINE COEFFICIENTS ARE IN B(1),B(2),...B(N),B(N+1).
C      NOTE THAT THE ARRAYS A AND B MUST HAVE DIMENSION N+1.
C      A TYPICAL CALLING SEQUENCE IS
C        CALL FFT  (A,B,N,N,N,-1)
C        CALL REALS(A,B,N,-1)
C
C IF ISN=1, THE INVERSE TRANSFORMATION IS DONE, THE FIRST
C      STEP IN EVALUATING A REAL FOURIER SERIES.
C      A TYPICAL CALLING SEQUENCE IS
C        CALL REALS(A,B,N,1)
C        CALL FFT  (A,B,N,N,N,1)
C      THE TIME DOMAIN RESULTS ALTERNATE IN ARRAYS A AND B,
C      I.E. A(1),B(1),A(2),B(2),...A(N),B(N).
C
C THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
C      ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
C      GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO
C      PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY
C      VALUES, E.G.
C        CALL FFT  (A,A(2),N,N,N,-2)
C        CALL REALS(A,A(2),N,-2)
C      IN THIS CASE, THE COSINE AND SINE COEFFICIENTS ALTERNATE IN A.
C
      implicit double precision (a-h,o-z)
      real a(1),b(1)
      INC = IABS(ISN)
      NF = IABS(N)
      IF (NF*ISN.NE.0) GO TO 10
      WRITE (*,9999) N, ISN
9999  FORMAT (33H ERROR - ZERO IN REALS PARAMETERS, 2I10)
      RETURN
C
  10  NK = NF*INC + 2
      NH = NK/2
      RAD = .78539816339744830961
      DR = -4.0/FLOAT(NF)
      CD = 2.0*DSIN(0.5*DR*RAD)**2
      SD = DSIN(DR*RAD)
C
C SIN,COS VALUES ARE RE-INITIALIZED EACH LIM STEPS
C
      LIM = 32
      MM = LIM
      ML = 0
      SN = 0.0
      IF (ISN.GT.0) GO TO 40
      CN = 1.0
      A(NK-1) = A(1)
      B(NK-1) = B(1)
  20  DO 30 J=1,NH,INC
        K = NK - J
        AA = A(J) + A(K)
        AB = A(J) - A(K)
        BA = B(J) + B(K)
        BB = B(J) - B(K)
        RE = CN*BA + SN*AB
        EM = SN*BA - CN*AB
        B(K) = (EM-BB)*0.5
        B(J) = (EM+BB)*0.5
        A(K) = (AA-RE)*0.5
        A(J) = (AA+RE)*0.5
        ML = ML + 1
C       IF (ML.EQ.MM) GO TO 50
        IF (ML.EQ.MM) THEN
          MM = MM + LIM
          SN = FLOAT(ML)*DR*RAD
          CN = DCOS(SN)
          IF (ISN.GT.0) CN = -CN
          SN = DSIN(SN)
          GO TO 30
        ENDIF
        AA = CN - (CD*CN+SD*SN)
        SN = (SD*CN-CD*SN) + SN
C
C THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE
C CN=AA
C
C       CN = 0.5/(AA**2+SN**2) + 0.5
C       SN = CN*SN
C       CN = CN*AA
        CN = AA
  30  CONTINUE
      RETURN
C
  40  CN = -1.0
      SD = -SD
      GO TO 20
      END
C-------------------------------------------------------------------
C SUBROUTINE:  REALT
C USED WITH 'FFT' OR ANY OTHER COMPLEX FOURIER TRANSFORM TO COMPUTE
C TRANSFORM OR INVERSE FOR REAL DATA
C THE DATA MAY BE EITHER SINGLE-VARIATE OR MULTI-VARIATE
C-------------------------------------------------------------------
C
      SUBROUTINE REALT(A, B, NSEG, N, NSPN, ISN)
C
C IF ISN=-1, THIS SUBROUTINE COMPLETES THE FOURIER TRANSFORM
C      OF 2*N REAL DATA VALUES, WHERE THE ORIGINAL DATA VALUES ARE
C      STORED ALTERNATELY IN ARRAYS A AND B, AND ARE FIRST
C      TRANSFORMED BY A COMPLEX FOURIER TRANSFORM OF DIMENSION N.
C      THE COSINE COEFFICIENTS ARE IN A(1),A(2),...A(N),A(N+1)
C      AND THE SINE COEFFICIENTS ARE IN B(1),B(2),...B(N),B(N+1).
C      NOTE THAT THE ARRAYS A AND B MUST HAVE DIMENSION N+1.
C      A TYPICAL CALLING SEQUENCE IS
C        CALL FFT  (A,B,1,N,1,-1)
C        CALL REALT(A,B,1,N,1,-1)
C
C IF ISN=1, THE INVERSE TRANSFORMATION IS DONE, THE FIRST
C      STEP IN EVALUATING A REAL FOURIER SERIES.
C      A TYPICAL CALLING SEQUENCE IS
C        CALL REALT(A,B,1,N,1,1)
C        CALL FFT  (A,B,1,N,1,1)
C      THE TIME DOMAIN RESULTS ALTERNATE IN ARRAYS A AND B,
C      I.E. A(1),B(1),A(2),B(2),...A(N),B(N).
C
C THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX
C       ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO
C       GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO
C       PASS THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY
C       VALUES, E.G.
C         CALL FFT  (A,A(2),1,N,1,-2)
C         CALL REALT(A,A(2),1,N,1,-2)
C      IN THIS CASE, THE COSINE AND SINE COEFFICIENTS ALTERNATE IN A.
C
C THIS SUBROUTINE IS SET UP TO DO THE ABOVE-DESCRIBED OPERATION ON
C      ALL SUB-VECTORS WITHIN ANY DIMENSION OF A MULTI-DIMENSIONAL
C      FOURIER TRANSFORM.  THE PARAMETERS NSEG, N, NSPN, AND INC
C      SHOULD AGREE WITH THOSE USED IN THE ASSOCIATED CALL OF 'FFT'.
C      THE FOLDING FREQUENCY COSINE COEFFICIENTS ARE STORED AT THE END
C      OF ARRAY A (WITH ZEROS IN CORRESPONDING LOCATIONS IN ARRAY B),
C      IN A SUB-MATRIX OF DIMENSION ONE LESS THAN THE MAIN ARRAY.  THE
C      DELETED DIMENSION IS THAT CORRESPONDING TO THE PARAMETER N IN
C      THE CALL OF REALT.  THUS ARRAYS A AND B MUST HAVE DIMENSION
C      NSEG*NSPN*(N+1).
C
      implicit double precision (a-h,o-z)
      real a(1),b(1)
      INC = IABS(ISN)
      KS = IABS(NSPN)*INC
      NF = IABS(N)
      NS = KS*NF
      NT = IABS(NS*NSEG)
      IF (ISN*NT.NE.0) GO TO 10
      WRITE (*,9999) NSEG, N, NSPN, ISN
9999  FORMAT (33H ERROR - ZERO IN REALT PARAMETERS, 3I10, I9)
      RETURN
C
  10  JC = KS
      K2 = IABS(KS*NSEG) - INC
      KD = NS
      NH = NS/2 + 1
      NN = NT - INC
      NT = NT + 1
      KK = 1
      RAD = .78539816339744830961
      DR = -4.0/FLOAT(NF)
      CD = 2.0*DSIN(0.5*DR*RAD)**2
      SD = DSIN(DR*RAD)
C
C SIN,COS VALUES ARE RE-INITIALIZED EACH LIM STEPS
C
      LIM = 32
      KLIM = LIM*KS
      MM = MIN0(NH,KLIM)
      SN = 0.0
      IF (ISN.GT.0) GO TO 70
C
  20  AA = A(KK)
      BA = B(KK)
      B(KK) = 0
      A(KK) = AA + BA
      A(NT) = AA - BA
      B(NT) = 0
      NT = NT + JC
      KK = KK + NS
      IF (KK.LE.NN) GO TO 20
      NT = NT - K2
      KK = KK - NN
      IF (KK.LE.JC) GO TO 20
      CN = 1.0
  30  IF (NF.EQ.1) RETURN
C
  40  AA = CN - (CD*CN+SD*SN)
      SN = (SD*CN-CD*SN) + SN
C
C THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION
C ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE
C CN=AA
C
C     CN = 0.5/(AA**2+SN**2) + 0.5
C     SN = CN*SN
C     CN = CN*AA
      CN = AA
  50  JC = JC + KS
      KD = KD - KS - KS
  60  K2 = KK + KD
      AA = A(KK) + A(K2)
      AB = A(KK) - A(K2)
      BA = B(KK) + B(K2)
      BB = B(KK) - B(K2)
      RE = CN*BA + SN*AB
      EM = SN*BA - CN*AB
      B(K2) = (EM-BB)*0.5
      B(KK) = (EM+BB)*0.5
      A(K2) = (AA-RE)*0.5
      A(KK) = (AA+RE)*0.5
      KK = KK + NS
      IF (KK.LE.NN) GO TO 60
      KK = KK - NN
      IF (KK.LE.JC) GO TO 60
      IF (KK.LE.MM) GO TO 40
      IF (KK.GT.NH) RETURN
      SN = FLOAT(JC/KS)*DR*RAD
      CN = DCOS(SN)
      IF (ISN.GT.0) CN = -CN
      SN = DSIN(SN)
      MM = MIN0(NH,MM+KLIM)
      GO TO 50
C
  70  AA = A(KK)
      BA = A(NT)
      A(KK) = (AA+BA)*0.5
      B(KK) = (AA-BA)*0.5
      NT = NT + JC
      KK = KK + NS
      IF (KK.LE.NN) GO TO 70
      NT = NT - K2
      KK = KK - NN
      IF (KK.LE.JC) GO TO 70
      CN = -1.0
      SD = -SD
      GO TO 30
      END
