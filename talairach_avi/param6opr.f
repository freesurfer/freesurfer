c $Header: /space/repo/1/dev/dev/talairach_avi/param6opr.f,v 1.1 2007/05/04 22:34:03 nicks Exp $
c $Log: param6opr.f,v $
c Revision 1.1  2007/05/04 22:34:03  nicks
c new talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  2007/04/17  05:34:45  avi
c gcc compliant
c
c Revision 1.1  1996/04/19  17:12:14  ty7777
c Initial revision
c
c     param6opr.f
c     Avi Snyder 22-Oct-93

      subroutine param6tst
      character*256 rcsid /'$Id: param6opr.f,v 1.1 2007/05/04 22:34:03 nicks Exp $'/
      real*4 param0(6),param1(6),param2(6)
      data param0/.1,.2,.3,.4,.5,.6/

      call param6inv(param0,param1)
      l=lnblnk(rcsid)
      write(*, "(a)")rcsid(1:l)
      write(*, "(6f10.4)")param1
      call param6mul(param0,param1,param2)
      write(*, "(6f10.4)")param2
      call param6inv(param1,param0)
      write(*, "(6f10.4)")param0
      do 1 i=1,6
    1 param0(i)=-2.*param0(i)
      call param6inv(param0,param1)
      write(*, "(6f10.4)")param1
      call param6mul(param0,param1,param2)
      write(*, "(6f10.4)")param2
      call param6inv(param1,param0)
      write(*, "(6f10.4)")param0
      end
      subroutine ang2rot(ang,rot)
      real*4 ang(3),rot(3,3)
      c1=cos(ang(1))
      s1=sin(ang(1))
      c2=cos(ang(2))
      s2=sin(ang(2))
      c3=cos(ang(3))
      s3=sin(ang(3))
      rot(1,1)= c3*c2
      rot(1,2)=-c3*s2*s1-s3*c1
      rot(1,3)=-c3*s2*c1+s3*s1
      rot(2,1)= s3*c2
      rot(2,2)=-s3*s2*s1+c3*c1
      rot(2,3)=-s3*s2*c1-c3*s1
      rot(3,1)= s2
      rot(3,2)= c2*s1
      rot(3,3)= c2*c1
      return
      end
      subroutine rot2ang(rot,ang)
      real*4 rot(3,3),ang(3)
      ang(1)=atan2(rot(3,2),rot(3,3))
      c1=cos(ang(1))
      s1=sin(ang(1))
      if(abs(c1).gt.abs(s1))then
        ang(2)=atan2(rot(3,1),rot(3,3)/c1)
      else
        ang(2)=atan2(rot(3,1),rot(3,2)/s1)
      endif
      ang(3)=atan2(rot(2,1),rot(1,1))
      return
      end
      subroutine param6mul(param1,param2,param3)
c     computes equivalent of param1 followed by param2
c     and leaves in param3
      real*4 param1(6),param2(6),param3(6)
      real*4 rot1(3,3),rot2(3,3),rot3(3,3)
      call ang2rot(param1(4),rot1)
      call ang2rot(param2(4),rot2)
      do 1 i=1,3
      do 1 j=1,3
      rot3(i,j)=0.
      do 1 k=1,3
    1 rot3(i,j)=rot3(i,j)+rot2(i,k)*rot1(k,j)
      do 2 i=1,3
      param3(i)=param2(i)
      do 2 j=1,3
    2 param3(i)=param3(i)+rot2(i,j)*param1(j)
      call rot2ang(rot3,param3(4))
      return
      end
      subroutine param6inv(param0,param1)
c     computes rotation-translation inverse of param0
c     and leaves in param1
      real*4 param0(6),param1(6)
      real*4 rot(3,3)
      call ang2rot(param0(4),rot)
      do 1 i=1,3
      param1(i)=0.
      do 1 j=1,3
    1 param1(i)=param1(i)-rot(j,i)*param0(j)
      param1(4)=atan2(rot(2,3),rot(3,3))
      param1(5)=atan2(rot(1,3),rot(3,3)/cos(param1(4)))
      param1(6)=atan2(rot(1,2),rot(1,1))
      return
      end
      SUBROUTINE img2vrt(MMPPIX,CENTER,C)
      REAL*4 MMPPIX(3),CENTER(3),C(4,4)
      DO 1 I=1,4
      DO 1 J=1,4
    1 C(I,J)=0.
      DO 2 I=1,3
      C(I,I)=MMPPIX(I)
    2 C(I,4)=-CENTER(I)
      C(4,4)=1.
      RETURN
      END
      SUBROUTINE vrt2img(MMPPIX,CENTER,C)
      REAL*4 MMPPIX(3),CENTER(3),C(4,4)
      DO 1 I=1,4
      DO 1 J=1,4
    1 C(I,J)=0.
      DO 2 I=1,3
      C(I,I)=1./MMPPIX(I)
    2 C(I,4)=CENTER(I)/MMPPIX(I)
      C(4,4)=1.
      RETURN
      END
      SUBROUTINE TROTCAL(IFWDREV,PARAM,TROT)
      REAL*4 PARAM(6),TROT(4,4)
      REAL*4 ROT(3,3)
      CALL ang2rot(PARAM(4),ROT)
      IF(IFWDREV.GE.0)THEN
        DO I=1,3
          DO J=1,3
            TROT(I,J)=ROT(I,J)
          ENDDO
          TROT(I,4)=PARAM(I)
        ENDDO
      ELSE
        DO I=1,3
          DO J=1,3
            TROT(I,J)=ROT(J,I)
          ENDDO
          TROT(I,4)=0.
          DO K=1,3
            TROT(I,4)=TROT(I,4)-ROT(K,I)*PARAM(K)
          ENDDO
        ENDDO
      ENDIF
      DO 7 J=1,3
    7 TROT(4,J)=0.
      TROT(4,4)=1.
      RETURN
      END
      SUBROUTINE TROTSET(DIST,ANG,TROT)
      REAL*4 DIST(3),ANG(3),TROT(4,4)
      REAL*4 ROT(3,3)
      CALL ang2rot(ANG,ROT)
      DO I=1,3
        DO J=1,3
          TROT(I,J)=ROT(I,J)
        ENDDO
        TROT(I,4)=DIST(I)
      ENDDO
      DO 7 J=1,3
    7 TROT(4,J)=0.
      TROT(4,4)=1.
      RETURN
      END
      subroutine stretchset(stretch,strang,t)
      real*4 stretch(3),strang(3),t(4,4)
      real*4 a(3,3),b(3,3),c(3,3)
      call ang2rot(strang,a)
      do 1 i=1,3
      do 2 j=1,3
    2 b(i,j)=0.
    1 b(i,i)=exp(stretch(i))
      do 3 i=1,3
      do 3 j=1,3
      c(i,j)=0.
      do 3 k=1,3
    3 c(i,j)=c(i,j)+a(k,i)*b(k,j)
      do 4 i=1,3
      do 5 j=1,3
      t(i,j)=0.
      do 5 k=1,3
    5 t(i,j)=t(i,j)+c(i,k)*a(k,j)
      t(i,4)=0.
    4 t(4,i)=0.
      t(4,4)=1.
      return
      end
