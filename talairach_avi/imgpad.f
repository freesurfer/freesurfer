cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Header: /space/repo/1/dev/dev/talairach_avi/imgpad.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c $Log: imgpad.f,v $
c Revision 1.1  2007/05/04 22:33:59  nicks
c new talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.4  2007/04/11  05:15:53  avi
c f90 compliant
c
c Revision 1.3  2000/12/13  18:46:45  avi
c copyright
c
c Revision 1.2  1996/04/25  22:02:02  ty7777
c Fix bug for case of margin equals zero.
c
c Revision 1.1  1996/04/19  17:10:57  ty7777
c Initial revision
c
      subroutine npad_test
      real*4 mmppix/2.086214/
c     margin computed as 2.*sd = 2.*0.1874/(mmppix*fhalf)
      fhalf=0.02
      margin=nint(2.*0.1874/(mmppix*fhalf))
      write(*, "('mmppix,fhalf,margin',2f10.4,i4)")mmppix,fhalf,margin
      do 1 n=1,1000
      np=npad(n,margin)
    1 write (*, "(2i10,f10.4)")n,np,float(np)/float(n)
      end

      function npad(n,margin)
      if(n.le.1)then
        npad=n
        return
      endif
      m=1
      do 2 j=1,12
      do 4 i=2,9
      npad=m*i
    4 if(i.ne.7.and.npad.ge.n+2*margin)return
    2 m=m*2
      npad=n
      return
      end

      subroutine imgpad(imag,nx,ny,nz,imagp,nxp,nyp,nzp)
      real*4 imag(nx,ny,nz),imagp(nxp,nyp,nzp)
      logical*4 ldebug/.false./

      pi=atan2(0.,-1.)

      mx=(nxp-nx)/2
      my=(nyp-ny)/2
      mz=(nzp-nz)/2
      l=nzp-nz-2*mz	! 0 or 1
      do 1 kp=1,nzp
      k=kp-mz
c     if(k.ge.-mz/2.and.k.lt.nz+mz/2)k=min0(nz,max0(1,kp-mz))
      t=0.
      if(k.lt.1)then
        if(mz.gt.0)t=float(k-1)/float(mz)
        k=1
      elseif(k.gt.nz)then
        if(mz.gt.0)t=float(nz-k+l)/float(mz)
        k=nz
      endif
      f=0.5*(1.+cos(pi*t))
      if(ldebug)write (*, "('kp,k,f ',2i6,f10.6)")kp,k,f
      do 1 jp=1,nyp
      j=jp-my
      do 1 ip=1,nxp
      i=ip-mx
      if(j.lt.1.or.j.gt.ny.or.i.lt.1.or.i.gt.nx)then
        imagp(ip,jp,kp)=0.
      else
        imagp(ip,jp,kp)=f*imag(i,j,k)
      endif
    1 continue

      return
      end

      subroutine imgdap(imag,nx,ny,nz,imagp,nxp,nyp,nzp)
      real*4 imag(nx,ny,nz),imagp(nxp,nyp,nzp)

      mx=(nxp-nx)/2
      my=(nyp-ny)/2
      mz=(nzp-nz)/2

      do 1 k=1,nz
      kp=k+mz
      do 1 j=1,ny
      jp=j+my
      do 1 i=1,nx
      ip=i+mx
    1 imag(i,j,k)=imagp(ip,jp,kp)
      return
      end

