c$Header: /space/repo/1/dev/dev/talairach_avi/ffind.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c$Log: ffind.f,v $
cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  2007/04/25  05:17:11  avi
c gcc v3 compliant
c
      subroutine ffind(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &                   img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode)
c     find optimal translation to achieve registration of img2 on img1
c     by brute force coarse search in translation space
      real*4 img1(nx1,ny1,nz1)
      integer*2 msk1(nx1,ny1,nz1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(nx2,ny2,nz2)
      integer*2 msk2(nx2,ny2,nz2)
      real*4 mmppix2(3),center2(3)
      real*4 param(13)
      parameter (dd=7.5)		! linear distance search increment in units of mmppix
      parameter (nterm=3)
      real*4 taram(nterm),garam(nterm)
      logical*4 lenable,l3d
      logical*4 ldebug/.true./

      if(nz2.lt.2)mode=iand(mode,not(2))	! prevent attempt to 3d align
      lenable=iand(mode,1).ne.0
      l3d=    iand(mode,2).ne.0
      write(*,"('image alignment mode ',i6,' decimal ',z8,' hex')")mode,mode
      if(.not.lenable)return

      do jj=1,nterm
        taram(jj)=param(jj)
      enddo

      if(l3d)then
        nk=6
      else
        nk=0
      endif
      nj=6
      ni=1

      etamax=0.
      do 81 k=-nk,nk
      param(3)=taram(3)+float(k)*dd
      do 81 j=-nj,nj
      param(2)=taram(2)+float(j)*dd
      do 81 i=-ni,ni
      param(1)=taram(1)+float(i)*dd
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,9.0,eta,q)
      write(*,"(6f10.4)")(param(jj),jj=1,nterm),eta,q
      if(eta.gt.etamax)then
        etamax=eta
        do jj=1,nterm
          garam(jj)=param(jj)
        enddo
      endif
   81 continue

      do jj=1,nterm
        param(jj)=garam(jj)
      enddo
      
      return
      end

