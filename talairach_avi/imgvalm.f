cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.1  1998/10/02  20:23:27  mcavoy
c Initial revision
c
c Revision 1.1  1997/11/14  17:16:24  tscull
c Initial revision
c
c Revision 1.2  1996/12/28  01:41:46  avi
c correct error in y bound checking and make algorithm
c more like imgvalx.f
c
c Revision 1.1  1996/12/27  23:26:08  avi
c Initial revision
c
      subroutine imgvalm(imag,mask,nx,ny,nz,t,x,v,lslice)
c     returns in v imag value at t transformed coordinate x
c     lslice=-1 coordinate not in imag array
c     lslice= 0 coordinate not in mask
      character*256 rcshdr
     &/'$Header: /space/repo/1/dev/dev/talairach_avi/imgvalm.f,v 1.1 2007/05/04 22:33:59 nicks Exp $'/
      real*4 imag(nx,ny,nz)
      integer*2 mask(nx,ny,nz)
      real*4 t(4,4),x(3)
      real*4 xp(3)
      logical*4 l2d

      l2d=nz.eq.1
      do l=1,3
        xp(l)=t(l,1)*x(1)+t(l,2)*x(2)+t(l,3)*x(3)+t(l,4)
      enddo
      ix=nint(xp(1)-0.5)
      wx=xp(1)-float(ix)
      if(ix.eq.0.and.wx.gt.0.999)then
        ix=1
        wx=0.
      endif
      if(ix.eq.nx.and.wx.lt.0.001)then
        ix=nx-1
        wx=1.
      endif
      if(ix.lt.1.or.ix.ge.nx)goto 9
      iy=nint(xp(2)-0.5)
      wy=xp(2)-float(iy)
      if(iy.eq.0.and.wy.gt.0.999)then
        iy=1
        wy=0.
      endif
      if(iy.eq.ny.and.wy.lt.0.001)then
        iy=ny-1
        wy=1.
      endif
      if(iy.lt.1.or.iy.ge.ny)goto 9
      if(l2d)goto 2
      iz=nint(xp(3)-0.5)
      wz=xp(3)-float(iz)
      if(iz.eq.0.and.wz.gt.0.999)then
        iz=1
        wz=0.
      endif
      if(iz.eq.nz.and.wz.lt.0.001)then
        iz=nz-1
        wz=1.
      endif
      if(iz.lt.1.or.iz.ge.nz)goto 9
      if(mask(ix,iy,iz).eq.0)then
        lslice=0
        v=0.
        return
      endif
      v=
     &  (1.-wz)*((1.-wy)*(
     &               (1.-wx)*imag(ix+0,iy+0,iz+0)
     &                   +wx*imag(ix+1,iy+0,iz+0))
     &               +wy*(
     &               (1.-wx)*imag(ix+0,iy+1,iz+0)
     &                   +wx*imag(ix+1,iy+1,iz+0)))
     &      +wz*((1.-wy)*(
     &               (1.-wx)*imag(ix+0,iy+0,iz+1)
     &                   +wx*imag(ix+1,iy+0,iz+1))
     &               +wy*(
     &               (1.-wx)*imag(ix+0,iy+1,iz+1)
     &                   +wx*imag(ix+1,iy+1,iz+1)))
      lslice=iz
      if(wz.gt.0.5)lslice=lslice+1
      return
    9 lslice=-1
      v=0.
      return

    2 if(mask(ix,iy,1).eq.0)then
        lslice=0
        v=0.
        return
      endif
      v=
     &      (1.-wy)*((1.-wx)*imag(ix+0,iy+0,1)
     &                   +wx*imag(ix+1,iy+0,1))
     &          +wy*((1.-wx)*imag(ix+0,iy+1,1)
     &                   +wx*imag(ix+1,iy+1,1))
      lslice=1
      return
      end

      subroutine imgvalmg(imag,mask,nx,ny,nz,t,x,eps,grad,lslice)
c     returns grad imag value at t transformed coordinate x
c     lslice=-1 coordinate not in imag array
c     lslice= 0 coordinate not in mask
      real*4 imag(nx,ny,nz)
      integer*2 mask(nx,ny,nz)
      real*4 t(4,4),x(3),grad(3)
      real*4 eps		! +/- interval over which to evaluate grad
      real*4 xe(3)
      logical*4 l2d

      l2d=nz.eq.1
      if(l2d)then
        nk=2
        grad(3)=0.
      else
        nk=3
      endif

      do 2 k=1,nk
    2 xe(k)=x(k)
      s=0.
      do 1 k=1,nk
      xe(k)=x(k)+eps
      call imgvalm(imag,mask,nx,ny,nz,t,xe,v1,lok)
      if(lok.le.0)goto 9
      s=s+float(lok)
      xe(k)=x(k)-eps
      call imgvalm(imag,mask,nx,ny,nz,t,xe,v0,lok)
      if(lok.le.0)goto 9
      s=s+float(lok)
      grad(k)=0.5*(v1-v0)/eps
    1 xe(k)=x(k)

      lslice=nint(0.5*s/float(nk))
      return
    9 lslice=-1
      return
      end
