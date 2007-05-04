c$Header: /space/repo/1/dev/dev/talairach_avi/ft4imgn.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c$Log: ft4imgn.f,v $
cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  2007/04/29  05:15:21  avi
c gcc v3 compliant (use set_rnan() instead of r_quiet_nan())
c
c Revision 1.1  2006/02/13  06:24:18  avi
c Initial revision
c
      subroutine ft4imgn(t4,imgt,nxt,nyt,nzt,centert,mmppixt,imgo,nxo,nyo,nzo,centero,mmppixo)
c     variant of ft4imgo for nearest neighbor resampling
      real*4 t4(4,4)
      real*4 imgt(nxt,nyt,nzt)
      real*4 centert(3),mmppixt(3)
      real*4 imgo(nxo,nyo,nzo)
      real*4 centero(3),mmppixo(3)
      character*256 rcsid/'$Id: ft4imgn.f,v 1.1 2007/05/04 22:33:59 nicks Exp $'/
      real*4 a(4,4),b(4,4),xk(3),xt(3)
      external matmul

      call img2vrt(mmppixo,centero,a)
      call matmul(t4,a,b,4)

      do 41 kz=1,nzo
      xk(3)=float(kz)
      do 41 ky=1,nyo
      xk(2)=float(ky)
      do 41 kx=1,nxo
      xk(1)=float(kx)
      do k=1,3
        xt(k)=b(k,1)*xk(1)+b(k,2)*xk(2)+b(k,3)*xk(3)+b(k,4)
      enddo
      call imgvaln(imgt,nxt,nyt,nzt,centert,mmppixt,xt,v,lslice)
      if(lslice.gt.0)then
        imgo(kx,ky,kz)=v
      else
        call set_rnan(imgo(kx,ky,kz))
      endif
   41 continue

      return
      end

      subroutine imgvaln(imgt,nx,ny,nz,center,mmppix,x,v,lslice)
c     variant of imgvalx
c     avi snyder 
c     returns in v the nearest neighbor value of imgt at locus x.
c     lslice returns the slice from which the data is taken.
c     lslice = 0 (and v = 0.) if imgt is undefined at locus x.

      real*4 imgt(nx,ny,nz)
      real*4 center(3),mmppix(3),x(3)
      logical*4 ldefined,lok
      real*4 xf(3)
      integer*4 ix(3),mx(3)

      mx(1)=nx
      mx(2)=ny
      mx(3)=nz

      ldefined=.true.
      do 1 k=1,3
      center(k)=sign(center(k),mmppix(k))
      xf(k)=(center(k)+x(k))/mmppix(k)
      ix(k)=nint(xf(k))
    1 ldefined=ldefined.and.ix(k).gt.0.and.ix(k).le.mx(k)
      if(.not.ldefined)then
        v=0.
        lslice=0
        return
      endif

      v=imgt(ix(1),ix(2),ix(3))
      lok=(v.ge.0.0.or.v.le.0.0)
      if(.not.lok)then
        v=0.
        lslice=0
        return
      endif

      lslice=ix(3)
      return
      end
