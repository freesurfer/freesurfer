c$Header: /space/repo/1/dev/dev/talairach_avi/ft4imgo.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c$Log: ft4imgo.f,v $
cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.3  2007/04/29  05:13:28  avi
c gcc v3 compliant (use set_rnan() instead of r_quite_nan())
c
c Revision 1.2  1999/01/04  01:47:37  avi
c subroutine vrtflip
c
c Revision 1.1  1997/10/11  00:17:27  avi
c Initial revision
c
c Revision 1.2  1997/06/18  05:32:25  avi
c fix mistake x'8000000' (negative zero) to correct r_quiet_nan()
c
c Revision 1.1  1997/01/04  09:02:09  avi
c Initial revision
c
      subroutine vrtflip(iori,imgdim,centeri,mmppixi,centert,mmppixt)
      integer*4 imgdim(3)
      real*4 mmppixi(3),centeri(3)
      real*4 mmppixt(3),centert(3)
      real*4 flips(3,3)/-1.,+1.,-1.,		! transverse
     &                  -1.,+1.,+1.,		! coronal
     &                  +1.,+1.,+1./		! sagittal

      k=iori-1					! ifh.orientation
      do 1 i=1,3
      mmppixt(i)=mmppixi(i)*flips(i,k)
      centert(i)=centeri(i)*flips(i,k)
      if(flips(i,k).lt.0.)centert(i)=mmppixt(i)*float(imgdim(i)+1)-centert(i)
    1 continue
      return
      end

      subroutine ft4imgo(t4,imgt,nxt,nyt,nzt,centert,mmppixt,imgo,nxo,nyo,nzo,centero,mmppixo)
      real*4 t4(4,4)
      real*4 imgt(nxt,nyt,nzt)
      real*4 centert(3),mmppixt(3)
      real*4 imgo(nxo,nyo,nzo)
      real*4 centero(3),mmppixo(3)
      character*256 rcsid/'$Id: ft4imgo.f,v 1.1 2007/05/04 22:33:59 nicks Exp $'/
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
      call imgvalx(imgt,nxt,nyt,nzt,centert,mmppixt,xt,v,lslice)
      if(lslice.gt.0)then
        imgo(kx,ky,kz)=v
      else
        call set_rnan(imgo(kx,ky,kz))
      endif
   41 continue

      return
      end
