c$Header: /space/repo/1/dev/dev/talairach_avi/ft4ixyz.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c$Log: ft4ixyz.f,v $
cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  1999/03/11  01:36:38  avi
c consolidate d2xt, d2yt, d2zt into imgt(,,,4)
c
c Revision 1.1  1999/01/07  07:19:33  avi
c Initial revision
c
      subroutine ft4ixyz_rcs
      write(*,"('$Id: ft4ixyz.f,v 1.1 2007/05/04 22:33:59 nicks Exp $')")
      return
      end
      
      subroutine pt4ixyz(imgdim,voxdim,centert,mmppixt,t4mat,centero,mmppixo,t4atl,t4)
      integer*4 imgdim(3)			! assumed positive
      real*4    voxdim(3)			! assumed positive
      real*4    centert(3),mmppixt(3),centero(3),mmppixo(3),t4atl(4,4),t4mat(4,4),t4(4,4)
      real*4 o(4,4),t(4,4),temp1(4,4),temp2(4,4),e(4,4),f(4,4),centere(3)
      real*4 f2c(4,4)/1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,-1.0,-1.0,-1.0,1./
      external matmul,img2vrt,vrt2img

      do 1 k=1,3
    1 centere(k)=voxdim(k)*float(imgdim(k)/2)	! as in cross_realign3d_4dfp

      call img2vrt(voxdim,centere,e)		! as in cross_realign3d_4dfp
      call vrt2img(voxdim,centere,f)		! as in cross_realign3d_4dfp
      call img2vrt(mmppixo,centero,o)		! as in t4imgs_4dfp
      call vrt2img(mmppixt,centert,t)		! as in t4imgs_4dfp

      call matmul(t4atl,o,temp1,4)		! target (atlas) FORTRAN index to virtual
      call matmul(t,temp1,temp2,4)		! to source FORTRAN index
      call matmul(f2c,temp2,temp1,4)		! to source C index
      call matmul(e,temp1,temp2,4)		! to virtual
      call matmul(t4mat,temp2,temp1,4)		! frame align as computed by cross_realign3d_4dfp
      call matmul(f,temp1,t4,4)			! to source C index
      return
      end

      subroutine ft4ixyz(mode,t4,imgt,nxt,nyt,nzt,imgo,nxo,nyo,nzo)
      real*4 t4(4,4)				! target FORTRAN index to source C index
      real*4 imgt(nxt,nyt,nzt,4)		! source image and spline coefficients
      real*4 imgo(nxo,nyo,nzo)			! target
      real*4 xo(3),xt(3)
      logical*4 lwrap,lendslc,leval
      logical*4 ldebug/.false./

      lwrap=  iand(mode,1024).ne.0
      lendslc=iand(mode,2048).ne.0

      do 41 kz=1,nzo
      xo(3)=float(kz)
      do 41 ky=1,nyo
      xo(2)=float(ky)
      do 41 kx=1,nxo
      xo(1)=float(kx)
      do k=1,3
        xt(k)=t4(k,1)*xo(1)+t4(k,2)*xo(2)+t4(k,3)*xo(3)+t4(k,4)
      enddo
      leval=.true.
      if(.not.lwrap)then
        if(xt(1).lt.-0.001.or.xt(1).gt.float(nxt-1)+0.001)leval=.false.
        if(xt(2).lt.-0.001.or.xt(2).gt.float(nyt-1)+0.001)leval=.false.
      endif
      if(leval.and.lendslc)then
        if(xt(3).lt.-0.001.or.xt(3).gt.float(nzt-1)+0.001)leval=.false.
      else
        if(xt(3).lt.-0.500.or.xt(3).gt.float(nzt-1)+0.500)leval=.false.
      endif
      if(leval)then
        call splint3dv(imgt,nxt,nyt,nzt,imgt(1,1,1,2),imgt(1,1,1,3),imgt(1,1,1,4),xt(1),xt(2),xt(3),v)
        imgo(kx,ky,kz)=v
      else
        call set_rnan(imgo(kx,ky,kz))
      endif
   41 continue

      return
      end
