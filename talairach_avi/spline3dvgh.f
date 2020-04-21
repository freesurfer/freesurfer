cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
c Washington University, Mallinckrodt Institute of Radiology.
c All Rights Reserved.
c This software may not be reproduced, copied, or distributed without written
c permission of Washington University. For further information contact A. Z. Snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Log: spline3dvgh.f,v $
c Revision 1.1  2007/05/04 22:34:03  nicks
c new talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.10  2007/04/22  02:24:17  avi
c remove pointers
c
c Revision 1.9  2007/04/01  00:27:11  avi
c further prepartion for transition to f90
c
c Revision 1.8  2007/03/31  03:38:06  avi
c eliminate call to rand(0) in splint3dvgh_testh and splint3dvgh_testg
c
c Revision 1.7  2006/03/09  06:10:10  avi
c correct splint3dvgh_testg()
c
c Revision 1.6  2000/12/13  02:55:55  avi
c copyright
c
c Revision 1.5  1999/03/09  06:17:12  avi
c splint3dvgl splint3dvl splint3dvl_test
c
c Revision 1.4  1998/05/20  10:01:14  avi
c minor changes to splinezf - dummy array length should not be numerically undefined
c
c Revision 1.3  1998/05/18  05:42:03  avi
c splinezf
c
c Revision 1.2  1998/04/20  01:09:20  avi
c correct evaluations for x and y wrpa addressing
c
c Revision 1.1  1998/04/20  00:47:49  avi
c Initial revision
c
c $Header: /space/repo/1/dev/dev/talairach_avi/spline3dvgh.f,v 1.1 2007/05/04 22:34:03 nicks Exp $

c     call splint3dv_test
c     call splint3dvl_test
c     call splint3dvgh_testh
c     call splint3dvgh_testg
c     call spline3dvgh_rcs
c     end

      subroutine spline3dvgh_rcs
      write (*,"('spline3dvgh.f')")
      return
      end

      subroutine splint3dvl_test
      parameter (nx=15)
      parameter (ny=9)
      parameter (nz=5)
      real*4 f(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2zf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d(-nx/2:nx/2)
      real*4 g(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2)
      real*4 h(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2)
      data del/0.2/

      write (*,"('splint3dvl_test')")
      ix0=nx/2
      iy0=ny/2
      iz0=nz/2
      do 2 iz=0,nz-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy,iz)=0.
      f(ix0,iy0,iz0)=1.

      call splinex (f,nx,ny,nz,d2xf)
      call spliney (f,nx,ny,nz,d2yf)
      call splineza(f,nx,ny,nz,d2zf)

      do 3 iz=-nz/2,nz/2
      do 3 iy=-ny/2,ny/2
      do 3 ix=-nx/2,nx/2
      z=float(iz0)+del*float(iz)
      y=float(iy0)+del*float(iy)
      x=float(ix0)+del*float(ix)
    3 call splint3dvl(f,nx,ny,nz,d2xf,d2yf,d2zf,x,y,z,g(ix,iy,iz),h(ix,iy,iz))

      do 4 ix=-nx/2,nx/2
    4 d(ix)=float(ix0)+del*float(ix)

      write (*,"('interpolation')")
      do 51 iz=-nz/2,nz/2
      write (*,"('iz=',i2,' z=',f5.2)")iz,float(iz0)+del*float(iz)
      write (*,"(6x,17f5.2)")d
      write (*,"(6x,17('-----'))")
      do iy=-ny/2,ny/2
        write (*,"(f5.2,'|',17f5.2)")float(iy0)+del*float(iy),(g(ix,iy,iz),ix=-nx/2,nx/2)
      enddo
   51 continue

      write (*,"('laplacian')")
      do 61 iz=-nz/2,nz/2
      write (*,"('iz=',i2,' z=',f5.2)")iz,float(iz0)+del*float(iz)
      write (*,"(6x,17f5.2)")d
      write (*,"(6x,17('-----'))")
      do iy=-ny/2,ny/2
        write (*,"(f5.2,'|',17f5.2)")float(iy0)+del*float(iy),(h(ix,iy,iz),ix=-nx/2,nx/2)
      enddo
   61 continue

      return
      end

      subroutine splint3dv_test
      parameter (nx=15)
      parameter (ny=9)
      parameter (nz=5)
      real*4 f(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2zf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yf(0:nx-1,0:ny-1,0:nz-1)
      real*4 g(-nx/2:nx/2,-ny/2:ny/2,-nz/2:nz/2)
      real*4 d(-nx/2:nx/2)
      data del/0.2/

      write (*,"('splint3dv_test')")
      ix0=nx/2
      iy0=ny/2
      iz0=nz/2
      do 2 iz=0,nz-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy,iz)=0.
      f(ix0,iy0,iz0)=1.

      call splinex (f,nx,ny,nz,d2xf)
      call spliney (f,nx,ny,nz,d2yf)
      call splineza(f,nx,ny,nz,d2zf)

      do 3 iz=-nz/2,nz/2
      do 3 iy=-ny/2,ny/2
      do 3 ix=-nx/2,nx/2
      z=float(iz0)+del*float(iz)
      y=float(iy0)+del*float(iy)
      x=float(ix0)+del*float(ix)
      call splint3dv(f,nx,ny,nz,d2xf,d2yf,d2zf,x,y,z,v)
    3 g(ix,iy,iz)=v

      do 4 ix=-nx/2,nx/2
    4 d(ix)=float(ix0)+del*float(ix)

      write (*,"('g')")
      do 51 iz=-nz/2,nz/2
      write (*,"('iz=',i2,' z=',f5.2)")iz,float(iz0)+del*float(iz)
      write (*,"(6x,17f5.2)")d
      write (*,"(6x,17('-----'))")
      do iy=-ny/2,ny/2
        write (*,"(f5.2,'|',17f5.2)")float(iy0)+del*float(iy),(g(ix,iy,iz),ix=-nx/2,nx/2)
      enddo
   51 continue

      write (*,"('d2xf')")
      do 21 iz=0,nz-1
      write (*,"('iz=',i1)")iz
      do iy=0,ny-1
        write (*,"(i3,17f5.2)")iy,(d2xf(ix,iy,iz),ix=0,nx-1)
      enddo
   21 continue

      write (*,"('d2yf')")
      do 31 iz=0,nz-1
      write (*,"('iz=',i1)")iz
      do iy=0,ny-1
        write (*,"(i3,17f5.2)")iy,(d2yf(ix,iy,iz),ix=0,nx-1)
      enddo
   31 continue

      write (*,"('d2zf')")
      do 41 iz=0,nz-1
      write (*,"('iz=',i1)")iz
      do iy=0,ny-1
        write (*,"(i3,17f5.2)")iy,(d2zf(ix,iy,iz),ix=0,nx-1)
      enddo
   41 continue

      return
      end

      subroutine splinex(imag,nx,ny,nz,d2xi)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:511),b(0:511),f(0:256)
c     real*4 a(0:nx-1),b(0:nx-1),f(0:nx)
c     pointer (pa,a),(pb,b),(pf,f)

c     pa=malloc(4*nx)
c     pb=malloc(4*nx)
c     pf=malloc(4*(nx+1))
      if (nx.gt.512) stop 'splinex: nx may not exceed 512'

      twopi=8.*atan(1.)
      do ix=0,(nx+1)/2
        c=cos(twopi*float(ix)/float(nx))
        f(ix)=6.*(c-1.)/(c+2.)
      enddo

      do 21 iz=0,nz-1
      do 21 iy=0,ny-1
      do ix=0,nx-1
        a(ix)=imag(ix,iy,iz)
        b(ix)=0.
      enddo
      call fft(a,b,1,nx,1,+1)
      a(0)=0.
      do ix=1,(nx+1)/2-1
        a(ix)=a(ix)*f(ix)
        b(ix)=b(ix)*f(ix)
        a(nx-ix)=a(nx-ix)*f(ix)
        b(nx-ix)=b(nx-ix)*f(ix)
      enddo
      if(mod(nx,2).eq.0)a(nx/2)=a(nx/2)*(-12.)
      call fft(a,b,1,nx,1,-1)
      do ix=0,nx-1
        d2xi(ix,iy,iz)=a(ix)
      enddo
   21 continue

c     call free(pa)
c     call free(pb)
c     call free(pf)
      return
      end

      subroutine spliney(imag,nx,ny,nz,d2yi)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2yi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:511),b(0:511),f(0:256)
c     real*4 a(0:ny-1),b(0:ny-1),f(0:ny)
c     pointer (pa,a),(pb,b),(pf,f)

c     pa=malloc(4*ny)
c     pb=malloc(4*ny)
c     pf=malloc(4*(ny+1))
      if (ny.gt.512) stop 'spliney: ny may not exceed 512'

      twopi=8.*atan(1.)
      do iy=0,(ny+1)/2
        c=cos(twopi*float(iy)/float(ny))
        f(iy)=6.*(c-1.)/(c+2.)
      enddo

      do 31 iz=0,nz-1
      do 31 ix=0,nx-1
      do iy=0,ny-1
        a(iy)=imag(ix,iy,iz)
        b(iy)=0.
      enddo
      call fft(a,b,1,ny,1,+1)
      a(0)=0.
      do iy=1,(ny+1)/2-1
        a(iy)=a(iy)*f(iy)
        b(iy)=b(iy)*f(iy)
        a(ny-iy)=a(ny-iy)*f(iy)
        b(ny-iy)=b(ny-iy)*f(iy)
      enddo
      if(mod(ny,2).eq.0)a(ny/2)=a(ny/2)*(-12.)
      call fft(a,b,1,ny,1,-1)
      do iy=0,ny-1
        d2yi(ix,iy,iz)=a(iy)
      enddo
   31 continue

c     call free(pa)
c     call free(pb)
c     call free(pf)
      return
      end

      subroutine splinezf(imag,nx,ny,nz,d2zi)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:511),b(0:511),f(0:256)
c     real*4 a(0:2*nz),b(0:2*nz),f(0:2*nz)
c     pointer (pa,a),(pb,b),(pf,f)
      logical*4 ldebug/.true./

      nzp=npad(nz,4)
      if(ldebug)write (*,"('nz',i4,' padded to',i4)")nz,nzp
c     pa=malloc(4*nzp)
c     pb=malloc(4*nzp)
c     pf=malloc(4*(nzp+1))
      if (nzp.gt.512) stop 'splinezf: nzp may not exceed 512'

      twopi=8.*atan(1.)
      do iz=0,(nzp+1)/2
        c=cos(twopi*float(iz)/float(nzp))
        f(iz)=6.*(c-1.)/(c+2.)
      enddo

      do 31 iy=0,ny-1
      do 31 ix=0,nx-1
      do iz=0,nzp-1
        a(iz)=0.
        b(iz)=0.
      enddo
      do iz=0,nz-1
        a(iz)=imag(ix,iy,iz)
      enddo
      call fft(a,b,1,nzp,1,+1)
      a(0)=0.
      do iz=1,(nzp+1)/2-1
        a(iz)=a(iz)*f(iz)
        b(iz)=b(iz)*f(iz)
        a(nzp-iz)=a(nzp-iz)*f(iz)
        b(nzp-iz)=b(nzp-iz)*f(iz)
      enddo
      if(mod(nzp,2).eq.0)a(nzp/2)=a(nzp/2)*(-12.)
      call fft(a,b,1,nzp,1,-1)
      do iz=0,nz-1
        d2zi(ix,iy,iz)=a(iz)
      enddo
   31 continue

c     call free(pa)
c     call free(pb)
c     call free(pf)
      return
      end

      subroutine splineza(imag,nx,ny,nz,d2zi)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 a(0:262143),u(0:511)
c     real*4 a(0:nz-1,0:nz-1),u(0:nz-1)
c     pointer (pa,a),(pu,u)

c     pa=malloc(4*nz*nz)
c     pu=malloc(4*nz)

c     do 1 iz=0,nz-1
c     do 1 jz=0,nz-1
c   1 a(iz,jz)=0.
c     a(0,0)=1.
c     a(nz-1,nz-1)=1.
c     do 2 iz=1,nz-2
c     a(iz,iz-1)=1./6.
c     a(iz,iz+0)=2./3.
c   2 a(iz,iz+1)=1./6.
      a(0)=1.
      a(nz**2-1)=1.
      do 1 iz=1,nz**2-2
    1 a(iz)=0.
      do 2 iz=1,nz-2
      a(iz+nz*(iz-1))=1./6.
      a(iz+nz*(iz+0))=2./3.
    2 a(iz+nz*(iz+1))=1./6.
      call matinv(a,nz,det)

      do 11 ix=0,nx-1
      do 11 iy=0,ny-1
      u(0)=0.
      do iz=1,nz-2
        u(iz)=imag(ix,iy,iz+1)-2.*imag(ix,iy,iz)+imag(ix,iy,iz-1)
      enddo
      u(nz-1)=0.
      do iz=0,nz-1
        d2zi(ix,iy,iz)=0.
        do jz=0,nz-1
c         d2zi(ix,iy,iz)=d2zi(ix,iy,iz)+a(iz,jz)*u(jz)
          d2zi(ix,iy,iz)=d2zi(ix,iy,iz)+a(iz+nz*jz)*u(jz)
        enddo
      enddo
   11 continue

c     call free(pa)
c     call free(pu)
      return
      end

      subroutine splint3dv(imag,nx,ny,nz,d2xi,d2yi,d2zi,x,y,z,v)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yi(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)

      ix=nint(x-.5)
      wx=x-float(ix)
      dowhile(ix.lt.0)
        ix=ix+nx
      enddo
      ix=mod(ix,nx)
      ix1=mod(ix+1,nx)
      ax=1.-wx
      cx=ax*(ax*ax-1.)/6.
      dx=wx*(wx*wx-1.)/6.

      iy=nint(y-.5)
      wy=y-float(iy)
      dowhile(iy.lt.0)
        iy=iy+ny
      enddo
      iy=mod(iy,ny)
      iy1=mod(iy+1,ny)
      ay=1.-wy
      cy=ay*(ay*ay-1.)/6.
      dy=wy*(wy*wy-1.)/6.

      iz=nint(z-.5)
      wz=z-float(iz)
      az=1.-wz
      if(iz.lt.0)then
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        iz=nz-2
        az=0.
        wz=1.
      endif
      cz=az*(az*az-1.)/6.
      dz=wz*(wz*wz-1.)/6.

      v=
     &+ax*(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))

      return
      end

      subroutine splint3dvgh_testh
      parameter (nx=15)
      parameter (ny=9)
      parameter (nz=5)
      real*4 f(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2zf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yf(0:nx-1,0:ny-1,0:nz-1)
      real*4 x(3),g0(3),gx(3),gy(3),gz(3),hess(3,3),test(3,3)
      data del/0.001/

      write (*,"('splint3dvgh_testh')")
      ix0=nx/2
      iy0=ny/2
      iz0=nz/2
      do 2 iz=0,nz-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy,iz)=0.
      f(ix0,iy0,iz0)=1.

      call splinex (f,nx,ny,nz,d2xf)
      call spliney (f,nx,ny,nz,d2yf)
      call splineza(f,nx,ny,nz,d2zf)

      do 11 l=-5,5
      x(1)=float(ix0)+float(l)/10.+0.05
      x(2)=float(iy0)+0.5
      x(3)=float(iz0)+0.5
      call splint3dvgh(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),x(2),x(3),v,g0,hess)

      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1)+del,x(2),    x(3),    v,gx)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2)+del,x(3),    v,gy)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2),    x(3)+del,v,gz)
      do k=1,3
        test(1,k)=(gx(k)-g0(k))/del
        test(2,k)=(gy(k)-g0(k))/del
        test(3,k)=(gz(k)-g0(k))/del
      enddo
      do i=1,3
        write (*,"(3(3f10.6,5x))")(hess(i,j),j=1,3),(test(i,j),j=1,3),(hess(i,j)-test(i,j),j=1,3)
      enddo
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1)-del,x(2),    x(3),    v,gx)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2)-del,x(3),    v,gy)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2),    x(3)-del,v,gz)
      do k=1,3
        test(1,k)=-(gx(k)-g0(k))/del
        test(2,k)=-(gy(k)-g0(k))/del
        test(3,k)=-(gz(k)-g0(k))/del
      enddo
      do i=1,3
        write (*,"(3(3f10.6,5x))")(hess(i,j),j=1,3),(test(i,j),j=1,3),(hess(i,j)-test(i,j),j=1,3)
      enddo
   11 write (*,"()")

      return
      end

      subroutine splint3dvgh_testg
      parameter (nx=15)
      parameter (ny=9)
      parameter (nz=5)
      real*4 f(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2zf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2xf(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yf(0:nx-1,0:ny-1,0:nz-1)
      real*4 x(3),grad(3),hess(3,3),test(3)
      data del/0.001/

      write (*,"('splint3dvgh_testg')")
      ix0=nx/2
      iy0=ny/2
      iz0=nz/2
      do 2 iz=0,nz-1
      do 2 iy=0,ny-1
      do 2 ix=0,nx-1
    2 f(ix,iy,iz)=0.
      f(ix0,iy0,iz0)=1.

      call splinex (f,nx,ny,nz,d2xf)
      call spliney (f,nx,ny,nz,d2yf)
      call splineza(f,nx,ny,nz,d2zf)

      do 11 l=-5,5
      x(1)=float(ix0)+float(l)/10.+0.05
      x(2)=float(iy0)+0.5
      x(3)=float(iz0)+0.5
      call splint3dvgh(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),x(2),x(3),v,grad,hess)

      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1)+del,x(2),    x(3),    gx,grad)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2)+del,x(3),    gy,grad)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2),    x(3)+del,gz,grad)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1)-del,x(2),    x(3),    ux,grad)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2)-del,x(3),    uy,grad)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2),    x(3)-del,uz,grad)
      call splint3dvg(f,nx,ny,nz,d2xf,d2yf,d2zf,x(1),    x(2),    x(3),     v,grad)

      test(1)=0.5*(gx-ux)/del
      test(2)=0.5*(gy-uy)/del
      test(3)=0.5*(gz-uz)/del
   11 write (*,"(10f10.4)")x,v,grad,(grad(i)-test(i),i=1,3)

      return
      end

      subroutine splint3dvg(imag,nx,ny,nz,d2xi,d2yi,d2zi,x,y,z,v,grad)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yi(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 grad(3)

      if(nx.lt.2.or.ny.lt.2.or.nz.lt.2)then
        write (*,"('splint3dvg: image dimensions must be at least 2 in x, y, and z')")
        call exit(-2)
      endif

      ix=nint(x-.5)
      wx=x-float(ix)
      dowhile(ix.lt.0)
        ix=ix+nx
      enddo
      ix=mod(ix,nx)
      ix1=mod(ix+1,nx)
      ax=1.-wx
      cx=ax*(ax*ax-1.)/6.
      dx=wx*(wx*wx-1.)/6.
      ex=(1.-3.*ax*ax)/6.
      fx=(3.*wx*wx-1.)/6.

      iy=nint(y-.5)
      wy=y-float(iy)
      dowhile(iy.lt.0)
        iy=iy+ny
      enddo
      iy=mod(iy,ny)
      iy1=mod(iy+1,ny)
      ay=1.-wy
      cy=ay*(ay*ay-1.)/6.
      dy=wy*(wy*wy-1.)/6.
      ey=(1.-3.*ay*ay)/6.
      fy=(3.*wy*wy-1.)/6.

      iz=nint(z-.5)
      wz=z-float(iz)
      az=1.-wz
      if(iz.lt.0)then
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        iz=nz-2
        az=0.
        wz=1.
      endif
      cz=az*(az*az-1.)/6.
      dz=wz*(wz*wz-1.)/6.
      ez=(1.-3.*az*az)/6.
      fz=(3.*wz*wz-1.)/6.

      v=
     &+ax*(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))

      grad(1)=
     &   -(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &   +(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+ex*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+fx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &   -(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &   +(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &   -(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &   +(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      grad(2)=
     &+ax*(  -(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))   +(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(  -(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))   +(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(  -(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))   +(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(  -(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))   +(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(ey*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+fy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(ey*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+fy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(  -(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))   +(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(  -(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))   +(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      grad(3)=
     &+ax*(ay*(  -imag(ix ,iy,iz)   +imag(ix ,iy,iz+1))+wy*(  -imag(ix ,iy1,iz)   +imag(ix ,iy1,iz+1)))
     &+wx*(ay*(  -imag(ix1,iy,iz)   +imag(ix1,iy,iz+1))+wy*(  -imag(ix1,iy1,iz)   +imag(ix1,iy1,iz+1)))
     &+cx*(ay*(  -d2xi(ix ,iy,iz)   +d2xi(ix ,iy,iz+1))+wy*(  -d2xi(ix ,iy1,iz)   +d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(  -d2xi(ix1,iy,iz)   +d2xi(ix1,iy,iz+1))+wy*(  -d2xi(ix1,iy1,iz)   +d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(  -d2yi(ix ,iy,iz)   +d2yi(ix ,iy,iz+1))+dy*(  -d2yi(ix ,iy1,iz)   +d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(  -d2yi(ix1,iy,iz)   +d2yi(ix1,iy,iz+1))+dy*(  -d2yi(ix1,iy1,iz)   +d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(ez*d2zi(ix ,iy,iz)+fz*d2zi(ix ,iy,iz+1))+wy*(ez*d2zi(ix ,iy1,iz)+fz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(ez*d2zi(ix1,iy,iz)+fz*d2zi(ix1,iy,iz+1))+wy*(ez*d2zi(ix1,iy1,iz)+fz*d2zi(ix1,iy1,iz+1)))

      return
      end

      subroutine splint3dvgh(imag,nx,ny,nz,d2xi,d2yi,d2zi,x,y,z,v,grad,hess)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yi(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 grad(3),hess(3,3)

      if(nx.lt.2.or.ny.lt.2.or.nz.lt.2)then
        write (*,"('splint3dvgh: image dimensions must be at least 2 in x, y, and z')")
        call exit(-2)
      endif

      ix=nint(x-.5)
      wx=x-float(ix)
      dowhile(ix.lt.0)
        ix=ix+nx
      enddo
      ix=mod(ix,nx)
      ix1=mod(ix+1,nx)
      ax=1.-wx
      cx=ax*(ax*ax-1.)/6.
      dx=wx*(wx*wx-1.)/6.
      ex=(1.-3.*ax*ax)/6.
      fx=(3.*wx*wx-1.)/6.

      iy=nint(y-.5)
      wy=y-float(iy)
      dowhile(iy.lt.0)
        iy=iy+ny
      enddo
      iy=mod(iy,ny)
      iy1=mod(iy+1,ny)
      ay=1.-wy
      cy=ay*(ay*ay-1.)/6.
      dy=wy*(wy*wy-1.)/6.
      ey=(1.-3.*ay*ay)/6.
      fy=(3.*wy*wy-1.)/6.

      iz=nint(z-.5)
      wz=z-float(iz)
      az=1.-wz
      if(iz.lt.0)then
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        iz=nz-2
        az=0.
        wz=1.
      endif
      cz=az*(az*az-1.)/6.
      dz=wz*(wz*wz-1.)/6.
      ez=(1.-3.*az*az)/6.
      fz=(3.*wz*wz-1.)/6.

      v=
     &+ax*(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))

      grad(1)=
     &   -(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &   +(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+ex*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+fx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &   -(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &   +(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &   -(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &   +(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      grad(2)=
     &+ax*(  -(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))   +(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(  -(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))   +(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(  -(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))   +(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(  -(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))   +(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(ey*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+fy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(ey*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+fy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(  -(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))   +(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(  -(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))   +(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      grad(3)=
     &+ax*(ay*(  -imag(ix ,iy,iz)   +imag(ix ,iy,iz+1))+wy*(  -imag(ix ,iy1,iz)   +imag(ix ,iy1,iz+1)))
     &+wx*(ay*(  -imag(ix1,iy,iz)   +imag(ix1,iy,iz+1))+wy*(  -imag(ix1,iy1,iz)   +imag(ix1,iy1,iz+1)))
     &+cx*(ay*(  -d2xi(ix ,iy,iz)   +d2xi(ix ,iy,iz+1))+wy*(  -d2xi(ix ,iy1,iz)   +d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(  -d2xi(ix1,iy,iz)   +d2xi(ix1,iy,iz+1))+wy*(  -d2xi(ix1,iy1,iz)   +d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(  -d2yi(ix ,iy,iz)   +d2yi(ix ,iy,iz+1))+dy*(  -d2yi(ix ,iy1,iz)   +d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(  -d2yi(ix1,iy,iz)   +d2yi(ix1,iy,iz+1))+dy*(  -d2yi(ix1,iy1,iz)   +d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(ez*d2zi(ix ,iy,iz)+fz*d2zi(ix ,iy,iz+1))+wy*(ez*d2zi(ix ,iy1,iz)+fz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(ez*d2zi(ix1,iy,iz)+fz*d2zi(ix1,iy,iz+1))+wy*(ez*d2zi(ix1,iy1,iz)+fz*d2zi(ix1,iy1,iz+1)))

      hess(1,1)=
     &+ax*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
      hess(1,2)=
     &   -(  -(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))   +(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &   +(  -(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))   +(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+ex*(  -(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))   +(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+fx*(  -(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))   +(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &   -(ey*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+fy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &   +(ey*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+fy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &   -(  -(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))   +(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &   +(  -(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))   +(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      hess(1,3)=
     &   -(ay*(  -imag(ix ,iy,iz)   +imag(ix ,iy,iz+1))+wy*(  -imag(ix ,iy1,iz)   +imag(ix ,iy1,iz+1)))
     &   +(ay*(  -imag(ix1,iy,iz)   +imag(ix1,iy,iz+1))+wy*(  -imag(ix1,iy1,iz)   +imag(ix1,iy1,iz+1)))
     &+ex*(ay*(  -d2xi(ix ,iy,iz)   +d2xi(ix ,iy,iz+1))+wy*(  -d2xi(ix ,iy1,iz)   +d2xi(ix ,iy1,iz+1)))
     &+fx*(ay*(  -d2xi(ix1,iy,iz)   +d2xi(ix1,iy,iz+1))+wy*(  -d2xi(ix1,iy1,iz)   +d2xi(ix1,iy1,iz+1)))
     &   -(cy*(  -d2yi(ix ,iy,iz)   +d2yi(ix ,iy,iz+1))+dy*(  -d2yi(ix ,iy1,iz)   +d2yi(ix ,iy1,iz+1)))
     &   +(cy*(  -d2yi(ix1,iy,iz)   +d2yi(ix1,iy,iz+1))+dy*(  -d2yi(ix1,iy1,iz)   +d2yi(ix1,iy1,iz+1)))
     &   -(ay*(ez*d2zi(ix ,iy,iz)+fz*d2zi(ix ,iy,iz+1))+wy*(ez*d2zi(ix ,iy1,iz)+fz*d2zi(ix ,iy1,iz+1)))
     &   +(ay*(ez*d2zi(ix1,iy,iz)+fz*d2zi(ix1,iy,iz+1))+wy*(ez*d2zi(ix1,iy1,iz)+fz*d2zi(ix1,iy1,iz+1)))
      hess(2,1)=hess(1,2)
      hess(2,2)=
     &+ax*(ay*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+wy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+wy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
      hess(2,3)=
     &+ax*(  -(  -imag(ix ,iy,iz)   +imag(ix ,iy,iz+1))   +(  -imag(ix ,iy1,iz)   +imag(ix ,iy1,iz+1)))
     &+wx*(  -(  -imag(ix1,iy,iz)   +imag(ix1,iy,iz+1))   +(  -imag(ix1,iy1,iz)   +imag(ix1,iy1,iz+1)))
     &+cx*(  -(  -d2xi(ix ,iy,iz)   +d2xi(ix ,iy,iz+1))   +(  -d2xi(ix ,iy1,iz)   +d2xi(ix ,iy1,iz+1)))
     &+dx*(  -(  -d2xi(ix1,iy,iz)   +d2xi(ix1,iy,iz+1))   +(  -d2xi(ix1,iy1,iz)   +d2xi(ix1,iy1,iz+1)))
     &+ax*(ey*(  -d2yi(ix ,iy,iz)   +d2yi(ix ,iy,iz+1))+fy*(  -d2yi(ix ,iy1,iz)   +d2yi(ix ,iy1,iz+1)))
     &+wx*(ey*(  -d2yi(ix1,iy,iz)   +d2yi(ix1,iy,iz+1))+fy*(  -d2yi(ix1,iy1,iz)   +d2yi(ix1,iy1,iz+1)))
     &+ax*(  -(ez*d2zi(ix ,iy,iz)+fz*d2zi(ix ,iy,iz+1))   +(ez*d2zi(ix ,iy1,iz)+fz*d2zi(ix ,iy1,iz+1)))
     &+wx*(  -(ez*d2zi(ix1,iy,iz)+fz*d2zi(ix1,iy,iz+1))   +(ez*d2zi(ix1,iy1,iz)+fz*d2zi(ix1,iy1,iz+1)))
      hess(3,1)=hess(1,3)
      hess(3,2)=hess(2,3)
      hess(3,3)=
     &+ax*(ay*(az*d2zi(ix ,iy,iz)+wz*d2zi(ix ,iy,iz+1))+wy*(az*d2zi(ix ,iy1,iz)+wz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2zi(ix1,iy,iz)+wz*d2zi(ix1,iy,iz+1))+wy*(az*d2zi(ix1,iy1,iz)+wz*d2zi(ix1,iy1,iz+1)))

      return
      end

      subroutine splint3dvgl(imag,nx,ny,nz,d2xi,d2yi,d2zi,x,y,z,v,grad,del2)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yi(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)
      real*4 grad(3)

      if(nx.lt.2.or.ny.lt.2.or.nz.lt.2)then
        write (*,"('splint3dvgh: image dimensions must be at least 2 in x, y, and z')")
        call exit(-2)
      endif

      ix=nint(x-.5)
      wx=x-float(ix)
      dowhile(ix.lt.0)
        ix=ix+nx
      enddo
      ix=mod(ix,nx)
      ix1=mod(ix+1,nx)
      ax=1.-wx
      cx=ax*(ax*ax-1.)/6.
      dx=wx*(wx*wx-1.)/6.
      ex=(1.-3.*ax*ax)/6.
      fx=(3.*wx*wx-1.)/6.

      iy=nint(y-.5)
      wy=y-float(iy)
      dowhile(iy.lt.0)
        iy=iy+ny
      enddo
      iy=mod(iy,ny)
      iy1=mod(iy+1,ny)
      ay=1.-wy
      cy=ay*(ay*ay-1.)/6.
      dy=wy*(wy*wy-1.)/6.
      ey=(1.-3.*ay*ay)/6.
      fy=(3.*wy*wy-1.)/6.

      iz=nint(z-.5)
      wz=z-float(iz)
      az=1.-wz
      if(iz.lt.0)then
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        iz=nz-2
        az=0.
        wz=1.
      endif
      cz=az*(az*az-1.)/6.
      dz=wz*(wz*wz-1.)/6.
      ez=(1.-3.*az*az)/6.
      fz=(3.*wz*wz-1.)/6.

      v=
     &+ax*(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))

      grad(1)=
     &   -(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &   +(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+ex*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+fx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &   -(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &   +(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &   -(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &   +(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      grad(2)=
     &+ax*(  -(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))   +(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(  -(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))   +(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(  -(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))   +(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(  -(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))   +(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(ey*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+fy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(ey*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+fy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(  -(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))   +(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(  -(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))   +(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))
      grad(3)=
     &+ax*(ay*(  -imag(ix ,iy,iz)   +imag(ix ,iy,iz+1))+wy*(  -imag(ix ,iy1,iz)   +imag(ix ,iy1,iz+1)))
     &+wx*(ay*(  -imag(ix1,iy,iz)   +imag(ix1,iy,iz+1))+wy*(  -imag(ix1,iy1,iz)   +imag(ix1,iy1,iz+1)))
     &+cx*(ay*(  -d2xi(ix ,iy,iz)   +d2xi(ix ,iy,iz+1))+wy*(  -d2xi(ix ,iy1,iz)   +d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(  -d2xi(ix1,iy,iz)   +d2xi(ix1,iy,iz+1))+wy*(  -d2xi(ix1,iy1,iz)   +d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(  -d2yi(ix ,iy,iz)   +d2yi(ix ,iy,iz+1))+dy*(  -d2yi(ix ,iy1,iz)   +d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(  -d2yi(ix1,iy,iz)   +d2yi(ix1,iy,iz+1))+dy*(  -d2yi(ix1,iy1,iz)   +d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(ez*d2zi(ix ,iy,iz)+fz*d2zi(ix ,iy,iz+1))+wy*(ez*d2zi(ix ,iy1,iz)+fz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(ez*d2zi(ix1,iy,iz)+fz*d2zi(ix1,iy,iz+1))+wy*(ez*d2zi(ix1,iy1,iz)+fz*d2zi(ix1,iy1,iz+1)))

      del2=
     &+ax*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(ay*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+wy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+wy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(az*d2zi(ix ,iy,iz)+wz*d2zi(ix ,iy,iz+1))+wy*(az*d2zi(ix ,iy1,iz)+wz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2zi(ix1,iy,iz)+wz*d2zi(ix1,iy,iz+1))+wy*(az*d2zi(ix1,iy1,iz)+wz*d2zi(ix1,iy1,iz+1)))

      return
      end

      subroutine splint3dvl(imag,nx,ny,nz,d2xi,d2yi,d2zi,x,y,z,v,del2)
      real*4 imag(0:nx-1,0:ny-1,0:nz-1),d2xi(0:nx-1,0:ny-1,0:nz-1)
      real*4 d2yi(0:nx-1,0:ny-1,0:nz-1),d2zi(0:nx-1,0:ny-1,0:nz-1)

      if(nx.lt.2.or.ny.lt.2.or.nz.lt.2)then
        write (*,"('splint3dvgh: image dimensions must be at least 2 in x, y, and z')")
        call exit(-2)
      endif

      ix=nint(x-.5)
      wx=x-float(ix)
      dowhile(ix.lt.0)
        ix=ix+nx
      enddo
      ix=mod(ix,nx)
      ix1=mod(ix+1,nx)
      ax=1.-wx
      cx=ax*(ax*ax-1.)/6.
      dx=wx*(wx*wx-1.)/6.

      iy=nint(y-.5)
      wy=y-float(iy)
      dowhile(iy.lt.0)
        iy=iy+ny
      enddo
      iy=mod(iy,ny)
      iy1=mod(iy+1,ny)
      ay=1.-wy
      cy=ay*(ay*ay-1.)/6.
      dy=wy*(wy*wy-1.)/6.

      iz=nint(z-.5)
      wz=z-float(iz)
      az=1.-wz
      if(iz.lt.0)then
        iz=0
        az=1.
        wz=0.
      endif
      if(iz.ge.nz-1)then
        iz=nz-2
        az=0.
        wz=1.
      endif
      cz=az*(az*az-1.)/6.
      dz=wz*(wz*wz-1.)/6.

      v=
     &+ax*(ay*(az*imag(ix ,iy,iz)+wz*imag(ix ,iy,iz+1))+wy*(az*imag(ix ,iy1,iz)+wz*imag(ix ,iy1,iz+1)))
     &+wx*(ay*(az*imag(ix1,iy,iz)+wz*imag(ix1,iy,iz+1))+wy*(az*imag(ix1,iy1,iz)+wz*imag(ix1,iy1,iz+1)))
     &+cx*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+dx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(cy*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+dy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(cy*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+dy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(cz*d2zi(ix ,iy,iz)+dz*d2zi(ix ,iy,iz+1))+wy*(cz*d2zi(ix ,iy1,iz)+dz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(cz*d2zi(ix1,iy,iz)+dz*d2zi(ix1,iy,iz+1))+wy*(cz*d2zi(ix1,iy1,iz)+dz*d2zi(ix1,iy1,iz+1)))

      del2=
     &+ax*(ay*(az*d2xi(ix ,iy,iz)+wz*d2xi(ix ,iy,iz+1))+wy*(az*d2xi(ix ,iy1,iz)+wz*d2xi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2xi(ix1,iy,iz)+wz*d2xi(ix1,iy,iz+1))+wy*(az*d2xi(ix1,iy1,iz)+wz*d2xi(ix1,iy1,iz+1)))
     &+ax*(ay*(az*d2yi(ix ,iy,iz)+wz*d2yi(ix ,iy,iz+1))+wy*(az*d2yi(ix ,iy1,iz)+wz*d2yi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2yi(ix1,iy,iz)+wz*d2yi(ix1,iy,iz+1))+wy*(az*d2yi(ix1,iy1,iz)+wz*d2yi(ix1,iy1,iz+1)))
     &+ax*(ay*(az*d2zi(ix ,iy,iz)+wz*d2zi(ix ,iy,iz+1))+wy*(az*d2zi(ix ,iy1,iz)+wz*d2zi(ix ,iy1,iz+1)))
     &+wx*(ay*(az*d2zi(ix1,iy,iz)+wz*d2zi(ix1,iy,iz+1))+wy*(az*d2zi(ix1,iy1,iz)+wz*d2zi(ix1,iy1,iz+1)))

      return
      end

