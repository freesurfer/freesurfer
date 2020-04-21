cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  2007/04/25  05:15:34  avi
c gcc v3 compliant
c
c Revision 1.1  1998/10/02  20:23:27  mcavoy
c Initial revision
c
c Revision 1.1  1997/11/14  17:16:24  tscull
c Initial revision
c
c Revision 1.3  1997/04/15  08:21:03  avi
c trap input image NaNs in subroutine imgvalx
c
c Revision 1.2  1996/12/28  01:40:42  avi
c minor changes - don't evaluate nonexistant memory if nz.eq.1
c
c Revision 1.1  1996/12/27  20:48:20  avi
c Initial revision
c
c revision 1.2  1993/10/25  17:37:33  avi

      subroutine imgvalx(imgt,nx,ny,nz,center,mmppix,x,v,lslice)
c     avi snyder 08-06-94
c     returns in v the interpolated value of imgt at locus x.
c     lslice returns the slice from which most of the data is taken.
c     lslice = 0 (and v = 0.) if imgt is undefined at locus x.

      real*4 imgt(nx,ny,nz)
      real*4 center(3),mmppix(3),x(3)
      logical*4 ldefined,lok
      real*4 wx(3),xf(3)
      integer*4 ix(3),mx(3)

      mx(1)=nx
      mx(2)=ny
      mx(3)=nz

      ldefined=.true.
      do 1 k=1,3
      center(k)=sign(center(k),mmppix(k))
      xf(k)=(center(k)+x(k))/mmppix(k)
      ix(k)=nint(xf(k)-0.5)
      wx(k)=xf(k)-float(ix(k))
      if((ix(k).eq.0.and.wx(k).gt.0.999).or.mx(k).eq.1)then
        ix(k)=1
        wx(k)=0.
      endif
      if(mx(k).gt.1.and.ix(k).eq.mx(k).and.wx(k).lt.0.001)then
        ix(k)=mx(k)-1
        wx(k)=1.
      endif
    1 ldefined=ldefined.and.((mx(k).eq.1).or.(ix(k).gt.0).and.(ix(k).lt.mx(k)))

      if(.not.ldefined)then
        v=0.
        lslice=0
        return
      endif

      lok=(imgt(ix(1)+0,ix(2)+0,ix(3)+0).ge.0.0.or.imgt(ix(1)+0,ix(2)+0,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+0,ix(3)+0).ge.0.0.or.imgt(ix(1)+1,ix(2)+0,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+0,ix(2)+1,ix(3)+0).ge.0.0.or.imgt(ix(1)+0,ix(2)+1,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+1,ix(3)+0).ge.0.0.or.imgt(ix(1)+1,ix(2)+1,ix(3)+0).le.0.0).and.
     &    (imgt(ix(1)+0,ix(2)+0,ix(3)+1).ge.0.0.or.imgt(ix(1)+0,ix(2)+0,ix(3)+1).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+0,ix(3)+1).ge.0.0.or.imgt(ix(1)+1,ix(2)+0,ix(3)+1).le.0.0).and.
     &    (imgt(ix(1)+0,ix(2)+1,ix(3)+1).ge.0.0.or.imgt(ix(1)+0,ix(2)+1,ix(3)+1).le.0.0).and.
     &    (imgt(ix(1)+1,ix(2)+1,ix(3)+1).ge.0.0.or.imgt(ix(1)+1,ix(2)+1,ix(3)+1).le.0.0)
      if(.not.lok)then
        v=0.
        lslice=0
        return
      endif

      v=
     & +(1.-wx(3))*((1.-wx(2))*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+0,ix(3)+0)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+0,ix(3)+0))
     &               +wx(2)*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+1,ix(3)+0)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+1,ix(3)+0)))
       if(mx(3).gt.1)v=v
     &      +wx(3)*((1.-wx(2))*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+0,ix(3)+1)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+0,ix(3)+1))
     &               +wx(2)*(
     &               (1.-wx(1))*imgt(ix(1)+0,ix(2)+1,ix(3)+1)
     &                  +wx(1) *imgt(ix(1)+1,ix(2)+1,ix(3)+1)))
      lslice=ix(3)
      if(wx(3).gt.0.5)lslice=lslice+1
      return
      end

      subroutine peak_find(mode,imgt,nx,ny,nz,mmppix,centert,coord)
c     iabs(mode) optimize
c     1          x
c     2          y
c     3          x y
c     4          z
c     5          x z
c     6          y z
c     7          x y z
c     sign of extremum found = sign of mode
      real*4 imgt(nx,ny,nz)
      real*4 mmppix(3),centert(3),coord(3)
      parameter (nterm=3)
      real*4 taram(nterm),curva(nterm)
      real*4 daram(nterm)/3*2./
      real*4 coef(0:2),array(3,3)
      parameter (ni=2)
      parameter (npts=2*ni+1)
      real*4 x(-ni:ni),y(-ni:ni)
      logical*4 liter,lfar

      call imgvalx(imgt,nx,ny,nz,centert,mmppix,coord,v,lslice)
      if(lslice.le.0)goto 99
      write(*,102)coord,v
  102 format(3f7.2,f10.1)
      sense=sign(1.,float(mode))
      jmask=iabs(mode)

      nnn=0
      niter=5
   83 liter=.false.
      lfar=.false.
      do 81 j=1,3
      if(iand(jmask,2**(j-1)).eq.0)goto 81
      taram(j)=coord(j)
      do 84 i=-ni,ni
      coord(j)=taram(j)+daram(j)*float(i)/float(ni)
      call imgvalx(imgt,nx,ny,nz,centert,mmppix,coord,v,lslice)
      if(lslice.le.0)goto 99
      x(i)=float(i)/float(ni)
   84 y(i)=-sense*v
      call polfit(x,y,npts,coef,array,3,chisqr)
      t=.5*coef(1)/coef(2)
      if(coef(2).lt.0..or.abs(t).gt..5)then
        t0=t
        t=sign(.5,coef(1))
c       write(*,"(' T ',f10.6,' ->',f10.6)")t0,t
        lfar=.true.
      endif
      liter=liter.or.abs(t).gt.0.02
      curva(j)=coef(2)
      coord(j)=taram(j)-daram(j)*t
   81 continue
      nnn=nnn+1
      if(nnn.gt.20)then
        write(*,"(' failed convergence')")
        goto 89
      endif
      if(.not.lfar)then
        niter=niter-1
      endif
      if(liter.and.niter.gt.0)goto 83

      call imgvalx(imgt,nx,ny,nz,centert,mmppix,coord,v,lslice)
      if(lslice.le.0)goto 99
      write(*,102)coord,v
   89 return
   99 write (*,"('peak_find: image boundary reached')")
      return
      end

      subroutine edge_trace(imgt,nx,ny,nz,mmppix,center,param,edge,npmax,np)
      real*4 imgt(nx,ny,nz)
      real*4 mmppix(3),center(3),edge(3,npmax)
      real*4 param(16)
c     param(01:03)	normal to osculating plane
c     param(04)		del arclength
c     param(05)		sigma eval
      real*4 n1(3)/0.,0.,1./		! for initial t, n guess
      real*4 n2(3)/0.,1.,0./		! for initial t, n guess
      real*4 t(3),n(3),b(3),x(3)
      real*4 coef(0:2),array(3,3)
      parameter (nk=3)
      parameter (npts=2*nk+1)
      real*4 u(-nk:nk),v(-nk-1:nk+1),dv(-nk:nk)
      parameter(nr=32)
      logical*4 ldebug/.false./

      pi=4.*atan(1.)
      a=sqrt(dot(param,param))		! compute binormal unit vector
      do i=1,3
        b(i)=param(i)/a
      enddo

c     initial guess for t, n
      if(abs(dot(b,n1)).lt.abs(dot(b,n2)))then
        do i=1,3
          n(i)=n1(i)
        enddo
      else
        do i=1,3
          n(i)=n2(i)
        enddo
      endif
      call cross(n,b,t)     
     
      np=1
   11 s=0.
      c=0.
      if(ldebug)write(*, "('t ',3f10.4)")t
      if(ldebug)write(*, "('n ',3f10.4)")n
      if(ldebug)write(*, "('np, edge(np)',i6,3f7.2)")np,(edge(i,np),i=1,3)
      do 21 k=0,nr-1
      a=2.*pi*float(k)/float(nr)
      do i=1,3
        x(i)=edge(i,np)+param(05)*(cos(a)*n(i)+sin(a)*t(i))
      enddo
      call imgvalx(imgt,nx,ny,nz,center,mmppix,x,g,lslice)
      if(lslice.le.0)goto 99
      s=s+sin(a)*g
      c=c+cos(a)*g
c     if(ldebug)write(*, "(i6,3f7.2,f7.0)")k,x,g
   21 continue
      theta=atan2(s,c)
      if(ldebug)write(*, "('s,c,atan2(s,c) (deg)',3f10.4)")s,c,theta*45./atan(1.)
      do i=1,3
        n(i)=sin(theta)*t(i)+cos(theta)*n(i)
      enddo
      call cross(n,b,t)
      if(ldebug)write(*, "('t ',3f10.4)")t
      if(ldebug)write(*, "('n ',3f10.4)")n

      iter=0
   10 iter=iter+1
      if(ldebug)write(*, "('iter',i4)")iter
      do 9 k=-nk-1,nk+1			! sample image along normal
      do i=1,3
        x(i)=edge(i,np)+param(05)*(float(k)/float(nk))*n(i)
      enddo
      call imgvalx(imgt,nx,ny,nz,center,mmppix,x,v(k),lslice)
      if(lslice.le.0)goto 99
c     if(ldebug)write(*, "('imgvalx x,v',3f7.2,f9.1)")x,v(k)
    9 continue
      do 8 k=-nk,nk			! differentiate
      u(k)=float(k)/float(nk)
    8 dv(k)=v(k+1)-v(k-1)
      if(ldebug)write(*, "('dv',10f9.1)")dv

      call polfit(u,dv,npts,coef,array,3,chisqr)
      if(ldebug)write(*, "('coef ',3g12.4)")coef
      q=.5*coef(1)/coef(2)		! find gradient extremum
      if(abs(q).gt.1.5.and.np.gt.1)then
        write(*, "('edge_trace: edge end encountered')")
        goto 98
      endif
      if(abs(q).gt.0.5)q=sign(0.5,q)
      do i=1,3				! move current edge point onto gradient extremum
        edge(i,np)=edge(i,np)-param(5)*q*n(i)
      enddo
      if(abs(q).gt.0.05.and.iter.lt.20)goto 10
      write(*, "('np, edge(np)',i6,3f7.2)")np,(edge(i,np),i=1,3)

      if(np.ge.npmax)return

      do j=1,np-2			! edge closed?
        do i=1,3
          x(i)=edge(i,np)-edge(i,j)
        enddo
        a=sqrt(dot(x,x))
        if(a.lt.param(4))then
          write(*, "('edge_trace: edge closed')")
          return
        endif
      enddo

      do i=1,3				! advance trace
        edge(i,np+1)=edge(i,np)+param(04)*t(i)
      enddo
      np=np+1
      goto 11

   99 write(*, "('edge_trace: image boundary reached')")
   98 np=np-1
      return
      end

      subroutine cross(a,b,c)
      real*4 a(3),b(3),c(3)
      c(1)= a(2)*b(3)-a(3)*b(2)
      c(2)=-a(1)*b(3)+a(3)*b(1)
      c(3)= a(1)*b(2)-a(2)*b(1)
      return
      end

      function dot(a,b)
      real*4 a(3),b(3)
      dot=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end

      function dist(a,b)
      real*4 a(3),b(3)
      dist=sqrt((a(1)-b(1))**2+(a(2)-b(2))**2+(a(3)-b(3))**2)
      return
      end

      subroutine reflect(a,b)
      real*4 a(3),b(3)
      b(1)=-a(1)
      b(2)=-a(2)
      b(3)=-a(3)
      return
      end
