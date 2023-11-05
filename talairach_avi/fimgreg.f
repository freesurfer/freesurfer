c$Header: /space/repo/1/dev/dev/talairach_avi/fimgreg.f,v 1.1 2007/05/04 22:33:59 nicks Exp $
c$Log: fimgreg.f,v $
cRevision 1.1  2007/05/04 22:33:59  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.19  2007/04/25  05:16:26  avi
c gcc v3 compliant
c
c Revision 1.18  2006/02/10  05:41:50  avi
c double index in mask memory (to 2097152) in imgrege()
c
c Revision 1.17  1999/01/25  03:45:41  avi
c Hessian mode iteration parameter change only on improvement
c
c Revision 1.16  1998/12/11  09:06:33  avi
c Revision 1.15  1998/12/02  07:50:11  avi
c protect Hessian algorithm post cndnum < 0 eta < 0.95*eta0
c
c Revision 1.14  1998/11/01  01:02:46  avi
c new eta,q format
c
c Revision 1.13  1998/08/20  06:52:29  avi
c retune superfine mode for baboon-baboon operation
c
c Revision 1.12  1997/10/05  21:46:31  avi
c call exit(-2) on termoinal condition no voxels in register
c
c Revision 1.11  1997/10/05  21:42:57  avi
c escale 80 -> 100 and 8 -> 20
c
c Revision 1.10  1997/06/11  23:25:02  avi
c in hessian mode make squeeze factor proportional to log(jmax)
c
c Revision 1.9  1997/06/06  03:16:16  avi
c new hessian parameter sample radius strategy
c
c Revision 1.8  1997/06/03  00:22:17  avi
c increase index in mask memory 500000 to 1048576
c
c Revision 1.7  1997/03/23  07:05:08  avi
c in hessian mode prevent additional parameter search radius constriction once
c indefinite condtion encountered
c
c Revision 1.6  1997/01/16  00:30:09  avi
c correct del=0. on first call to fimgrege
c
c Revision 1.5  1997/01/15  09:28:06  avi
c increase niter from 3 to 4 in hessian method
c define new mode bit lsuperfine=mode.and.512.ne.0
c
c Revision 1.4  1997/01/02  07:57:16  avi
c hessian mode convergence failure escape
c
c Revision 1.3  1997/01/01  23:24:41  avi
c improve hessian mode algorithm convergence
c
c Revision 1.2  1997/01/01  02:17:56  avi
c add rscale to hessian algorithm
c
      subroutine fimgreg_rcs
      write(*,"('fimgreg.f')")
      return
      end

      subroutine fimgreg(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &                   img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode)
c     finds optimal translation/rotation to achieve registration of img2 on img1
      real*4 img1(nx1,ny1,nz1)
      integer*2 msk1(nx1,ny1,nz1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(nx2,ny2,nz2)
      integer*2 msk2(nx2,ny2,nz2)
      real*4 mmppix2(3),center2(3)
      real*4 param(13)
      parameter (dd=5.0)		! linear distance search increment in units of mmppix
      parameter (da=0.0872665)		! angle in radians = +/- 5 deg search range
      parameter (ds=0.06)		! stretch increment
      real*4 raram(12,3)/dd,dd,dd,da,da,da,ds,ds,ds,0.,0.,0.,	! range
     &                   dd,dd,da,ds,ds,0.,0.,0.,0.,0.,0.,0.,
     &                   2.,2.,2.,ds,ds,ds,ds,ds,ds,ds,ds,ds/
      parameter (sr=57.29578)		! degrees/radian
      parameter (ss=100.)		! converts stretch to percent
      real*4 saram(12,3)/1.,1.,1.,sr,sr,sr,ss,ss,ss,0.,0.,0.,	! scale
     &                   1.,1.,sr,ss,ss,0.,0.,0.,0.,0.,0.,0.,
     &                   1.,1.,1.,ss,ss,ss,ss,ss,ss,ss,ss,ss/
      real*4 qaram(12)						! search range
      real*4 garam(12)						! computed gradient
      real*4 baram(12)						! best parameters
      real*4 hessian(144),w(144),hesstmp(144)
      parameter (nterm=12)
      real*4 taram(nterm),curva(nterm)
      real*4 coef(0:2),array(3,3)
      parameter (ni=3)
      parameter (npts=2*ni+1)
      real*4 x(-ni:ni),y(-ni:ni)
      real*4 report(-ni:ni,12)
      logical*4 liter,lhessian
      logical*4 lenable,l3d,lstretch,lvoxsiz,lsuperfine,lfast,lfine,lfind,lcorrel
      logical*4 lrestart/.false./
      logical*4 ldebug/.false./

      if(ldebug)then
        write(*, "(6i6)")nx1,ny1,nz1,nx2,ny2,nz2
        write(*, "(6f10.4)")mmppix1,center1
        write(*, "(6f10.4)")mmppix2,center2
      endif
      if(nz2.lt.2)mode=iand(mode,not(2))	! prevent attempt to 3d align
      lenable=   iand(mode,1).ne.0
      l3d=       iand(mode,2).ne.0
      lstretch=  iand(mode,4).ne.0
      lvoxsiz=   iand(mode,8).ne.0
      lcorrel=   iand(mode,256).ne.0
      lsuperfine=iand(mode,512).ne.0
      lfast=     iand(mode,1024).ne.0
      lfine=     iand(mode,2048).ne.0
      lfind=     iand(mode,4096).ne.0
      lhessian=  iand(mode,8192).ne.0

      call fimgreg_rcs
      write(*,"('image alignment mode ',i6,' decimal ',z8,' hex')")mode,mode

      if(lenable)then
        if(lstretch)then
          jndex=3
          if(l3d)then
            jmax=12
          else
            jmax=6
          endif
        else			! .not.lstretch
          if(l3d)then
            jndex=1
            jmax=6
            if(lvoxsiz)then
              mskbit=16
              do i=1,3
                if(iand(mode,mskbit).eq.0)jmax=jmax+1
                mskbit=mskbit*2
              enddo
            endif
          else			! 2d
            jndex=2
            jmax=3
            if(lvoxsiz)then
              mskbit=16
              do i=1,2
                if(iand(mode,mskbit).eq.0)jmax=jmax+1
                mskbit=mskbit*2
              enddo
            endif
          endif
        endif
      else
        jmax=0
      endif

      write(*,"('parameters')")
      write(*,102)(param(k)*saram(k,jndex),k=1,jmax)
      if(lhessian)goto 70		! hessian branch point
      del=5.0				! 3d grid sampling interval in mm
      if(lsuperfine)then
        del=2.5
        rscale=0.5
        do j=1,3
          raram(j,jndex)=1.
        enddo
      endif
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,del,eta,q)
      write(*,"('eta,q',f8.5,f9.5)")eta,q
      niter=5
      nnn=0
      do 82 j=1,jmax
   82 qaram(j)=1.
   83 liter=.false.
      if(lfast)then
        rscale=1.2
        del=12.0			! 3d grid sampling interval in mm
      elseif(.not.lsuperfine)then
        rscale=0.1*float(niter+3)
        del=4.+float(niter)		! 3d grid sampling interval in mm
      endif
      if(lfine)del=5.0
      mod1=ior(mode,z'20000')		! tell imgrege to set index in mask memory
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mod1,del,eta,q)
      write(*,"('niter,ni,del,rscale ',2i5,2f10.4)")niter,ni,del,rscale
      do 81 j=1,jmax
      taram(j)=param(j)
      do 84 i=-ni,ni
      param(j)=taram(j)+rscale*qaram(j)*raram(j,jndex)*float(i)/float(ni)
      mod1=ior(mode,z'10000')		! tell imgrege to use index in mask memory
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mod1,del,eta,q)
      x(i)=float(i)/float(ni)
      y(i)=1.-eta
      if(ldebug)then
        write(*,"('parameters')")
        write(*,"(12f8.4)")(param(k)*saram(k,jndex),k=1,jmax)
        write(*,"('eta,q',f8.5,f9.5)")eta,q
      endif
   84 report(i,j)=y(i)
      call polfit(x,y,npts,coef,array,3,chisqr)
      t=.5*coef(1)/coef(2)
      if(coef(2).lt.0..or.abs(t).gt..5)then
        call polfit(x,y,npts,coef,array,2,chisqr)
        t0=t
        t=sign(.5,coef(1))
        niter=min0(niter+1,8)
        qaram(j)=qaram(j)*1.15
        write(*,"('parameter',i3,f10.6,' ->',f10.6,' qaram',f10.4)")j,t0,t,qaram(j)
      endif
      liter=liter.or.abs(t).gt.0.01
      curva(j)=coef(2)
   81 param(j)=taram(j)-rscale*qaram(j)*raram(j,jndex)*t
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,del,eta,q)
      write(*,"('parameters')")
      write(*,102)(param(k)*saram(k,jndex),k=1,jmax)
      write(*,"('eta,q',f8.5,f9.5)")eta,q
      param(13)=q
  102 format(6f10.4)
      nnn=nnn+1
      if(nnn.gt.40)then
        write(*,"('fimgreg: failed convergence mode ',z8,' hex')")mode
        goto 89
      endif
      niter=niter-1
      if(.true.)then
        do i=-ni,ni
          write(*,"(12f8.4)")(report(i,j),j=1,jmax)
        enddo
      endif
      if(liter.and.niter.gt.0)goto 83
   89 write(*,"('registration optimization eta report')")
      write(*,"('parameter search radius')")
      write(*,102)(rscale*qaram(j)*raram(j,jndex)*saram(j,jndex),j=1,jmax)
      write(*,"('100000*second partial in parameter space')")
      write(*,"(6f10.0)")(100000.*curva(j)/(rscale*qaram(j)*raram(j,jndex)*saram(j,jndex))**2,j=1,jmax)
      do 88 i=-ni,ni
   88 write(*,"(12f8.4)")(report(i,j),j=1,jmax)
      write(*,102)(param(j)*saram(j,jndex),j=1,jmax)
      write(*,"('eta,q',f8.5,f9.5)")eta,q
      return

   70 if(lcorrel)then
        escale=100.
      else
        escale=20.
      endif
      do j=1,jmax
        qaram(j)=raram(j,jndex)
        baram(j)=param(j)
      enddo
      del=7.5
      if(lfine)del=5.0
      if(lsuperfine)del=2.5
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,del,eta,q)
      write(*,"('eta,q',f8.5,f9.5)")eta,q

      etamax=eta
      rscale=1.
      if(lsuperfine)then
        rscale=0.5*rscale
        escale=0.5*escale
      endif
      nnn=0
      niter=4
   79 err=1.-eta
   78 if(lrestart)then
        do j=1,jmax
          param(j)=baram(j)
        enddo
        lrestart=.false.
        write(*,"('restart rscale ',f7.2,' ->',f7.2)")rscale,rscale*1.5
        rscale=rscale*1.5
      endif
      write(*,"('niter,del,escale ',i5,f10.4,f10.2)")niter,del,escale
      write(*,"('parameter search radius')")
      write(*,102)(rscale*qaram(j)*saram(j,jndex),j=1,jmax)

      do 71 j=1,jmax
      do k=1,jmax
        taram(k)=param(k)
      enddo
      g=0.
      h=-2.*err
      do 72 i=-1,1,2
      taram(j)=param(j)+float(i)*qaram(j)*rscale
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,taram,mode,del,eta,q)
      g=g+float(i)*(1.-eta)
   72 h=h+(1.-eta)
      garam(j)=0.5*g*escale
   71 call array_w(hessian,jmax,j,j,h*escale)

      do 73 j=1,jmax-1
      do 73 jj=j+1,jmax
      do k=1,jmax
        taram(k)=param(k)
      enddo
      h=0.
      do 75 i=-1,1,2
      do 75 ii=-1,1,2
      taram(j) =param(j) +float(i) *qaram(j) *rscale
      taram(jj)=param(jj)+float(ii)*qaram(jj)*rscale
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,taram,mode,del,eta,q)
   75 h=h+float(ii*i)*(1.-eta)
      call array_w(hessian,jmax,j,jj,0.25*h*escale)
   73 call array_w(hessian,jmax,jj,j,0.25*h*escale)

      if(ldebug)then
        write(*, "('linear system')")
        do i=1,jmax
          do j=1,jmax
            call array_r(hessian,jmax,i,j,taram(j))
          enddo
          write(*,"(f8.4,2x,12f8.4)")garam(i),(taram(j),j=1,jmax)
        enddo
      endif

      call matcop(hessian,hesstmp,jmax)
      call eigen(hesstmp,w,jmax)
      cndnum=hesstmp(1)/hesstmp(jmax**2)
      call matcop(hessian,hesstmp,jmax)
      call matinv(hessian,jmax,det)
      write(*,"('hessian determinant ',e10.4,'  condition number ',e10.4)")det,cndnum
      if(cndnum.lt.0.)then
        nnn=nnn+1
	if(nnn.gt.2)then
          write(*,"('fimgreg: failed convergence mode ',z8,' hex')")mode
          goto 76
        endif
        lrestart=.true.
        goto 78
      endif

      do 77 j=1,jmax
      do 77 i=1,jmax
      call array_r(hessian,jmax,j,i,h)
   77 param(j)=param(j)-h*garam(i)*qaram(j)*rscale
      write(*,"('parameters')")
      write(*,102)(param(j)*saram(j,jndex),j=1,jmax)
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,del,eta,q)
      write(*,"('eta,q',f8.5,f9.5)")eta,q
      if(eta.lt.0.95*etamax)then
        lrestart=.true.
        goto 78
      endif
      if(eta.gt.etamax)then
        do j=1,jmax
          baram(j)=param(j)
        enddo
        etamax=eta
      else
        do j=1,jmax
          param(j)=baram(j)
        enddo
        eta=etamax
      endif

c     if(nnn.eq.0)then
        squeeze=alog10(float(jmax-1)*15./cndnum)
        squeeze=amax1(squeeze,0.5)
        squeeze=amin1(squeeze,2.0)
        do j=1,jmax
          call array_r(hesstmp,jmax,j,j,t)
          qaram(j)=qaram(j)/sqrt(squeeze*t)
        enddo
        escale=escale*squeeze
c     endif
      niter=niter-1
      if(niter.gt.0)goto 79      

   76 if(.not.lsuperfine)then
        do j=1,jmax
          param(j)=baram(j)
        enddo
      endif
      call imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &             img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,del,eta,q)
      param(13)=q
      write(*,"('parameters')")
      write(*,102)(param(j)*saram(j,jndex),j=1,jmax)
      write(*,"('eta,q',f8.5,f9.5)")eta,q

      return
      end

      subroutine array_w(a,n,i,j,t)
      real*4 a(n,n)
      a(i,j)=t
      return
      end

      subroutine array_r(a,n,i,j,t)
      real*4 a(n,n)
      t=a(i,j)
      return
      end

      subroutine imgrege(img1,msk1,nx1,ny1,nz1,mmppix1,center1,
     &                   img2,msk2,nx2,ny2,nz2,mmppix2,center2,param,mode,del,eta,q)
      real*4 img1(nx1,ny1,nz1)
      integer*2 msk1(nx1,ny1,nz1)
      real*4 mmppix1(3),center1(3)
      real*4 img2(nx2,ny2,nz2)
      integer*2 msk2(nx2,ny2,nz2)
      real*4 mmppix2(3),center2(3)
      real*4 param(12)
      real*4 t1(4,4),t2(4,4),grad1(3),grad2(3)
      real*4 x(3),b(4,4),c(4,4)
      real*4 eps/1./		! +/- interval (mm) over which to evaluate grad
      parameter (kmax=2097152)	! index in mask memory
      integer*2 kmem(kmax)
      integer*4 nk/0/
      logical*4 lcorrel,lsetk,lusek,lfind
      logical*4 ldebug/.false./
      external matmul

      lcorrel=iand(mode,256).ne.0
      lfind=  iand(mode,z'1000').ne.0
      lsetk=  iand(mode,z'20000').ne.0
      lusek=  iand(mode,z'10000').ne.0.and.(nk.gt.0).and..not.lsetk

      call vrt2img(mmppix1,center1,t1)
      call vrt2img(mmppix2,center2,b)
      call tparam2warp(mode,param,c)
      call matmul(b,c,t2,4)

      kz2=ifix(float(nz2)*abs(mmppix2(3))/del)
      kz1=ifix(-abs(center2(3))/del)
      kz2=kz2+kz1
      ky2=ifix(float(ny2)*abs(mmppix2(2))/del)
      ky1=ifix(-abs(center2(2))/del)
      ky2=ky2+ky1
      kx2=ifix(float(nx2)*abs(mmppix2(1))/del)
      kx1=ifix(-abs(center2(1))/del)
      kx2=kx2+kx1

      if(lsetk)then
        nk=(kx2-kx1+1)*(ky2-ky1+1)*(kz2-kz1+1)
        if(kmax.lt.nk)then
          write(*,"('imgrege:',i8,' exceeds index in mask memory capacity=',i8)")nk,kmax
          call exit(2)
        endif
        do k=1,kmax
          kmem(k)=0
        enddo
        nk=0
      endif

      n=0
      u10=0.
      u20=0.
      u11=0.
      u12=0.
      u22=0.
      k=1
      do 41 kz=kz1,kz2
      x(3)=float(kz)*del
      do 41 ky=ky1,ky2
      x(2)=float(ky)*del
      do 41 kx=kx1,kx2
      x(1)=float(kx)*del
      if(lusek.and.kmem(k).gt.0)goto 41
      if(lcorrel)then
        call imgvalm(img1,msk1,nx1,ny1,nz1,t1,x,v1,lslice)
        if(lslice.le.0)goto 49
        call imgvalm(img2,msk2,nx2,ny2,nz2,t2,x,v2,lslice)
        if(lslice.le.0)goto 49
        u10=u10+v1
        u20=u20+v2
        u12=u12+v1*v2
        u11=u11+v1*v1
        u22=u22+v2*v2
        n=n+1
      else
        call imgvalmg(img1,msk1,nx1,ny1,nz1,t1,x,eps,grad1,lslice)
        if(lslice.le.0)goto 49
        call imgvalmg(img2,msk2,nx2,ny2,nz2,t2,x,eps,grad2,lslice)
        if(lslice.le.0)goto 49
        a11=grad1(1)*grad1(1)+grad1(2)*grad1(2)+grad1(3)*grad1(3)
        a22=grad2(1)*grad2(1)+grad2(2)*grad2(2)+grad2(3)*grad2(3)
        a12=grad1(1)*grad2(1)+grad1(2)*grad2(2)+grad1(3)*grad2(3)
        q=a11*a22
        if(q.gt.0.)then
c         u12=u12+abs(a12)
          u12=u12+a12**2/sqrt(q)
          u11=u11+a11
          u22=u22+a22
          n=n+1
        endif
      endif
      goto 41
   49 if(lsetk)then
        kmem(k)=1
        nk=nk+1
      endif
   41 k=k+1
      if(ldebug)write(*, "('imrege: voxels used ',i6)")n
      if(n.gt.0)then
        eta=u12/sqrt(u11*u22)
      else
        eta=0.
        if(.not.lfind)then
          write(*,"('imgrege: no voxels in register')")
          call exit(-2)
        endif
      endif
      if(lfind)eta=eta*float(n)
      if(lcorrel)then
        q=u10/u20
      else
        q=0.
      endif
      return
      end
