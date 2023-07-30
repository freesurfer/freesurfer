c$Header: /space/repo/1/dev/dev/talairach_avi/t4_sub.f,v 1.2 2009/08/21 19:14:07 nicks Exp $
c$Log: t4_sub.f,v $
cRevision 1.2  2009/08/21 19:14:07  nicks
cadded width specifiers to allow compilation with gfortran
c
cRevision 1.1  2007/05/04 22:34:03  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.8  2007/04/25  05:04:07  avi
c gcc v3 compliant
c
c Revision 1.7  1999/10/08  01:12:15  avi
c iget_t4_scale()
c
c Revision 1.6  1999/01/05  04:27:40  avi
c t4_list()
c
c Revision 1.5  1998/12/26  04:50:17  avi
c Solaris
c do not echo t4 matrix on t4_read
c
c Revision 1.4  1998/02/05  13:42:00  avi
c fix disordered voxsiz adjust relative to rotation computation
c
c Revision 1.3  1998/02/05  10:38:56  avi
c correct setting of sag and cor bits when on together -
c still have to correct order of stretch operation
c
c Revision 1.2  1997/08/28  10:50:40  avi
c tparam2warp and warp2tparam modified to treat sagittal (-) images specially
c in rigid body transform mode keeping track of sagittal state with mode bit '80000'x
c
c Revision 1.1  1997/01/15  09:26:57  avi
c Initial revision

      subroutine t4_sub_rcsid
      write (*,"('t4_sub.f')")
      return
      end

      subroutine t4file2param(mode,t4_file,param)
      character*256 t4_file
      real*4 param(13)
      real*4 t4(4,4)
      character*256 string
      logical*4 lcorrel

      lcorrel=iand(mode,256).ne.0

      call t4_read(t4_file,t4)
      call warp2tparam(mode,t4,param)
      param(13)=1.

      l=index(t4_file,char(0))-1
      if(l.le.0)l=lnblnk(t4_file)
      open(8,file=t4_file(1:l),status='old',err=9)
    1 read(8,"(a)",end=9)string
      if(string(1:6).ne.'scale:')goto 1
      k=lnblnk(string)
      read(string(7:k),"(f10.6)",err=9)param(13)
    9 close(8)
      end

      integer*4 function iget_t4_scale(t4_file,scale)
      character*256 t4_file
      character*80 string

      scale=1.
      iget_t4_scale=1
      l=index(t4_file,char(0))-1
      if(l.le.0)l=lnblnk(t4_file)
      open(8,file=t4_file(1:l),status='old',err=9)
    1 read(8,"(a)",end=9)string
      if(string(1:6).ne.'scale:')goto 1
      k=lnblnk(string)
      read(string(7:k),"(f10.6)",err=9)scale
      iget_t4_scale=0
    9 close(8)
      end

      subroutine t4_list(t4)
      real*4 t4(4,4)

      write(*,"('t4')")
      do i=1,4
        write(*,"(3f10.6,f10.4)")(t4(i,j),j=1,4)
      enddo
      return
      end

      subroutine t4_read(t4_file,t4)
      character*256 t4_file
      real*4 t4(4,4)
      character*256 string
      logical*4 ldebug/.false./

      l=index(t4_file,char(0))-1
      if(l.le.0)l=lnblnk(t4_file)
      open(8,file=t4_file(1:l),status='old',err=10)
      do i=1,4
        read(8,"(3f10.6,f10.4)",err=7)(t4(i,j),j=1,4)
      enddo
      close (8)
      if(ldebug)call t4_list(t4)
      return
    7 rewind 8
    5 read(8,"(a)",end=9)string
      if(index(string,'t4').ne.1)goto 5
      do i=1,4
        read(8,"(3f10.6,f10.4)",err=9)(t4(i,j),j=1,4)
      enddo
      close(8)
      return
    9 write(*,"('t4_read: read error in t4_sub.f l:',i4,' t4_file(1:l)',a,' ')")l,t4_file(1:l)
      call t4_init(t4)
      write(*,"('t4_read: transform initialized to I4')")
      close(8)
      return
   10 write(*,"('t4_read: open error in t4_sub.f l:',i4,' t4_file(1:l)',a,' ')")l,t4_file(1:l)
      call t4_init(t4)
      write(*,"('t4_read: transform initialized to I4')")
      close(8)
      return
      end

      subroutine b6_read(b6_file,b6)
      character*256 b6_file
      real*4 b6(6,3)

      l=index(b6_file,char(0))-1
      if(l.le.0)l=lnblnk(b6_file)
      open(8,file=b6_file(1:l),status='old',err=99)
      do i=1,3
        read(8,"(6f10.6)")(b6(j,i),j=1,6)
      enddo
      close(8)
      write(*,"('b6')")
      do i=1,3
        write(*,"(6f10.6)")(b6(j,i),j=1,6)
      enddo
      return
   99 write(*,"('b6_read: ',a,' not found')")b6_file(1:l)
      do 1 i=1,3
      do 1 j=1,6
    1 b6(j,i)=0.
      write(*,"('t6_read: transform parameters cleared')")
      return
      end

      subroutine param2t4file0(mode,param,t4_file,string)
      real*4 param(13)
      character*256 t4_file,string
      real*4 t4(4,4)
      logical*4 lcorrel

      call tparam2warp(mode,param,t4)

      l=index(t4_file,char(0))-1
      if(l.le.0)l=lnblnk(t4_file)
      open(9,file=t4_file(1:l))
      k=index(string,char(0))-1
      if(k.eq.0)k=lnblnk(string)
      write(9,"(a)")string(1:k)
      do i=1,4
        write(9,"(3f10.6,f10.4)")(t4(i,j),j=1,4)
      enddo
      lcorrel=iand(mode,256).ne.0
      if(lcorrel)then
        write(9,"('scale:    ',f10.6)")param(13)
      endif
      close(9)

      return
      end

      subroutine t4_write(t4,t4_file)
      real*4 t4(4,4)
      character*256 t4_file

      l=index(t4_file,char(0))-1
      if(l.le.0)l=lnblnk(t4_file)
      open(9,file=t4_file(1:l))
      do i=1,4
        write(9,"(3f10.6,f10.4)")(t4(i,j),j=1,4)
      enddo
      close(9)
      return
      end

      subroutine t4_init(t4)
      real*4 t4(4,4)
      do i=1,4
        do j=1,4
          t4(i,j)=0.
        enddo
        t4(i,i)=1.
      enddo
      return
      end

      subroutine tparam2warp(mode,param,a)
      real*4 param(12),a(4,4)
      logical*4 lenable,l3d,lstretch,lvoxsiz,lreflect,lsag

      lenable= iand(mode,1).ne.0
      l3d=     iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0
      lvoxsiz= iand(mode,8).ne.0
      lreflect=.not.lstretch.and.iand(mode,z'40000').ne.0
      lsag=    .not.lstretch.and.iand(mode,z'80000').ne.0

      do l=1,4					! initialize a(4,4) to I
        do m=1,4
          a(l,m)=0.
        enddo
        a(l,l)=1.
      enddo
      if(.not.lenable)return

      if(lstretch)then				! stretch mode
        if(l3d)then
          lmax=3
        else
          lmax=2
        endif
        jj=1
        do l=1,lmax
          a(l,4)=param(jj)
          jj=jj+1
        enddo
        do l=1,lmax
          a(l,l)=param(jj)+1.
          jj=jj+1
        enddo
        do l=1,lmax
          do m=1,lmax
            if(m.ne.l)then
              a(l,m)=param(jj)
              jj=jj+1
            endif
          enddo
        enddo
      else					! .not.lstretch
        if(l3d)then				! 3d
          call trotset(param(1),param(4),a)
          if(lsag)then
            do j=1,3
              t     = a(1,j)
              a(1,j)=-a(2,j)
              a(2,j)= a(3,j)
              a(3,j)=-t
            enddo
          endif
          if(lreflect)then
            do j=1,3
              t     =a(2,j)
              a(2,j)=a(3,j)
              a(3,j)=t
            enddo
          endif
          if(lvoxsiz)then			! 3d voxsiz processing
            mskbit=16				! object image voxel size adjust dis-enable bits
            k=0
            do i=1,3
              if(iand(mode,mskbit).eq.0)then
                k=k+1
                do j=1,3
                  a(i,j)=a(i,j)*(1.+param(6+k))
                enddo
              endif
              mskbit=mskbit*2
            enddo
          endif
        else					! 2d
          a(1,4)=param(1)
          a(2,4)=param(2)
          a(1,1)= cos(param(3))
          a(1,2)=-sin(param(3))
          a(2,1)= sin(param(3))
          a(2,2)= cos(param(3))
          if(lvoxsiz)then			! 2d voxsiz processing
            mskbit=16				! object image voxel size adjust dis-enable bits
            k=0
            do i=1,2
              if(iand(mode,mskbit).eq.0)then
                k=k+1
                do j=1,3
                  a(i,j)=a(i,j)*(1.+param(3+k))
                enddo
              endif
              mskbit=mskbit*2
            enddo
          endif
        endif
      endif					! end .not.lstretch 

      return
      end

      subroutine warp2tparam(mode,t4,param)
      real*4 t4(4,4),param(12)
      real*4 rot(3,3),g(3,3),h(3,3),a(4,4)
      logical*4 lenable,l3d,lstretch,lvoxsiz,lreflect,lsag
      logical*4 ldebug/.true./
      external matmul

      lenable= iand(mode,1).ne.0
      l3d=     iand(mode,2).ne.0
      lstretch=iand(mode,4).ne.0
      lvoxsiz= iand(mode,8).ne.0

      if(.not.lenable)return

      do 1 i=1,4
      do 1 j=1,4
    1 a(i,j)=t4(i,j)

      if(l3d)then
        lmax=3
      else
        lmax=2
      endif
      jj=1
      do l=1,lmax
        param(jj)=a(l,4)
        jj=jj+1
      enddo

      if(.not.lstretch)then
        do i=1,3
          do j=1,3
            g(i,j)=0.
            h(i,j)=a(i,j)
          enddo
        enddo
        do i=1,3
          g(i,i)=1./sqrt(a(i,1)**2+a(i,2)**2+a(i,3)**2)
        enddo
        call matmul(g,h,rot,3)
        if(lvoxsiz)then				! start voxsiz processing
          mskbit=16
          k=0
          if(l3d)then				! process 3d voxsiz
            do i=1,3
              if(iand(mode,mskbit).eq.0)then
                k=k+1
                param(6+k)=1./g(i,i)-1.
              endif
              mskbit=mskbit*2
            enddo
          else					! process 2d voxsiz
            do i=1,2
              if(iand(mode,mskbit).eq.0)then
                k=k+1
                param(3+k)=1./g(i,i)-1.
              endif
              mskbit=mskbit*2
            enddo
          endif
        endif					! end voxsiz processing
        det=(rot(1,1)*(rot(2,2)*rot(3,3)-rot(3,2)*rot(2,3))
     &      -rot(1,2)*(rot(2,1)*rot(3,3)-rot(3,1)*rot(2,3))
     &      +rot(1,3)*(rot(2,1)*rot(3,2)-rot(3,1)*rot(2,2)))
        if(ldebug)write(*, "('rotation matrix determinant',f10.6)")det
        lreflect=l3d.and.(det.lt.0.)
        if(lreflect)then
          do j=1,3
            t       =rot(2,j)
            rot(2,j)=rot(3,j)
            rot(3,j)=t
          enddo
          mode=ior(mode,z'40000')		! set negative determinant bit
        endif
        lsag=l3d.and.(abs(rot(3,1)).gt.sqrt(0.5))
        if(lsag)then
          do j=1,3
            t       = rot(1,j)
            rot(1,j)=-rot(3,j)
            rot(3,j)= rot(2,j)
            rot(2,j)=-t
          enddo
          mode=ior(mode,z'80000')		! set sagittal bit
        endif
        if(l3d)then				! 3d rot->ang
          call rot2ang(rot,param(4))
        else					! 2d rot->ang
          param(1)=a(1,4)
          param(2)=a(2,4)
          param(3)=atan2(rot(2,1),rot(2,2))
        endif
      else					! stretch mode
        do l=1,lmax
          param(jj)=a(l,l)-1.
          jj=jj+1
        enddo
        do l=1,lmax
          do m=1,lmax
            if(m.ne.l)then
              param(jj)=a(l,m)
              jj=jj+1
            endif
          enddo
        enddo
      endif

      return
      end

      subroutine trns30(mode,t4,b6,pnti,pnta)
      real*4 t4(4,4),b6(6,3),pnti(3),pnta(3)
      real*4 tcenter(3)/0.,-23.5,5./		! Talairach space approximate geometric center
      real*4 pntn(3),pntc(3)			! domain normalized input coordinates
      logical*4 linear

      linear=iand(mode,128).eq.0
      do k=1,3
        pnta(k)=t4(k,1)*pnti(1)+t4(k,2)*pnti(2)+t4(k,3)*pnti(3)+t4(k,4)
      enddo
      if(linear)return

      pntc(1)=pnti(1)-tcenter(1)
      pntc(2)=pnti(2)-tcenter(2)
      pntc(3)=pnti(3)-tcenter(3)
      pntn(1)=pntc(1)/60.
      pntn(2)=pntc(2)/85.
      pntn(3)=pntc(3)/60.
      do k=1,3
        pnta(k)=pnta(k)+b6(1,k)*pntn(1)*pntc(1)
     &                 +b6(2,k)*pntn(2)*pntc(2)
     &                 +b6(3,k)*pntn(3)*pntc(3)
     &                 +b6(4,k)*pntn(1)*pntc(2)
     &                 +b6(5,k)*pntn(1)*pntc(3)
     &                 +b6(6,k)*pntn(2)*pntc(3)
      enddo

      return
      end
