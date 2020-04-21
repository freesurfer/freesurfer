cRevision 1.1  2007/05/04 22:34:03  nicks
cnew talairach alignment utility, using Avi Snyders registration tools
c
c Revision 1.2  2007/04/29  05:14:20  avi
c gcc v3 compliant
c
c Revision 1.1  1997/10/11  07:02:48  avi
c Initial revision
c
      subroutine to_711_2b(t4)
c     convert atlas transform from target 711-2A to 711-2B
      real*4 t4(4,4)
      real*4 bt(4,4)/
     &  1.051053, -0.002200,  0.018579,   -0.4981,
     &  0.000148,  1.040993,  0.105308,    5.5848,
     & -0.019619, -0.108215,  1.004962,    0.9322,
     &  0.000000,  0.000000,  0.000000,    1.0000/
      real*4 b(4,4),a(4,4)
      external transpos,matmul,matcop

      call transpos(bt,b,4)
      call matmul(t4,b,a,4)
      call matcop(a,t4,4)
      return
      end
