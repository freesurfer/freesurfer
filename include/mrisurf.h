#ifndef MRISURF_H
#define MRISURF_H


#include "volume_io.h"

typedef struct face_type_
{
  int v[4];                      /* vertex numbers */
  int ripflag;                   /* ripped face */
#if 0
  float logshear,shearx,sheary;  /* compute_shear */
#endif
} face_type, FACE ;

typedef struct vertex_type_
{
  float x,y,z;           /* curr position */
  float nx,ny,nz;        /* curr normal */
#if 0
  float ox,oy,oz;        /* last position */
#endif
  float curv;            /* curr curvature */
#if 0
  float mx,my,mz;        /* last movement */
  float dipx,dipy,dipz;  /* dipole position */
  float dipnx,dipny,dipnz; /* dipole orientation */
  float nc;              /* curr length normal comp */
  float onc;             /* last length normal comp */
  float val;             /* scalar data value (file: rh.val, sig2-rh.w) */
  float val2;            /* complex comp data value (file: sig3-rh.w) */
  float valbak;          /* scalar data stack */
  float val2bak;         /* complex comp data stack */
  float oval;            /* last scalar data (for smooth_val) */
  float stat;            /* statistic */
  int undefval;          /* [previously dist=0] */
  int old_undefval;      /* for smooth_val_sparse */
  int fixedval;          /* [previously val=0] */
  float dist;            /* dist from sampled point [defunct: or 1-cos(a)] */
  float fieldsign;       /* fieldsign--final: -1,0,1 (file: rh.fs) */
  float fsmask;          /* significance mask (file: rh.fm) */
#endif
  int num;               /* number neighboring faces */
  int *f;                /* array neighboring face numbers */
  int *n;                /* [0-3, num long] */
  int vnum;              /* number neighboring vertices */
  int *v;                /* array neighboring vertex numbers, vnum long */
#if 0
  int   annotation;      /* area label (defunct--now from label file name!) */
  float stress;          /* explosion */
  float logarat,ologarat,sqrtarat; /* for area term */
  float logshear,shearx,sheary,oshearx,osheary;  /* for shear term */
  float bnx,bny,obnx,obny;                       /* boundary normal */
  float smx,smy,smz,osmx,osmy,osmz;              /* smoothed curr,last move */
  int   marked;            /* cursor */
  int   oripflag,origripflag;  /* cuts flags */
  float coord[3];
#endif
  int   ripflag ;
  int   border;            /* flag */
  float area,origarea ;
  int   tethered ;
} vertex_type, VERTEX ;

#if 0
typedef struct face2_type_
{
  int v[4];                      /* vertex numbers */
  int ripflag;                   /* ripped face */
} face2_type;

typedef struct vertex2_type_
{
  float x,y,z;           /* curr position */
  float nx,ny,nz;        /* curr normal */
  float curv;            /* curr curvature */
  int num;               /* number neighboring faces */
  int *f;                /* array neighboring face numbers */
  int *n;                /* [0-3, num long] */
  int vnum;              /* number neighboring vertices */
  int *v;                /* array neighboring vertex numbers, vnum long */
  int ripflag;           /* cut flags */
  int border;            /* flag */
} vertex2_type;
#endif

typedef struct
{
  int          nvertices ;      /* # of vertices on surface */
  int          nfaces ;         /* # of faces on surface */
  VERTEX       *vertices ;
  FACE         *faces ;
  float        xctr ;
  float        yctr ;
  float        zctr ;
  float        xlo ;
  float        ylo ;
  float        zlo ;
  float        xhi ;
  float        yhi ;
  float        zhi ;
  VERTEX       *v_temporal_pole ;
  VERTEX       *v_frontal_pole ;
  VERTEX       *v_occipital_pole ;
  float        max_curv ;
  float        min_curv ;
  int          total_area ;
  int          hemisphere ;            /* which hemisphere */

  General_transform transform ;   /* the next two are from this struct */
  Transform         *linear_transform ;
  Transform         *inverse_linear_transform ;
  int               free_transform ;
} MRI_SURFACE, MRIS ;




MRI_SURFACE  *MRISread(char *fname) ;
int          MRISwrite(MRI_SURFACE *mris) ;
MRI_SURFACE  *MRISalloc(int nvertices, int nfaces) ;
int          MRISfree(MRI_SURFACE **pmris) ;
MRI_SURFACE  *MRISprojectOntoEllipsoid(MRI_SURFACE *mris_src, 
                                       MRI_SURFACE *mris_dst, 
                                       float a, float b, float c) ;
MRI_SURFACE  *MRISclone(MRI_SURFACE *mris_src) ;
MRI_SURFACE  *MRIScenter(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
MRI_SURFACE  *MRIStalairachTransform(MRI_SURFACE *mris_src, 
                                    MRI_SURFACE *mris_dst);
MRI_SURFACE  *MRISunfold(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst,
                        int niterations, float l_area, float l_angle,
                        float l_corr);



/* constants for vertex->tethered */
#define TETHERED_NONE           0
#define TETHERED_FRONTAL_POLE   1
#define TETHERED_OCCIPITAL_POLE 2
#define TETHERED_TEMPORAL_POLE  3

#define MAX_TALAIRACH_Y         100.0f

/* constants for mris->hemisphere */
#define LEFT_HEMISPHERE         0
#define RIGHT_HEMISPHERE        1
#endif
