/**
 * @file  mrisurf.h
 * @brief MRI_SURFACE utilities.
 *
 * Utilities, constants and structure definitions for manipulation
 * and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/01/08 15:41:50 $
 *    $Revision: 1.374 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef MRISURF_H
#define MRISURF_H

#include "minc_volume_io.h"
#include "const.h"

#define MAX_SURFACES 20
#define TALAIRACH_COORDS     0
#define SPHERICAL_COORDS     1
#define ELLIPSOID_COORDS     2

#define VERTICES_PER_FACE    3
#define ANGLES_PER_TRIANGLE  3

#define INFLATED_NAME        "inflated"
#define SMOOTH_NAME          "smoothwm"
#define SPHERE_NAME          "sphere"
#define ORIG_NAME            "orig"
#define WHITE_MATTER_NAME    "white"
#define GRAY_MATTER_NAME     "gray"
#define LAYERIV_NAME         "graymid"
#define GRAYMID_NAME         LAYERIV_NAME
#define MAX_CMDS 1000

#define NEW_VERSION_MAGIC_NUMBER  16777215 // was in mrisurf.c

typedef struct _area_label
{
  char     name[STRLEN] ;     /* name of region */
  float    cx ;               /* centroid x */
  float    cy ;               /* centroid y */
  float    cz ;               /* centroid z */
  int      label ;            /* an identifier (used as an index) */
}
MRIS_AREA_LABEL ;

/*
  the vertices in the face structure are arranged in
  counter-clockwise fashion when viewed from the outside.
*/
typedef struct face_type_
{
  int    v[VERTICES_PER_FACE];           /* vertex numbers of this face */
  float  nx ;
  float  ny ;
  float  nz ;
  float  area ;
  float  orig_area ;
  float  angle[ANGLES_PER_TRIANGLE] ;
  float  orig_angle[ANGLES_PER_TRIANGLE]  ;
  char   ripflag;                        /* ripped face */
  char   oripflag;                       /* stored version */
  int    marked;                         /* marked face */
#if 0
  float logshear,shearx,sheary;  /* compute_shear */
#endif
  float  cx, cy, cz ;         // coordinates of centroid
}
face_type, FACE ;

#ifndef uchar
#define uchar  unsigned char
#endif

#include "colortab.h" // 'COLOR_TABLE'

typedef struct vertex_type_
{
  float x,y,z;           /* curr position */
  float nx,ny,nz;        /* curr normal */
  float pnx,pny,pnz;     /* pial normal */
  float wnx,wny,wnz;     /* white normal */
  float onx,ony,onz;     /* original normal */
  float dx, dy, dz ;     /* current change in position */
  float odx, ody, odz ;  /* last change of position (for momentum) */
  float tdx, tdy, tdz ;  /* temporary storage for averaging gradient */
  float curv;            /* curr curvature */
  float curvbak ;
  float val;             /* scalar data value (file: rh.val, sig2-rh.w) */
  float imag_val ;       /* imaginary part of complex data value */
  float cx, cy, cz ;     /* coordinates in canonical coordinate system */
  float tx, ty, tz ;     /* tmp coordinate storage */
  float tx2, ty2, tz2 ;  /* tmp coordinate storage */
  float origx, origy, origz ;   /* original coordinates */
  float targx, targy, targz ;   // target coordinates
  float pialx, pialy, pialz ;   /* pial surface coordinates */
  float whitex, whitey, whitez ;/* white surface coordinates */
  float l4x, l4y, l4z ;   /* layerIV surface coordinates */
  float infx, infy, infz; /* inflated coordinates */
  float fx, fy, fz ;      /* flattened coordinates */
  int   px,qx, py,qy, pz,qz; /* rational coordinates for exact calculations */
  float e1x, e1y, e1z ;  /* 1st basis vector for the local tangent plane */
  float e2x, e2y, e2z ;  /* 2nd basis vector for the local tangent plane */
  float pe1x, pe1y, pe1z ;  /* 1st basis vector for the local tangent plane */
  float pe2x, pe2y, pe2z ;  /* 2nd basis vector for the local tangent plane */
#if 0
  float dipx,dipy,dipz;  /* dipole position */
  float dipnx,dipny,dipnz; /* dipole orientation */
#endif
  float nc;              /* curr length normal comp */
  float val2;            /* complex comp data value (file: sig3-rh.w) */
  float valbak;          /* scalar data stack */
  float val2bak;         /* complex comp data stack */
  float stat;            /* statistic */
#if 1
  int undefval;          /* [previously dist=0] */
  int old_undefval;      /* for smooth_val_sparse */
  int fixedval;          /* [previously val=0] */
#endif
  float fieldsign;       /* fieldsign--final: -1,0,1 (file: rh.fs) */
  float fsmask;          /* significance mask (file: rh.fm) */
  uchar num;             /* number neighboring faces */
  int   *f;              /* array neighboring face numbers */
  uchar *n;              /* [0-3, num long] */
  uchar vnum;            /* number neighboring vertices */
  int   *v;              /* array neighboring vertex numbers, vnum long */
  int   *e;              /* edge state for neighboring vertices */
  int    v2num ;         /* number of 2-connected neighbors */
  int    v3num ;         /* number of 3-connected neighbors */
  short  vtotal ;        /* total # of neighbors,
                                    will be same as one of above*/
  float d ;              /* for distance calculations */
  uchar nsize ;          /* size of neighborhood (e.g. 1, 2, 3) */
#if 0
  float *tri_area ;      /* array of triangle areas - num long */
  float *orig_tri_area ; /* array of original triangle areas - num long */
  float dist;            /* dist from sampled point [defunct: or 1-cos(a)] */
  float ox,oy,oz;        /* last position (for undoing time steps) */
  float mx,my,mz;        /* last movement */
  float onc;             /* last length normal comp */
  float oval;            /* last scalar data (for smooth_val) */
  float *fnx ;           /* face normal - x component */
  float *fny ;           /* face normal - y component */
  float *fnz ;           /* face normal - z component */
  float bnx,bny,obnx,obny; /* boundary normal */
  float *tri_angle ;     /* angles of each triangle this vertex belongs to */
  float *orig_tri_angle ;/* original values of above */
  float stress;          /* explosion */
  float logshear,shearx,sheary,oshearx,osheary;  /* for shear term */
  float ftmp ;           /* temporary floating pt. storage */
  float logarat,ologarat,sqrtarat; /* for area term */
  float smx,smy,smz,osmx,osmy,osmz;/* smoothed curr,last move */
#endif
  int   annotation;      /* area label (defunct--now from label file name!) */
  char   oripflag,origripflag;  /* cuts flags */
#if 0
  float coords[3];
#endif
  float theta, phi ;     /* parameterization */
  short  marked;         /* for a variety of uses */
  short  marked2 ;
  short  marked3 ;
  char   ripflag ;
  char   border;         /* flag */
  float area, origarea, group_avg_area ;
  float K ;              /* Gaussian curvature */
  float H ;              /* mean curvature */
  float k1, k2 ;         /* the principal curvatures */
  float *dist ;          /* original distance to neighboring vertices */
  float *dist_orig ;     /* original distance to neighboring vertices */
  char   neg ;           /* 1 if the normal vector is inverted */
  float mean ;
  float mean_imag ;      /* imaginary part of complex statistic */
  float std_error ;
  unsigned int flags ;
  void *vp; /* to store user's information */
  int   linked ;         // is this vertex linked to some others?
  int   fno ;            // face that this vertex is in
}
vertex_type, VERTEX ;

#define VERTEX_SULCAL  0x00000001L

typedef struct
{
  int nvertices;
  unsigned int *vertex_indices;
}
STRIP;

#include "transform.h" // TRANSFORM, LTA

typedef struct
{
  int          nvertices ;      /* # of vertices on surface */
  int          nfaces ;         /* # of faces on surface */
  int          nstrips;
  VERTEX       *vertices ;
  FACE         *faces ;
  STRIP        *strips;
  float        xctr ;
  float        yctr ;
  float        zctr ;
  float        xlo ;
  float        ylo ;
  float        zlo ;
  float        xhi ;
  float        yhi ;
  float        zhi ;
  float        x0 ;   // center of spherical expansion
  float        y0 ;
  float        z0 ;
  VERTEX       *v_temporal_pole ;
  VERTEX       *v_frontal_pole ;
  VERTEX       *v_occipital_pole ;
  float        max_curv ;
  float        min_curv ;
  float        total_area ;
  double       avg_vertex_area;
  double       avg_vertex_dist;
  double       std_vertex_dist;
  float        orig_area ;
  float        neg_area ;
  float        neg_orig_area ;   /* amount of original surface in folds */
  int          zeros ;
  int          hemisphere ;      /* which hemisphere */
  int          initialized ;
#if 0
  General_transform transform ;   /* the next two are from
                                             this struct (MNI transform) */
  Transform         *linear_transform ;
  Transform         *inverse_linear_transform ;
#endif
  LTA         *lta;
  MATRIX      *SRASToTalSRAS_;
  MATRIX      *TalSRASToSRAS_;
  int          free_transform ;
  double       radius ;           /* radius (if status==MRIS_SPHERE) */
  float        a, b, c ;          /* ellipsoid parameters */
  char         fname[STRLEN] ;    /* file it was originally loaded from */
  float        Hmin ;             /* min mean curvature */
  float        Hmax ;             /* max mean curvature */
  float        Kmin ;             /* min Gaussian curvature */
  float        Kmax ;             /* max Gaussian curvature */
  double       Ktotal ;           /* total Gaussian curvature */
  int          status ;           /* type of surface (e.g. sphere, plane) */
  int          patch ;            /* if a patch of the surface */
  int          nlabels ;
  MRIS_AREA_LABEL *labels ;       /* nlabels of these (may be null) */
  int          nsize ;            /* size of neighborhoods */
  float        avg_nbrs ;         /* mean # of vertex neighbors */
  void         *vp ;              /* for misc. use */
  float        alpha ;            /* rotation around z-axis */
  float        beta ;             /* rotation around y-axis */
  float        gamma ;            /* rotation around x-axis */
  float        da, db, dg ;       /* old deltas */
  int          type ;             /* what type of surface was this initially*/
  int          max_vertices ;     /* may be bigger than nvertices */
  int          max_faces ;        /* may be bigger than nfaces */
  char         subject_name[STRLEN] ;/* name of the subject */
  float        canon_area ;
  int          noscale ;          /* don't scale by surface area if true */
  float        *dx2 ;             /* an extra set of gradient
                                     (not always alloced) */
  float        *dy2 ;
  float        *dz2 ;
  COLOR_TABLE  *ct ;
  int          useRealRAS;        /* if 0, vertex position is a
                                     conformed volume RAS with c_(r,a,s)=0 */
  /* if 1, verteix position is a
     real RAS (volume stored RAS)         */
  /* The default is 0.        */
  VOL_GEOM     vg;                /* volume information from which
                                     this surface is created.
                                     check validity by vg.valid = 1 or not */
  char   *cmdlines[MAX_CMDS] ;
  int    ncmds;
  float  group_avg_surface_area ;  // average of total surface area for group
  int    group_avg_vtxarea_loaded; /* average vertex area for group
                                              at each vertex */
  int    triangle_links_removed ;  // for quad surfaces
  void   *user_parms ;             // for whatever the user wants to hang here 
  MATRIX *m_sras2vox ;             // for converting surface ras to voxel 
  MRI    *mri_sras2vox ;           // volume that the above matrix is for
  void   *mht ;
}
MRI_SURFACE, MRIS ;

#define IPFLAG_HVARIABLE                0x0001 /* for parms->flags */
#define IPFLAG_NO_SELF_INT_TEST         0x0002
#define IPFLAG_QUICK                    0x0004 /* sacrifice qualty for speed */
#define IPFLAG_ADD_VERTICES             0x0008
#define IP_USE_CURVATURE                0x0010         /* was 0x0001 !!! */
#define IP_NO_RIGID_ALIGN               0x0020         /* was 0x0002 !!! */
#define IP_RETRY_INTEGRATION            0x0040         /* was 0x0004 !!! */
/* VECTORIAL_REGISTRATION*/
#define IP_USE_MULTIFRAMES              0x0080
#define IP_NO_SULC                      0x0100
#define IP_USE_INFLATED                 0x0200
#define IP_NO_FOLD_REMOVAL              0x0400
#define IP_DONT_ALLOW_FOLDS             0x0800
/* MRIScorrectTopology : topology preserving patch deformation */
#define IPFLAG_PRESERVE_TOPOLOGY_CONVEXHULL 0x1000 /* apply topology
preserving gradient */
#define IPFLAG_PRESERVE_SPHERICAL_POSITIVE_AREA 0x2000 /* apply gradients
that preserve
positive areas */
#define IPFLAG_MAXIMIZE_SPHERICAL_POSITIVE_AREA 0x4000  /* apply  gradients
that will
maximize the
positive areas */
#define IPFLAG_NOSCALE_TOL            0x8000   // don't scale tol with navgs
#define IPFLAG_FORCE_GRADIENT_OUT    0x10000
#define IPFLAG_FORCE_GRADIENT_IN     0x20000

#define INTEGRATE_LINE_MINIMIZE    0  /* use quadratic fit */
#define INTEGRATE_MOMENTUM         1
#define INTEGRATE_ADAPTIVE         2
#define INTEGRATE_LM_SEARCH        3  /* binary search for minimum */

#define MRIS_SURFACE               0
#define MRIS_PATCH                 1
#define MRIS_CUT                   MRIS_PATCH
#define MRIS_PLANE                 2
#define MRIS_ELLIPSOID             3
#define MRIS_SPHERE                4
#define MRIS_PARAMETERIZED_SPHERE  5
#define MRIS_RIGID_BODY            6
#define MRIS_SPHERICAL_PATCH       7
#define MRIS_UNORIENTED_SPHERE     8
#define MRIS_PIAL_SURFACE          9

// different Hausdorff distance modes
#define HDIST_MODE_SYMMETRIC_MEAN 0
double MRIScomputeHausdorffDistance(MRI_SURFACE *mris, int mode) ;


/*
  the following structure is built after the surface has been morphed
  into a parameterizable surface (such as a sphere or an ellipsoid). It
  contains relevant data (such as curvature or sulc) which represented
  in the plane of the parameterization. This then allows operations such
  as convolution to be carried out directly in the parameterization space
*/

/*
   define it for a sphere, the first axis (x or u) goes from 0 to pi,
   the second axis (y or v) goes from 0 to 2pi.
*/
#define PHI_MAX                   M_PI
#define U_MAX                     PHI_MAX
#define X_DIM(mrisp)              (mrisp->Ip->cols)
#define U_DIM(mrisp)              X_DIM(mrisp)
#define PHI_DIM(mrisp)            X_DIM(mrisp)
#define PHI_MAX_INDEX(mrisp)      (X_DIM(mrisp)-1)
#define U_MAX_INDEX(mrisp)        (PHI_MAX_INDEX(mrisp))

#define THETA_MAX                 (2.0*M_PI)
#define Y_DIM(mrisp)              (mrisp->Ip->rows)
#define V_DIM(mrisp)              (Y_DIM(mrisp))
#define THETA_DIM(mrisp)          (Y_DIM(mrisp))
#define THETA_MAX_INDEX(mrisp)    (Y_DIM(mrisp)-1)
#define V_MAX_INDEX(mrisp)        (THETA_MAX_INDEX(mrisp))

#include "image.h" // IMAGE
typedef struct
{
  MRI_SURFACE  *mris ;        /* surface it came from (if any) */
  IMAGE        *Ip ;
#if 0
  /* 2-d array of curvature, or sulc in parms */
  float        data[X_DIM][Y_DIM] ;
  float        distances[X_DIM][Y_DIM] ;
  int          vertices[X_DIM][Y_DIM] ;   /* vertex numbers */
#endif
  float        sigma ;                    /* blurring scale */
  float        radius ;
}
MRI_SURFACE_PARAMETERIZATION, MRI_SP ;

#define L_ANGLE              0.25f /*was 0.01*/ /* coefficient of angle term */
#define L_AREA               1.0f    /* coefficient of angle term */
#define N_AVERAGES           4096
#define WRITE_ITERATIONS     10
#define NITERATIONS          1
#define NO_PROJECTION        0
#define PROJECT_ELLIPSOID    1
#define ELLIPSOID_PROJECTION PROJECT_ELLIPSOID
#define PROJECT_PLANE        2
#define PLANAR_PROJECTION    PROJECT_PLANE
#define PROJECT_SPHERE       3

#define MOMENTUM             0.8
#define EPSILON              0.25

#define TOL                  1e-6  /* minimum error tolerance for unfolding */
#define DELTA_T              0.1

/* VECTORIAL_REGISTRATION */
#include "field_code.h"

typedef struct
{
  double  tol ;               /* tolerance for terminating a step */
  float   l_angle ;           /* coefficient of angle term */
  float   l_area ;            /* coefficient of (negative) area term */
  float   l_parea ;           /* coefficient of (all) area term */
  float   l_nlarea ;          /* coefficient of nonlinear area term */
  float   l_nldist ;          /* coefficient of nonlinear distance term */
  float   l_thick_normal ;    // coefficient to keep vector field close to surface normal
  float   l_thick_spring ;    // coefficient to keep vector field spaced on pial surface
  float   l_ashburner_triangle ;   // coefficient for (Ashburner, 1999) invertibility term
  float   l_ashburner_lambda ;  // amount of regularization
  float   l_corr ;            /* coefficient of correlation term */
  float   l_ocorr ;           // overlay correlation weight
  float   l_pcorr ;           /* polar correlation for rigid body */
  float   l_curv ;            /* coefficient of curvature term */
  float   l_scurv ;           /* coefficient of curvature term */
  float   l_lap ;             // coefficient of laplacian term
  float   l_link ;            /* coefficient of link term to keep
                                 white and pial vertices approximately
                                 along normal */
  float   l_spring ;          /* coefficient of spring term */
  float   l_nlspring ;        /* nonlinear spring term */
  float   l_max_spring ;      /* draws vertices towards the most distant neighbor */
  float   l_spring_norm ;     /* coefficient of normalize spring term */
  float   l_tspring ;         /* coefficient of tangential spring term */
  float   l_nltspring ;       /* coefficient of nonlinear tangential spring term */
  float   l_nspring ;         /* coefficient of normal spring term */
  float   l_repulse ;         /* repulsize force on tessellation */
  float   l_repulse_ratio ;   /* repulsize force on tessellation */
  float   l_boundary ;        /* coefficient of boundary term */
  float   l_dist ;            /* coefficient of distance term */
  float   l_location ;        // target location term
  float   l_neg ;
  float   l_intensity ;       /* for settling surface at a specified val */
  float   l_sphere ;          /* for expanding the surface to a sphere */
  float   l_expand ;          /* for uniformly expanding the surface */
  float   l_grad ;            /* gradient term */
  float   l_convex ;          /* convexity term */
  float   l_tsmooth ;         /* thickness smoothness term */
  float   l_surf_repulse ;    /* repulsive orig surface (for white->pial) */
  float   l_osurf_repulse ;   /* repulsive outer surface (for layer IV) */
  float   l_external ;        /* external (user-defined) coefficient */
  float   l_thick_parallel ;  // term that encourages thickness vectors to be parallel
  float   l_thick_min ;       // term that encourages thickness vectors to be minimal
  float   l_shrinkwrap ;      /* move in if MRI=0 and out otherwise */
  float   l_expandwrap ;      /* move out */
  float   l_unfold ;          /* move inwards along normal */
  float   l_dura ;            // move away from dura
  float   l_histo ;           // increase the likelihood of the entire volume given the surfaces
  double  l_map ;             // for MAP deformation
  double  l_map2d ;             // for 2D MAP deformation (intensity x distance)
  double  dura_thresh ;
  MRI     *mri_dura ;         /* ratio of early to late echo -
                                         dura shows up bright */
  int     n_averages ;        /* # of averages */
  int     min_averages ;
  int     first_pass_averages;// # of averages to use in the first pass
  int     nbhd_size ;
  int     max_nbrs ;
  int     write_iterations ;  /* # of iterations between saving movies */
  char    base_name[STRLEN] ; /* base name of movie files */
  int     projection ;        /* what kind of projection to do */
  int     niterations ;       /* max # of time steps */
  float   a ;
  float   b ;
  float   c ;                 /* ellipsoid parameters */
  int     start_t ;           /* starting time step */
  int     t ;                 /* current time */
  FILE    *fp ;               /* for logging results */
  float   Hdesired ;          /* desired (mean) curvature */
  int     integration_type ;  /* line minimation or momentum */
  double  momentum ;
  double  dt ;                /* current time step (for momentum only) */
  double  base_dt ;           /* base time step (for momentum only) */
  int     flags ;
  double  dt_increase ;       /* rate at which time step increases */
  double  dt_decrease ;       /* rate at which time step decreases */
  double  error_ratio ;       /* ratio at which to undo step */
  double  epsilon ;           /* constant in Sethian inflation */
  double  desired_rms_height; /* desired height above tangent plane */
  double  starting_sse ;
  double  ending_sse ;
  double  scale ;             /* scale current distances to mimic spring */
  MRI_SP  *mrisp ;            /* parameterization  of this surface */
  int     frame_no ;          /* current frame in template parameterization */
  MRI_SP  *mrisp_template ;   /* parameterization of canonical surface */
  MRI_SP  *mrisp_blurred_template ; /* parameterization of canonical
                                       surface convolve with Gaussian */
  double  area_coef_scale ;
  float   sigma ;             /* blurring scale */

  /* VECTORIAL_REGISTRATION
     The average template mrisp is assumed to be composed of several
     different fields (look in 'field_code.h').
     MRISvectorRegistration will use 'nfields' fields,
     with their corresponding
     location at fields[n].frames, 0 <= n< nfields in the mrisp structure.
     The field code is in fields[n].field (look in 'field_code.h').
     Corresponding correlation terms are in
     fields[n].l_corrs and fields[n].l_pcorrs.
     MRISvectorRegistration will use the structure VALS_VP in v->vp
  */
#define MAX_NUMBER_OF_FIELDS_IN_VECTORIAL_REGISTRATION  50
#define MNOFIV MAX_NUMBER_OF_FIELDS_IN_VECTORIAL_REGISTRATION
  int nfields;                  /* the number of fields in mrisp */
  FIELD_LABEL fields[MNOFIV];   /* information for each field */
  /*    int     ncorrs;            /\* the number of fields in mrisp *\/ */
  /*    int     corrfields[MNOFIV];/\* field code (see below) *\/ */
  /*    int     frames[MNOFIV];    /\* corresponding frame in mrisp *\/  */
  /*    int     types[MNOFIV];     /\* the field type
    (default,distance field...) *\/ */
  /*    float   l_corrs[MNOFIV];   /\* correlation coefficient *\/ */
  /*    float   l_pcorrs[MNOFIV];  /\* polar correlation coefficient *\/ */
  /*    float   sses[MNOFIV];      /\* corresponding sse *\/ */

  MRI     *mri_brain ;        /* for settling surfaace to e.g. g/w border */
  MRI     *mri_smooth ;       /* smoothed version of mri_brain */
  void    *user_parms ;       /* arbitrary spot for user to put stuff */
  MRI     *mri_dist ;         /* distance map for repulsion term */
  float   target_radius ;
  int     ignore_energy ;     // when no valid energy func availb...integrate
  int     check_tol ;         // to avoid changing mris_make_surfaces
  char    *overlay_dir;       // subject/overlay_dir/parms->fields[n].fname
  int     nsurfaces ;         // if 0 use default
  MRI     *mri_ll ;           // log-likelihood image
  double  rmin ;
  double  rmax ;              // for nonlinear spring term
  int     var_smoothness ;    // for space/time varying weights on
                              // metric distortion/likelihood
  float   *vsmoothness ;      // variable smoothness coefficients
                              // (one per vertex)
  float   *dist_error ;       // the values for each vertex of the
                              // various error type
  float   *area_error ;       // needed for the variable smoothness
                              // gradient calculation
  float   *geometry_error ;
  int     which_norm ;        // mean or median normalization
  int     abs_norm ;
  int     grad_dir ;          // use this instead of gradient direction
  int     fill_interior ;     // use filled interior to constrain gradient to not leave surface
  double  rms ;
  int     complete_dist_mat ; //whether to sample or use complete dist mat
  int     nsubjects ;
  int     nlabels ;
  void       **mht_array;
  MRI_SURFACE **mris_array ;
  MRI_SURFACE *mris_ico ;     // for sampling from a spherical template
  void         *mht ;        // hash table for surface vertex/face lookup
  int          smooth_averages ;
  int          ico_order ;  // which icosahedron to use
  int          remove_neg ;
  MRI          *mri_hires ;
  MRI          *mri_hires_smooth ;
  MRI          *mri_vno ;
  MRI          *mri_template ;
  int          which_surface ;
  float        trinarize_thresh;   // for trinarizing curvature in registration
  int          smooth_intersections ;  // run soap bubble smoothing during surface positioning
  int          uncompress ;            // run code to remove compressions in tessellation
  double       min_dist ;
  HISTOGRAM    *h_wm ;
  HISTOGRAM    *h_gm ;
  HISTOGRAM    *h_nonbrain ;
  MRI          *mri_labels ;   // hires labeling of interior of WM, GM and nonbrain
  MRI          *mri_white ;
  MRI          *mri_aseg ;
  HISTOGRAM    *hwm ;
  HISTOGRAM    *hgm ;
  HISTOGRAM    *hout ;

  HISTOGRAM2D  *h2d_wm ;
  HISTOGRAM2D  *h2d_gm ;
  HISTOGRAM2D  *h2d_out ;
  HISTOGRAM2D  *h2d ;
  MRI          *mri_volume_fractions ;  // the partial volume fractions associated with the boundaries in this mris
  MRI          *mri_dtrans ;   // distance to surface
  float        resolution ;  // at which to compute distance transforms and such
}
INTEGRATION_PARMS ;

extern double (*gMRISexternalGradient)(MRI_SURFACE *mris,
                                         INTEGRATION_PARMS *parms) ;
extern double (*gMRISexternalSSE)(MRI_SURFACE *mris,
                                    INTEGRATION_PARMS *parms) ;
extern double (*gMRISexternalRMS)(MRI_SURFACE *mris,
                                    INTEGRATION_PARMS *parms) ;
extern int (*gMRISexternalTimestep)(MRI_SURFACE *mris,
                                      INTEGRATION_PARMS *parms) ;
extern int (*gMRISexternalRipVertices)(MRI_SURFACE *mris,
                                         INTEGRATION_PARMS *parms);
extern int (*gMRISexternalClearSSEStatus)(MRI_SURFACE *mris) ;
extern int (*gMRISexternalReduceSSEIncreasedGradients)(MRI_SURFACE *mris,
      double pct) ;

// These are for backwards compatibility for when we don't want to fix
// the vertex area. Default is to fix, but this can be changed
// by setting to 0.
int MRISsetFixVertexAreaValue(int value);
int MRISgetFixVertexAreaValue(void);

/* The following structure is used in MRISvectorRegistration
   The original values are loaded in orig_vals
   The used values are in vals (for instance the  blurred orig_vals)
*/
typedef struct
{
  float *vals;
  float *orig_vals;
  int nvals;
}
VALS_VP;

const char *MRISurfSrcVersion(void);

/* new functions */
int MRISvectorRegister(MRI_SURFACE *mris,
                       MRI_SP *mrisp_template,
                       INTEGRATION_PARMS *parms, int max_passes,
                       float min_degrees, float max_degrees,
                       int nangles);
int MRISrigidBodyAlignVectorLocal(MRI_SURFACE *mris,
                                  INTEGRATION_PARMS *old_parms);
int MRISrigidBodyAlignVectorGlobal(MRI_SURFACE *mris,
                                   INTEGRATION_PARMS *parms,
                                   float min_degrees,
                                   float max_degrees, int nangles);
void MRISnormalizeField(MRIS *mris , int distance_field, int which_norm);
int  MRISsmoothCurvatures(MRI_SURFACE *mris, int niterations) ;
void MRISsetCurvaturesToValues(MRIS *mris,int fno);
void MRISsetCurvaturesToOrigValues(MRIS *mris,int fno);
void MRISsetOrigValuesToCurvatures(MRIS *mris,int fno);
void MRISsetOrigValuesToValues(MRIS *mris,int fno);

MRI *MRISbinarizeVolume(MRI_SURFACE *mris,
                        MRI_REGION *region,
                        float resolution,
                        float distance_from_surface);

int MRISfindClosestCanonicalVertex(MRI_SURFACE *mris, float x, float y,
                                   float z) ;
int MRISfindClosestOriginalVertex(MRI_SURFACE *mris, float x, float y,
                                  float z) ;
double MRIScomputeWhiteVolume(MRI_SURFACE *mris, 
                              MRI *mri_aseg, 
                              double resolution) ;
MRI *MRISfillWhiteMatterInterior(MRI_SURFACE *mris,
                                 MRI *mri_aseg,
                                 MRI *mri_filled,
                                 double resolution,
                                 int wm_val, int gm_val, int csf_val);
int MRISfindClosestWhiteVertex(MRI_SURFACE *mris, float x, float y,
                               float z) ;
int MRISfindClosestVertex(MRI_SURFACE *mris,
                          float x, float y, float z,
                          float *dmin);
double MRIScomputeSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
double MRIScomputeSSEExternal(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                              double *ext_sse) ;
double       MRIScomputeCorrelationError(MRI_SURFACE *mris,
    MRI_SP *mrisp_template, int fno) ;
int          MRISallocExtraGradients(MRI_SURFACE *mris) ;
MRI_SURFACE  *MRISread(const char *fname) ;
MRI_SURFACE  *MRISreadOverAlloc(const char *fname, double pct_over) ;
MRI_SURFACE  *MRISfastRead(const char *fname) ;
int          MRISreadOriginalProperties(MRI_SURFACE *mris,const  char *sname) ;
int          MRISreadCanonicalCoordinates(MRI_SURFACE *mris,
                                          const  char *sname) ;
int          MRISreadInflatedCoordinates(MRI_SURFACE *mris,
                                         const  char *sname) ;
int          MRISreadFlattenedCoordinates(MRI_SURFACE *mris,
                                          const  char *sname) ;
int          MRISreadPialCoordinates(MRI_SURFACE *mris, const char *sname);
int          MRISreadWhiteCoordinates(MRI_SURFACE *mris,const  char *sname) ;
int          MRIScomputeCanonicalCoordinates(MRI_SURFACE *mris) ;
int          MRIScanonicalToWorld(MRI_SURFACE *mris, double phi, double theta,
                                  double *pxw, double *pyw, double *pzw) ;
int          MRISreadPatch(MRI_SURFACE *mris,const  char *pname) ;
int          MRISreadPatchNoRemove(MRI_SURFACE *mris,const  char *pname) ;
int          MRISreadTriangleProperties(MRI_SURFACE *mris,
                                        const  char *mris_fname) ;
int          MRISreadBinaryCurvature(MRI_SURFACE *mris,
                                     const  char *mris_fname) ;
int          MRISreadCurvatureFile(MRI_SURFACE *mris,const char *fname) ;
float        *MRISreadNewCurvatureVector(MRI_SURFACE *mris,
                                         const  char *sname) ;
int          MRISreadNewCurvatureIntoArray(const char *fname,
                                           int in_array_size,
                                           float** out_array) ;
float        *MRISreadCurvatureVector(MRI_SURFACE *mris,const  char *sname) ;
int          MRISreadCurvatureIntoArray(const char *fname,
                                        int in_array_size,
                                        float** out_array) ;
int          MRISreadFloatFile(MRI_SURFACE *mris,const char *fname) ;
#define MRISreadCurvature MRISreadCurvatureFile

MRI *MRISloadSurfVals(const char *srcvalfile,
                      const char *typestring,
                      MRI_SURFACE *Surf,
                      const char *subject,
                      const char *hemi,
                      const char *subjectsdir);
int          MRISreadValues(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadValuesIntoArray(const char *fname,
                                     int in_array_size,
                                     float** out_array) ;
int          MRISreadAnnotation(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteAnnotation(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadAnnotationIntoArray(const char *fname,
                                         int in_array_size,
                                         int** out_array);
int          MRISreadCTABFromAnnotationIfPresent(const char *fname,
                                                 COLOR_TABLE** out_table);
int          MRISisCTABPresentInAnnotation(const char *fname, int* present);
int          MRISreadValuesBak(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadImagValues(MRI_SURFACE *mris,const  char *fname) ;
int          MRIScopyImagValuesToValues(MRI_SURFACE *mris) ;
int          MRIScopyMarksToAnnotation(MRI_SURFACE *mris) ;
int          MRIScopyValuesToImagValues(MRI_SURFACE *mris) ;
int          MRIScopyStatsToValues(MRI_SURFACE *mris) ;

int          MRIScopyValToVal2(MRI_SURFACE *mris) ;
int          MRIScopyValToVal2Bak(MRI_SURFACE *mris) ;
int          MRIScopyValToValBak(MRI_SURFACE *mris) ;
int          MRISsqrtVal(MRI_SURFACE *mris) ;
int          MRISmulVal(MRI_SURFACE *mris, float mul) ;

int          MRISwrite(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteWhiteNormals(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteNormalsAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadNormals(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteNormals(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwritePrincipalDirection(MRI_SURFACE *mris, int dir_index, const  char *fname) ;
int          MRISwriteVTK(MRI_SURFACE *mris,const  char *fname);
int          MRISwriteCurvVTK(MRI_SURFACE *mris, const char *fname);
int          MRISwriteGeo(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteICO(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteSTL(MRI_SURFACE *mris, const char *fname) ;
int          MRISwritePatchAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteDists(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteCurvature(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadNewCurvatureFile(MRI_SURFACE *mris,const  char *fname) ;
int          MRISrectifyCurvature(MRI_SURFACE *mris) ;
#define NORM_NONE  -1
#define NORM_MEAN   0
#define NORM_MEDIAN 1
#define NORM_MAX    2

int          MRISnormalizeCurvature(MRI_SURFACE *mris, int norm_type) ;
int          MRISnormalizeCurvatureVariance(MRI_SURFACE *mris) ;
int          MRISzeroMeanCurvature(MRI_SURFACE *mris) ;
int          MRISnonmaxSuppress(MRI_SURFACE *mris) ;
int          MRISscaleCurvatures(MRI_SURFACE *mris,
                                 float min_curv, float max_curv) ;
int          MRISwriteAreaError(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteAngleError(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwritePatch(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteValues(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteD(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteCurvatureToWFile(MRI_SURFACE *mris, const char *fname) ;
int          MRISwriteTriangleProperties(MRI_SURFACE *mris,
                                         const char *mris_fname);
int          MRISaverageCurvatures(MRI_SURFACE *mris, int navgs) ;
int          MRISminFilterCurvatures(MRI_SURFACE *mris, int niter) ;
int          MRISmaxFilterCurvatures(MRI_SURFACE *mris, int niter) ;
int          MRISaverageMarkedCurvatures(MRI_SURFACE *mris, int navgs) ;
double       MRIScomputeAverageCurvature(MRI_SURFACE *mris, double *psigma) ;
int          MRISaverageVertexPositions(MRI_SURFACE *mris, int navgs) ;
int          MRIScomputeNormal(MRIS *mris, int which, int vno,
                               double *pnx, double *pny, double *pnz) ;

MRI_SURFACE  *MRISoverAlloc(int max_vertices, int max_faces,
                            int nvertices, int nfaces) ;
MRI_SURFACE  *MRISalloc(int nvertices, int nfaces) ;
int          MRISfreeDists(MRI_SURFACE *mris) ;
int          MRISfree(MRI_SURFACE **pmris) ;
int   MRISintegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_avgs);
int   mrisLogIntegrationParms(FILE *fp, MRI_SURFACE *mris,
			      INTEGRATION_PARMS *parms) ;
MRI_SURFACE  *MRISprojectOntoSphere(MRI_SURFACE *mris_src,
                                    MRI_SURFACE *mris_dst, double r) ;
MRI_SURFACE  *MRISprojectOntoEllipsoid(MRI_SURFACE *mris_src,
                                       MRI_SURFACE *mris_dst,
                                       float a, float b, float c) ;
int          MRISsetNeighborhoodSize(MRI_SURFACE *mris, int nsize) ;
int          MRISresetNeighborhoodSize(MRI_SURFACE *mris, int nsize) ;
int          MRISsampleDistances(MRI_SURFACE *mris, int *nbr_count,int n_nbrs);
int          MRISsampleAtEachDistance(MRI_SURFACE *mris, int nbhd_size,
                                      int nbrs_per_distance) ;
int          MRISscaleDistances(MRI_SURFACE *mris, float scale) ;
MRI_SURFACE  *MRISradialProjectOntoEllipsoid(MRI_SURFACE *mris_src,
    MRI_SURFACE *mris_dst,
    float a, float b, float c);
MRI_SURFACE  *MRISclone(MRI_SURFACE *mris_src) ;
MRI_SURFACE  *MRIScenter(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
int MRISorigVertexToVoxel(MRI_SURFACE *,
                          VERTEX *v,
                          MRI *mri,
                          double *pxv, double *pyv, double *pzv) ;
int MRISwhiteVertexToVoxel(MRI_SURFACE *,
                           VERTEX *v,
                           MRI *mri,
                           double *pxv, double *pyv, double *pzv) ;
int          MRISvertexToVoxel(MRI_SURFACE *, VERTEX *v, MRI *mri,
                               double *pxv, double *pyv,
                               double *pzv) ;
int          MRISvertexCoordToVoxel(MRI_SURFACE *, VERTEX *v, MRI *mri,
                                    int coord,
                                    double *pxv, double *pyv,
                                    double *pzv) ;
#if 0
int          MRISworldToTalairachVoxel(MRI_SURFACE *mris, MRI *mri,
                                       double xw, double yw, double zw,
                                       double *pxv, double *pyv, double *pzv) ;
#endif
int          MRISsurfaceRASToVoxel(MRI_SURFACE *mris, MRI *mri, double r, 
                                   double a, double s, 
                                   double *px, double *py, double *pz) ;
int          MRISsurfaceRASToVoxelCached(MRI_SURFACE *mris,
                                         MRI *mri,
                                         double r, double a, double s, 
                                         double *px, double *py, double *pz) ;

// these are the inverse of the previous two
int          MRISsurfaceRASFromVoxel(MRI_SURFACE *mris, MRI *mri, 
                                     double x, double y, double z, 
                                     double *pr, double *pa, double *ps) ;
int          MRISsurfaceRASFromVoxelCached(MRI_SURFACE *mris, MRI *mri, 
                                           double x, double y, double z, 
                                           double *pr, double *pa, double *ps);

int          MRISsurfaceRASToTalairachVoxel(MRI_SURFACE *mris, MRI *mri,
    double xw, double yw, double zw,
    double *pxv, double *pyv, double *pzv) ;

int          MRIStalairachToVertex(MRI_SURFACE *mris,
                                   double xt, double yt, double zt) ;
int           MRIScanonicalToVertex(MRI_SURFACE *mris,
                                    double phi,
                                    double theta) ;
MRI_SURFACE  *MRIStalairachTransform(MRI_SURFACE *mris_src,
                                     MRI_SURFACE *mris_dst);
MRI_SURFACE  *MRISunfold(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                         int max_passes) ;
MRI_SURFACE  *MRISquickSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                              int max_passes) ;
int          MRISregister(MRI_SURFACE *mris,
                          MRI_SP *mrisp_template,
                          INTEGRATION_PARMS *parms,
                          int max_passes,
                          float min_degrees,
                          float max_degrees,
                          int nangles) ;
int          MRISrigidBodyAlignLocal(MRI_SURFACE *mris,
                                     INTEGRATION_PARMS *parms) ;
int          MRISrigidBodyAlignGlobal(MRI_SURFACE *mris,
                                      INTEGRATION_PARMS *parms,
                                      float min_degrees,
                                      float max_degrees,
                                      int nangles) ;
MRI_SURFACE  *MRISunfoldOnSphere(MRI_SURFACE *mris,
                                 INTEGRATION_PARMS *parms,
                                 int max_passes);
MRI_SURFACE  *MRISflatten(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
MRI_SURFACE  *MRISremoveNegativeVertices(MRI_SURFACE *mris,
    INTEGRATION_PARMS *parms,
    int min_neg, float min_neg_pct) ;
int          MRIScomputeFaceAreas(MRI_SURFACE *mris) ;
int          MRISupdateEllipsoidSurface(MRI_SURFACE *mris) ;
MRI_SURFACE  *MRISrotate(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst,
                         float alpha, float beta, float gamma) ;

MRI          *MRISwriteIntoVolume(MRI_SURFACE *mris, MRI *mri, int type) ;
MRI_SURFACE  *MRISreadFromVolume(MRI *mri, MRI_SURFACE *mris) ;



int          MRIScomputeTriangleProperties(MRI_SURFACE *mris) ;
int          MRISpaintVolume(MRI_SURFACE *mris, LTA *lta, MRI *mri) ;
int          MRISsampleStatVolume(MRI_SURFACE *mris, void *sv,int time,
                                  int use_talairach_xform);

int          MRISflattenPatch(MRI_SURFACE *mris) ;
int          MRISflattenPatchRandomly(MRI_SURFACE *mris) ;
int          MRIScomputeMeanCurvature(MRI_SURFACE *mris) ;

int          MRIScomputeEulerNumber(MRI_SURFACE *mris, int *pnvertices,
                                    int *pnfaces, int *pnedges) ;
int          MRIStopologicalDefectIndex(MRI_SURFACE *mris) ;
int          MRISremoveTopologicalDefects(MRI_SURFACE *mris,
                                          float curv_thresh);

int          MRIScomputeSecondFundamentalFormAtVertex(MRI_SURFACE *mris,
                                                      int vno,
                                                      int *vertices,
                                                      int vnum) ;
int          MRIScomputeSecondFundamentalFormThresholded(MRI_SURFACE *mris,
                                                         double thresh) ;
int          MRIScomputeSecondFundamentalForm(MRI_SURFACE *mris) ;
int          MRISuseCurvatureDifference(MRI_SURFACE *mris) ;
int          MRISuseCurvatureStretch(MRI_SURFACE *mris) ;
int          MRISuseCurvatureMax(MRI_SURFACE *mris) ;
int          MRISuseCurvatureMin(MRI_SURFACE *mris) ;
int          MRISuseNegCurvature(MRI_SURFACE *mris) ;
int          MRISuseAreaErrors(MRI_SURFACE *mris) ;
int          MRISuseGaussianCurvature(MRI_SURFACE *mris) ;
int          MRISclearCurvature(MRI_SURFACE *mris) ;
int          MRISclearDistances(MRI_SURFACE *mris) ;
int          MRISusePrincipalCurvature(MRI_SURFACE *mris) ;
int          MRISuseMeanCurvature(MRI_SURFACE *mris) ;
int          MRISuseK1Curvature(MRI_SURFACE *mris);
int          MRISuseK2Curvature(MRI_SURFACE *mris);
int          MRISusePrincipalCurvatureFunction(MRI_SURFACE*		pmris, 
                                               float 			(*f)(float k1, float k2));


int          MRIScomputeCurvatureIndices(MRI_SURFACE *mris,
                                         double *pici, double *pfi) ;
int          MRISuseCurvatureRatio(MRI_SURFACE *mris) ;
int          MRISuseCurvatureContrast(MRI_SURFACE *mris) ;

double MRIScomputeFolding(MRI_SURFACE *mris) ;
double MRISavgVetexRadius(MRIS *Surf, double *StdDev);

int          MRISprojectOntoCylinder(MRI_SURFACE *mris, float radius) ;
double       MRISaverageRadius(MRI_SURFACE *mris) ;
double       MRISaverageCanonicalRadius(MRI_SURFACE *mris) ;
double       MRISmaxRadius(MRI_SURFACE *mris) ;
int          MRISinflateToSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
int          MRISinflateBrain(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
double       MRISrmsTPHeight(MRI_SURFACE *mris) ;
double       MRIStotalVariation(MRI_SURFACE *mris) ;
double       MRIStotalVariationDifference(MRI_SURFACE *mris) ;
double       MRIScurvatureError(MRI_SURFACE *mris, double Kd) ;
MRI_SURFACE  *MRISscaleBrain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst,
                             float scale) ;
int          MRISstoreMetricProperties(MRI_SURFACE *mris) ;
int          MRISrestoreMetricProperties(MRI_SURFACE *mris) ;
double       MRIScomputeAnalyticDistanceError(MRI_SURFACE *mris, int which,
    FILE *fp);
int          MRISzeroNegativeAreas(MRI_SURFACE *mris) ;
int          MRIScountNegativeTriangles(MRI_SURFACE *mris) ;
int          MRIScountMarked(MRI_SURFACE *mris, int mark_threshold) ;
int          MRIScountTotalNeighbors(MRI_SURFACE *mris, int nsize) ;
int          MRISstoreMeanCurvature(MRI_SURFACE *mris) ;
int          MRISreadTetherFile(MRI_SURFACE *mris,
                                const char *fname,
                                float radius) ;
int          MRISreadVertexPositions(MRI_SURFACE *mris,const  char *fname) ;
int          MRISspringTermWithGaussianCurvature(MRI_SURFACE *mris,
                                                 double gaussian_norm,
                                                 double l_spring) ;
int          MRISmarkedSpringTerm(MRI_SURFACE *mris, double l_spring) ;
double       MRISmomentumTimeStep(MRI_SURFACE *mris,
                                  float momentum,
                                  float dt,
                                  float tol,
                                  float n_averages) ;
int          MRISapplyGradient(MRI_SURFACE *mris, double dt) ;
int          MRIScomputeNormals(MRI_SURFACE *mris) ;
int          MRIScomputeSurfaceNormals(MRI_SURFACE *mris,
                                       int which,
                                       int navgs) ;
int          MRIScomputeMetricProperties(MRI_SURFACE *mris) ;
double       MRISrescaleMetricProperties(MRIS *surf);
int          MRISrestoreOldPositions(MRI_SURFACE *mris) ;
int          MRISstoreCurrentPositions(MRI_SURFACE *mris) ;
int          MRISupdateSurface(MRI_SURFACE *mris) ;
double       MRISpercentDistanceError(MRI_SURFACE *mris) ;
int          MRISscaleBrainArea(MRI_SURFACE *mris) ;


int          MRISPsetFrameVal(MRI_SP *mrisp, int frame, float val) ;
MRI_SP       *MRISPcombine(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno);
MRI_SP       *MRISPaccumulate(MRI_SP *mrisp, MRI_SP *mrisp_template, int fno);
int          MRISPcoordinate(MRI_SP *mrisp, float x, float y, float z,
                             int *pu, int *pv) ;
double       MRISPfunctionVal(MRI_SURFACE_PARAMETERIZATION *mrisp,
                              MRI_SURFACE *mris,
                              float x, float y, float z, int fno) ;
MRI_SP       *MRIStoParameterization(MRI_SURFACE *mris, MRI_SP *mrisp,
                                     float scale, int fno) ;
MRI_SURFACE  *MRISfromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris,
                                       int fno) ;
MRI_SURFACE  *MRISnormalizeFromParameterization(MRI_SP *mrisp,
                                                MRI_SURFACE *mris, int fno) ;
MRI_SP       *MRISgradientToParameterization(MRI_SURFACE *mris, MRI_SP *mrisp,
                                             float scale) ;
MRI_SURFACE  *MRISgradientFromParameterization(MRI_SP*mrisp,MRI_SURFACE *mris);

MRI_SP       *MRIScoordsToParameterization(MRI_SURFACE *mris, MRI_SP *mrisp,
                                           float scale, int which_vertices) ;
MRI_SURFACE  *MRIScoordsFromParameterization(MRI_SP *mrisp, MRI_SURFACE *mris, int which_vertices);


float         MRISPsample(MRI_SP *mrisp, float x, float y, float z, int fno) ;
MRI_SP       *MRISPblur(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, float sigma,
                        int fno) ;
MRI_SP       *MRISPconvolveGaussian(MRI_SP *mrisp_src, MRI_SP *mrisp_dst,
                                    float sigma, float radius, int fno) ;
MRI_SP       *MRISPalign(MRI_SP *mrisp_orig, MRI_SP *mrisp_src,
                         MRI_SP *mrisp_tmp, MRI_SP *mrisp_dst) ;
MRI_SP       *MRISPtranslate(MRI_SP *mrisp_src, MRI_SP *mrisp_dst, int du,
                             int dv) ;
MRI_SP       *MRISPclone(MRI_SP *mrisp_src) ;
MRI_SP       *MRISPalloc(float scale, int nfuncs) ;
int          MRISPfree(MRI_SP **pmrisp) ;
MRI_SP       *MRISPread(char *fname) ;
int          MRISPwrite(MRI_SP *mrisp, char *fname) ;

int          MRISwriteArea(MRI_SURFACE *mris,const  char *sname) ;
int          MRISwriteMarked(MRI_SURFACE *mris,const  char *sname) ;
int          MRISmarkedToCurv(MRI_SURFACE *mris) ;

double       MRISParea(MRI_SP *mrisp) ;


#define ORIGINAL_VERTICES   0
#define ORIG_VERTICES       ORIGINAL_VERTICES
#define GOOD_VERTICES       1
#define TMP_VERTICES        2
#define CANONICAL_VERTICES  3
#define CURRENT_VERTICES    4
#define INFLATED_VERTICES   5
#define FLATTENED_VERTICES  6
#define PIAL_VERTICES       7
#define TMP2_VERTICES       8
#define WHITE_VERTICES      9
#define TARGET_VERTICES     10
#define LAYERIV_VERTICES    11
#define LAYER4_VERTICES     LAYERIV_VERTICES

int MRISsaveVertexPositions(MRI_SURFACE *mris, int which) ;
int MRISrestoreVertexPositions(MRI_SURFACE *mris, int which) ;
int MRISrestoreNormals(MRI_SURFACE *mris, int which) ;
int MRISsaveNormals(MRI_SURFACE *mris, int which) ;

/* constants for vertex->tethered */
#define TETHERED_NONE           0
#define TETHERED_FRONTAL_POLE   1
#define TETHERED_OCCIPITAL_POLE 2
#define TETHERED_TEMPORAL_POLE  3

#define MAX_TALAIRACH_Y         100.0f

/* constants for mris->hemisphere */
#define LEFT_HEMISPHERE         0
#define RIGHT_HEMISPHERE        1
#define NO_HEMISPHERE           2
#define BOTH_HEMISPHERES        3

#if 0
#define DEFAULT_A  44.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  70.0f
#else
#define DEFAULT_A  122.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  122.0f
#endif
#define DEFAULT_RADIUS  100.0f

#define MAX_DIM    DEFAULT_B

#if 1
#define DT_INCREASE  1.1 /* 1.03*/
#define DT_DECREASE  0.5
#else
#define DT_INCREASE  1.0 /* 1.03*/
#define DT_DECREASE  1.0
#endif
#define DT_MIN       0.01
#ifdef ERROR_RATIO
#undef ERROR_RATIO
#endif
#define ERROR_RATIO  1.03  /* 1.01 then 1.03 */

#define NO_LABEL     -1


#define WHITE_SURFACE 0
#define GRAY_SURFACE  1
#define GRAY_MID      2

#define REVERSE_X     0
#define REVERSE_Y     1
#define REVERSE_Z     2

int   MRIStranslate(MRI_SURFACE *mris, float dx, float dy, float dz) ;
int   MRISpositionSurface(MRI_SURFACE *mris, MRI *mri_brain,
                          MRI *mri_smooth, INTEGRATION_PARMS *parms);
int   MRISpositionSurfaces(MRI_SURFACE *mris, MRI **mri_flash,
                           int nvolumes, INTEGRATION_PARMS *parms);

int MRISpositionSurface_mef(MRI_SURFACE *mris,
                            MRI *mri_30,
                            MRI *mri_5,
                            INTEGRATION_PARMS *parms,
                            float weight30,
                            float weight5);

int   MRISmoveSurface(MRI_SURFACE *mris, MRI *mri_brain,
                      MRI *mri_smooth, INTEGRATION_PARMS *parms);
int   MRISscaleVals(MRI_SURFACE *mris, float scale) ;
int   MRISsetVals(MRI_SURFACE *mris, float val) ;
int   MRISaverageD(MRI_SURFACE *mris, int navgs) ;
int   MRISgaussianFilterD(MRI_SURFACE *mris, double wt) ;
int   MRISaverageVals(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageVal2baks(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageVal2s(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageMarkedVals(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageMarkedStats(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageEveryOtherVertexPositions(MRI_SURFACE *mris,
                                           int navgs,
                                           int which) ;
int   MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISsoapBubbleOrigVertexPositions(MRI_SURFACE *mris, int navgs) ;
MRI   *MRISwriteSurfaceIntoVolume(MRI_SURFACE *mris, MRI *mri_template,
                                  MRI *mri) ;
#if 0
int   MRISmeasureCorticalThickness(MRI_SURFACE *mris, MRI *mri_brain,
                                   MRI *mri_wm, float nsigma) ;
#else
int   MRISmeasureCorticalThickness(MRI_SURFACE *mris, int nbhd_size,
                                   float max_thickness) ;
#endif

int MRISfindClosestOrigVertices(MRI_SURFACE *mris, int nbhd_size) ;
int MRISfindClosestPialVerticesCanonicalCoords(MRI_SURFACE *mris, int nbhd_size) ;

int   MRISmarkRandomVertices(MRI_SURFACE *mris, float prob_marked) ;
int   MRISmarkNegativeVertices(MRI_SURFACE *mris, int mark) ;
int   MRISripNegativeVertices(MRI_SURFACE *mris) ;
int   MRISclearGradient(MRI_SURFACE *mris) ;
int   MRISclearMarks(MRI_SURFACE *mris) ;
int   MRISclearMark2s(MRI_SURFACE *mris) ;
int   MRISclearFaceMarks(MRI_SURFACE *mris) ;
int   MRISclearFixedValFlags(MRI_SURFACE *mris) ;
int   MRIScopyFixedValFlagsToMarks(MRI_SURFACE *mris) ;
int   MRISclearAnnotations(MRI_SURFACE *mris) ;
int   MRISsetMarks(MRI_SURFACE *mris, int mark) ;
int   MRISsequentialAverageVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISreverseCoords(MRI_SURFACE *mris, int which_direction, int reverse_face_order, int which_coords) ;
int   MRISreverse(MRI_SURFACE *mris, int which, int reverse_face_order) ;
int   MRISreverseFaceOrder(MRIS *mris);
int   MRISdisturbOriginalDistances(MRI_SURFACE *mris, double max_pct) ;
double MRISstoreAnalyticDistances(MRI_SURFACE *mris, int which) ;
int   MRISnegateValues(MRI_SURFACE *mris) ;
int   MRIScopyMeansToValues(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureToValues(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureToImagValues(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureFromValues(MRI_SURFACE *mris) ;
int   MRIScopyVal2ToVal(MRI_SURFACE *mris) ;
int   MRIScopyVal2BakToVal(MRI_SURFACE *mris) ;
int   MRIScopyCurvatureFromImagValues(MRI_SURFACE *mris) ;
int   MRIScopyImaginaryMeansToValues(MRI_SURFACE *mris) ;
int   MRIScopyStandardErrorsToValues(MRI_SURFACE *mris) ;
int   MRIScopyImagValuesToCurvature(MRI_SURFACE *mris) ;
int   MRIScopyValuesToCurvature(MRI_SURFACE *mris) ;
int   MRISaccumulateMeansInVolume(MRI_SURFACE *mris, MRI *mri, int mris_dof,
                                  int mri_dof, int coordinate_system,int sno);
int   MRISaccumulateStandardErrorsInVolume(MRI_SURFACE *mris, MRI *mri,
    int mris_dof, int mri_dof,
    int coordinate_system, int sno) ;
int   MRISaccumulateMeansOnSurface(MRI_SURFACE *mris, int total_dof,
                                   int new_dof);
int   MRISaccumulateImaginaryMeansOnSurface(MRI_SURFACE *mris, int total_dof,
    int new_dof);
int   MRISaccumulateStandardErrorsOnSurface(MRI_SURFACE *mris,
    int total_dof,int new_dof);
#define GRAY_WHITE     1
#define GRAY_CSF       2
int
MRIScomputeMaxGradBorderValuesPial(MRI_SURFACE *mris,MRI *mri_brain,
                                   MRI *mri_smooth, double sigma,
                                   float max_thickness,
                                   float dir, FILE *log_fp,
                                   int callno, MRI *mri_mask) ;
int
MRIScomputeMaxGradBorderValues(MRI_SURFACE *mris,MRI *mri_brain,
                               MRI *mri_smooth, double sigma,
                               float max_thickness, float dir, FILE *log_fp,
                               MRI *mri_wm, int callno) ;
int MRIScomputeInvertedGrayWhiteBorderValues(MRI_SURFACE *mris,
                                             MRI *mri_brain,
                                             MRI *mri_smooth,
                                             double inside_hi,
                                             double border_hi,
                                             double border_low,
                                             double outside_low,
                                             double outside_hi,
                                             double sigma,
                                             float max_thickness,
                                             FILE *log_fp);
int MRIScomputeInvertedPialBorderValues(MRI_SURFACE *mris,
                                        MRI *mri_brain,
                                        MRI *mri_smooth,
                                        double inside_hi,
                                        double border_hi,
                                        double border_low,
                                        double outside_low,
                                        double outside_hi,
                                        double sigma,
                                        float max_thickness,
                                        FILE *log_fp);
int   MRIScomputeBorderValues(MRI_SURFACE *mris,
                              MRI *mri_brain,
                              MRI *mri_smooth,
                              double inside_hi,
                              double border_hi,
                              double border_low,
                              double outside_low,
                              double outside_hi,
                              double sigma,
                              float max_dist,
                              FILE *log_fp,
                              int white,
                              MRI *mri_mask, double thresh);
int  MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                                   MRI *mri_smooth);
int  MRIScomputeGraySurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                                  MRI *mri_smooth, float gray_surface);
int  MRIScomputeDistanceErrors(MRI_SURFACE *mris, int nbhd_size,int max_nbrs);
int  MRISsetAllMarks(MRI_SURFACE *mris, int mark) ;
int  MRISscaleCurvature(MRI_SURFACE *mris, float scale) ;
int  MRISwriteTriangularSurface(MRI_SURFACE *mris,const char *fname) ;
int  MRISripFaces(MRI_SURFACE *mris) ;
int  MRISremoveRipped(MRI_SURFACE *mris) ;
int  MRISbuildFileName(MRI_SURFACE *mris, const char *sname, char *fname) ;
int  MRISsmoothSurfaceNormals(MRI_SURFACE *mris, int niter) ;
int  MRISsoapBubbleD(MRI_SURFACE *mris, int niter) ;
int  MRISsoapBubbleVals(MRI_SURFACE *mris, int niter) ;
int  MRISmodeFilterVals(MRI_SURFACE *mris, int niter) ;
int  MRISmodeFilterAnnotations(MRI_SURFACE *mris, int niter) ;
int  MRISmodeFilterZeroVals(MRI_SURFACE *mris) ;
int  MRISreadBinaryAreas(MRI_SURFACE *mris,const  char *mris_fname) ;
int  MRISwriteAreaErrorToValFile(MRI_SURFACE *mris,const  char *name) ;
int  MRIStransform(MRI_SURFACE *mris, MRI *mri,
                   TRANSFORM *transform, MRI *mri_dst) ;
int  MRISmatrixMultiply(MRIS *mris, MATRIX *M);

int  MRISanisotropicScale(MRI_SURFACE *mris, float sx, float sy, float sz) ;
double MRIScomputeVertexSpacingStats(MRI_SURFACE *mris, double *psigma,
                                     double *pmin, double *pmax, int *pvno,
                                     int *pvno2, int which_vertices);
double MRIScomputeTotalVertexSpacingStats(MRI_SURFACE *mris, double *psigma,
                                          double *pmin, double *pmax,
                                          int *pvno,
                                          int *pvno2);
double MRIScomputeFaceAreaStats(MRI_SURFACE *mris, double *psigma,
                                double *pmin, double *pmax);
int MRISprintTessellationStats(MRI_SURFACE *mris, FILE *fp) ;
int MRISprintVertexStats(MRI_SURFACE *mris, int vno, FILE *fp, int which_vertices) ;
int MRISmergeIcosahedrons(MRI_SURFACE *mri_src, MRI_SURFACE *mri_dst) ;
int MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico) ;

////////////////////////////////////////////////////
/* for mris_topo_fixer */
typedef struct
{
  MRIS *mris;
  int n_vertices,n_faces; //number of vertices & faces before top. correction
  int *vtrans_to,*ftrans_to,*vtrans_from,*ftrans_from;
  MRIS *mris_source;
}
MRI_PATCH, MRIP;

/* Different verbose mode for mris_fix_topology and mris_topo_fix */
#define VERBOSE_MODE_DEFAULT 1
#define VERBOSE_MODE_LOW 2
#define VERBOSE_MODE_MEDIUM 3
#define VERBOSE_MODE_HIGH 4

#include "histo.h" // HISTOGRAM

typedef struct
{
  int verbose; // verbose mode
  int smooth;  // smoothing defect
  int match;   // using local intensity to match surface
  int defect_number; // the defect_number
  int mode; // which mode to use (not used so far)
  int no_self_intersections; //to prevent self-intersection
  int contrast; //direction of the contrast
  int write; //writing out temporary surfaces using fname
  char fname[STRLEN] ;
  double nattempts_percent;
  int nattempts;
  int minimal_mode;
  int nminattempts;
  double minimal_loop_percent;
  int nminimal_attempts;
  int max_face;
  void* patchdisk; // the patching surfaces
  // loglikelihood
  double  l_mri ;
  double  l_curv ;
  double  l_qcurv;
  double  l_unmri ;
  int volume_resolution; /* used if l_unmri is on */
  MRI *mri; //brain volume
  MRI *mri_wm; //wm volume
  HISTOGRAM *h_k1, *h_k2,*h_gray,*h_white,*h_dot,*h_border, *h_grad;
  MRI *mri_gray_white, *mri_k1_k2;
  MATRIX *transformation_matrix;
  //defect info
  MRIP *mrip;
  MRIS *mris_defect;
  void *defect_list;
  void   *dp;
  //statistics
  float fitness,initial_fitness;
  int ngeneratedpatches,nselfintersectingpatches;
}
TOPOFIX_PARMS;

void MRIScomputeInitialFitness(MRIS *mris, TOPOFIX_PARMS *parms);
void MRIScomputeDistanceVolume(TOPOFIX_PARMS *parms,
                               float distance_to_surface);
void MRISsaveLocal(MRIS *mris, TOPOFIX_PARMS *parms, char *name); //floflo
int MRISisSurfaceValid(MRIS *mris, int patch,int verbose);
void MRISinitTopoFixParameters(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISinitDefectParameters(MRIS *mris, TOPOFIX_PARMS *parms);
void TOPOFIXfreeDP(TOPOFIX_PARMS *parms);
void MRISinitDefectPatch(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISdefectMatch(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISprintInfo(TOPOFIX_PARMS *parms);
double MRIScomputeFitness(MRIS* mris,TOPOFIX_PARMS *parms,int verbose);
int IsMRISselfIntersecting(MRI_SURFACE *mris);
void MRISmapOntoSphere(MRIS *mris);
void MRISidentifyDefects(MRIS *mris, TOPOFIX_PARMS *parms);
void MRISmarkPatchVertices(MRIS *mris,
                           TOPOFIX_PARMS *parms,
                           int ninitvertices);
void MRISmarkBorderVertices(MRIS *mris,
                            TOPOFIX_PARMS *parms,
                            int mark);

#define GREEDY_SEARCH 0
#define GENETIC_SEARCH 1
#define RANDOM_SEARCH 2

typedef struct
{
  int     max_patches ;
  int     max_unchanged ;
  int     niters;  /* stop the genetic algorithm after n iterations */
  int     genetic; /* to use the genetic algorithm */
  int     search_mode; /* which topology correction (greedy,genetic,random)*/
  int     retessellation_mode; /* which retessellation (ordering dependent)*/
  int     vertex_eliminate; /* kill less used vertices */
  int     initial_selection; /* find set of initial chromosomes */
  double  l_mri ;
  double  l_curv ;
  double  l_qcurv;
  double  l_unmri ;
  int volume_resolution; /* used if l_unmri is on */
  int keep; /* keep every vertex in the defect */
  int edge_table; /*using edge table or not */
  char *save_fname; /* save results into folder name save_fname */
  int defect_number; /* save only for one specific defect */
  int verbose; /* outputs information */
  int smooth; /* smoothing the patch */
  int match; /* matching the patch onto the surface using local parameters */
  int movie; /* save interesting files for movie purpose */
  int correct_defect; /* correct only one single defect */
  int check_surface_intersection; /* check if self-intersection happens */
  int optimal_mapping; /* find optmal mapping by genrating sevral mappings */
}
TOPOLOGY_PARMS ;

int MRIScenterSphere(MRI_SURFACE *mris);
int MRISmarkOrientationChanges(MRI_SURFACE *mris);
MRIS* MRISextractMainComponent(MRI_SURFACE *mris,
                               int do_not_extract,
                               int verbose,
                               int *ncpts);
MRIS* MRISextractMarkedVertices(MRIS *mris);
MRIS* MRISremoveRippedSurfaceElements(MRIS *mris);

MRI_SURFACE *MRIScorrectTopology(MRI_SURFACE *mris,
                                 MRI_SURFACE *mris_corrected,
                                 MRI *mri,
                                 MRI *mri_wm,
                                 int nsmooth,
                                 TOPOLOGY_PARMS *parms) ;
int MRISsmoothOnSphere(MRIS* mris, int niters);
int mrisCountIntersectingFaces(MRIS *mris, int*flist , int nfaces);
int MRISripDefectiveFaces(MRI_SURFACE *mris) ;
int MRISunrip(MRI_SURFACE *mris) ;
int MRISdivideLongEdges(MRI_SURFACE *mris, double thresh) ;
int MRISdivideEdges(MRI_SURFACE *mris, int npoints) ;
int MRISremoveTriangleLinks(MRI_SURFACE *mris) ;
int MRISsetOriginalFileName(char *orig_name) ;
int MRISsetSulcFileName(const char *sulc_name) ;
int MRISsetCurvatureName(int nth, char *name);
int MRISprintCurvatureNames(FILE *fp);
int MRISsetInflatedFileName(char *inflated_name) ;
int MRISsetRegistrationSigmas(float *sigmas, int nsigmas) ;

int MRISextractVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices) ;
int MRISimporttVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices) ;
int MRISextractCurvatureVector(MRI_SURFACE *mris, float *curvs) ;
int MRISextractCurvatureDoubleVector(MRI_SURFACE *mris, double *curvs) ;
int MRISextractCurvatureVectorDouble(MRI_SURFACE *mris, double *curvs, int offset) ;
#define MRISexportCurvatureVector  MRISextractCurvatureVector

int MRISimportCurvatureVector(MRI_SURFACE *mris, float *curvs) ;
int MRISimportValVector(MRI_SURFACE *mris, float *vals) ;
int MRISexportValVector(MRI_SURFACE *mris, float *vals) ;
int MRISexportValVectorDouble(MRI_SURFACE *mris, double *vals, int offset) ;
int MRISimportValFromMatrixColumn(MRI_SURFACE *mris, MATRIX *m, int col) ;


/* multi-timepoint (or stc) files */
int  MRISwriteStc(char *fname, MATRIX *m_data, float epoch_begin_lat,
                  float sample_period, int *vertices) ;

#define MAX_LINKS 10
typedef struct
{
  double     tx, ty, tz ;   // target coords
  int        l_in ;         // inside label
  int        l_out ;        // outside label
  HISTOGRAM *h_in ;         // inside label histogram
  HISTOGRAM *h_out ;        // inside label histogram
  double     mag ;          // directional derivative in normal dir initially
  int        linked_vno[MAX_LINKS] ;   // is it the same as another vertex
  int        linked_sno[MAX_LINKS] ;   // surface that it's linked to
  int        nlinks ;
  double     p ;            // p-value for user to fill in
}
VERTEX_INFO ;

typedef struct
{
  MRI_SURFACE *mris[MAX_SURFACES] ;
  int         nsurfaces ;
  int         nvertices ;   // sum of individual surface nvertices
  int         nfaces ;      // sum of individual surface faces
  int         labels[MAX_SURFACES] ;  // interior label of the structure
  MRI_SURFACE *mris_total ; // single surface with all other surfs in it
  VERTEX_INFO *vi ;                   // one/vertex in mris_total
  int         vstart[MAX_SURFACES] ;  /* starting index of this
                                         surface in mris_total vertices */
}
MRI_SURFACE_ARRAY, MSA ;
#if 1
#include "mrishash.h"
float  MRISdistanceToSurface(MRI_SURFACE *mris, MHT *mht,
                             float x0, float y0, float z0,
                             float nx, float ny, float nz) ;
int    MRISexpandSurface(MRI_SURFACE *mris,
                         float distance,
                         INTEGRATION_PARMS *parms, int use_thickness, int nsurfaces) ;
int MRISripZeroThicknessRegions(MRI_SURFACE *mris) ;

#endif

/* cortical ribbon */
MRI   *MRISribbon(MRI_SURFACE *inner_mris,
                  MRI_SURFACE *outer_mris,
                  MRI *mri_src,
                  MRI *mri_dst);
MRI   *MRISaccentuate(MRI *mri_src,
                      MRI *mri_dst,
                      int lo_thresh,
                      int hi_thresh);
MRI *MRISfillInterior(MRI_SURFACE *mris,
                      double resolution,
                      MRI *mri_interior) ;
MRI   *MRISshell(MRI *mri_src,
                 MRI_SURFACE *mris,
                 MRI *mri_dst,
                 int clearflag);
MRI   *MRISfloodoutside(MRI *mri_src,MRI *mri_dst);
// src never used. dst is modified

/* high resolution cortical ribbon */
MRI   *MRISpartialribbon(MRI_SURFACE *inner_mris_lh,
                         MRI_SURFACE *outer_mris_lh,
                         MRI_SURFACE *inner_mris_rh,
                         MRI_SURFACE *outer_mris_rh,
                         MRI *mri_src,
                         MRI *mri_dst,
                         MRI *mri_mask);
MRI   *MRISpartialaccentuate(MRI *mri_src,
                             MRI *mri_dst,
                             int lo_thresh,
                             int hi_thresh);
MRI   *MRISpartialshell(MRI *mri_src,
                        MRI_SURFACE *mris,
                        MRI *mri_dst,
                        int clearflag);
MRI   *MRISpartialfloodoutside(MRI *mri_src,
                               MRI *mri_dst);
// src is used. dst must be initialized

#define MRIS_BINARY_QUADRANGLE_FILE    0    /* homegrown */
#define MRIS_ASCII_TRIANGLE_FILE       1    /* homegrown */
#define MRIS_ASCII_FILE MRIS_ASCII_TRIANGLE_FILE
#define MRIS_GEO_TRIANGLE_FILE         2    /* movie.byu format */
#define MRIS_ICO_SURFACE               3
#define MRIS_TRIANGULAR_SURFACE        MRIS_ICO_SURFACE
#define MRIS_ICO_FILE                  4
#define MRIS_VTK_FILE                  5
#define MRIS_STL_FILE                  6
#define MRIS_GIFTI_FILE                7
#define MRIS_ANNOT_FILE                8
#define MRIS_VOLUME_FILE               9


#define IS_QUADRANGULAR(mris)                   \
  (mris->type == MRIS_BINARY_QUADRANGLE_FILE)


/* actual constants are in mri.h */
#define RH_LABEL           MRI_RIGHT_HEMISPHERE
#define RH_LABEL2          MRI_RIGHT_HEMISPHERE2
#define LH_LABEL           MRI_LEFT_HEMISPHERE

#define MRISPvox(m,u,v)   (*IMAGEFpix(m->Ip,u,v))


/* VECTORIAL_REGISTRATION*/
int MRISPfunctionVectorVals(MRI_SURFACE_PARAMETERIZATION *mrisp,
                            MRI_SURFACE *mris,
                            float x, float y, float z,
                            int *frames, int nframes, double *vals);
MRI_SURFACE * MRISfromParameterizations(MRI_SP *mrisp,
                                        MRI_SURFACE *mris, int *frames,
                                        int *indices,int nframes);
MRI_SP * MRIStoParameterizations(MRI_SURFACE *mris,
                                 MRI_SP *mrisp, float scale,int *frames,
                                 int *indices,int nframes);
MRI_SP * MRISPblurFrames(MRI_SP *mrisp_src,
                         MRI_SP *mrisp_dst, float sigma, int *frames,
                         int nframes);

typedef struct
{
  float   x ;
  float   y ;
  float   z ;
}
SMALL_VERTEX ;

typedef struct
{
  int   nvertices ;
  SMALL_VERTEX *vertices ;
}
SMALL_SURFACE ;


int   MRISpositionSurfaceArray(MRI_SURFACE_ARRAY *msa, MRI *mri_brain,
                               MRI *mri_smooth, INTEGRATION_PARMS *parms);

#define MRISSread  MRISreadVerticesOnly
SMALL_SURFACE *MRISreadVerticesOnly(char *fname) ;
int           MRISSfree(SMALL_SURFACE **pmriss) ;


int  MRISsubsampleDist(MRI_SURFACE *mris, float spacing) ;
int  MRISwriteDecimation(MRI_SURFACE *mris, char *fname) ;
int  MRISreadDecimation(MRI_SURFACE *mris, char *fname) ;


#define VERTEX_COORDS      0
#define VERTEX_VALS        1
#define VERTEX_VAL         VERTEX_VALS
#define VERTEX_AREA        2
#define VERTEX_CURV        3
#define VERTEX_CURVATURE   VERTEX_CURV
#define VERTEX_LABEL       4
#define VERTEX_ANNOTATION  5
#define VERTEX_DX          6
#define VERTEX_DY          7
#define VERTEX_DZ          8
#define VERTEX_STATS       9
#define VERTEX_LOGODDS     10


int MRISclearOrigArea(MRI_SURFACE *mris) ;
int MRISclearOrigDistances(MRI_SURFACE *mris) ;
int MRIScombine(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total,
                MRIS_HASH_TABLE *mht, int which) ;
int MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total,
                      MRIS_HASH_TABLE *mht, int which) ;
int   MRISorigAreaToCurv(MRI_SURFACE *mris) ;
int   MRISareaToCurv(MRI_SURFACE *mris) ;
int   MRISclear(MRI_SURFACE *mris, int which) ;
int   MRISnormalize(MRI_SURFACE *mris, int dof, int which) ;

int  MRIScopyMRI(MRIS *Surf, MRI *Src, int Frame, char *Field);
MRI *MRIcopyMRIS(MRI *mri, MRIS *surf, int Frame, char *Field);

MRI *MRISsmoothMRI(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask, MRI *Targ);
MRI *MRISsmoothMRIFast(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask,  MRI *Targ);
MRI *MRISsmoothMRIFastD(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *IncMask,  MRI *Targ);
int MRISsmoothMRIFastCheck(int nSmoothSteps);
int MRISsmoothMRIFastFrame(MRIS *Surf, MRI *Src, int frame, int nSmoothSteps, MRI *IncMask);


int  MRISclearFlags(MRI_SURFACE *mris, int flags) ;
int  MRISsetCurvature(MRI_SURFACE *mris, float val) ;
int  MRISsetFlags(MRI_SURFACE *mris, int flags) ;

int MRISgaussianFilterD(MRI_SURFACE *mris, double wt) ;
int MRISmedianFilterD(MRI_SURFACE *mris, int nmedians, int vtotal) ;
int MRISmedianFilterVals(MRI_SURFACE *mris, int nmedians) ;
int MRISmedianFilterVal2s(MRI_SURFACE *mris, int nmedians) ;
int MRISmedianFilterVal2baks(MRI_SURFACE *mris, int nmedians) ;
int MRISmedianFilterCurvature(MRI_SURFACE *mris, int nmedians);
int MRISfileNameType(const char *fname) ;
unsigned long MRISeraseOutsideOfSurface(float h,
                                        MRI* mri_dst,
                                        MRIS *mris,
                                        unsigned char val) ;

/* Some utility functions to handle reading and writing annotation
   values. MRISRGBToAnnot stuffs an r,g,b tuple into an annotation
   value and MRISAnnotToRGB separates an annotation value into an
   r,g,b tuple. */
#define MRISAnnotToRGB(annot,r,g,b)             \
  r = annot & 0xff ;                            \
  g = (annot >> 8) & 0xff ;                     \
  b = (annot >> 16) & 0xff ;
#define MRISRGBToAnnot(r,g,b,annot)                                     \
  annot = ((r) & 0xff) | (((g) & 0xff) << 8) | (((b) & 0xff) << 16);

int MRISextendedNeighbors(MRIS *SphSurf,int TargVtxNo, int CurVtxNo,
                          double DotProdThresh, int *XNbrVtxNo,
                          double *XNbrDotProd, int *nXNbrs,
                          int nXNbrsMax, int DistType);
MRI *MRISgaussianSmooth(MRIS *Surf, MRI *Src, double GStd, MRI *Targ,
                        double TruncFactor);
MRI *MRISdistSphere(MRIS *surf, double dmax);
int MRISgaussianWeights(MRIS *surf, MRI *dist, double GStd);
MRI *MRISspatialFilter(MRI *Src, MRI *wdist, MRI *Targ);

MATRIX *surfaceRASToSurfaceRAS_(MRI *src, MRI *dst, LTA *lta);

int MRISsurf2surf(MRIS *mris, MRI *dst, LTA *lta);
// convert all vertex positions
int MRISsurf2surfAll(MRIS *mris, MRI *dst, LTA *lta);
void MRISsetReadFrame(int frame);
int MRISaddCommandLine(MRI_SURFACE *mris, char *cmdline) ;
int MRISgetReadFrame(void);
int MRISabsCurvature(MRI_SURFACE *mris) ;
int MRISabsVals(MRI_SURFACE *mris) ;
int MRISsmoothFrames(MRI_SURFACE *mris, MRI *mri, int navgs) ;
int MRISwriteFrameToValues(MRI_SURFACE *mris, MRI *mri, int frame) ;
int MRISreadFrameFromValues(MRI_SURFACE *mris, MRI *mri, int frame) ;
MRI *MRISar1(MRIS *surf, MRI *src, MRI *mask, MRI *ar1);
int **MRIScrsLUT(MRIS *surf, MRI *src);
int MRIScrsLUTFree(int **crslut);
int MRISremoveOverlapWithSmoothing(MRI_SURFACE *mris,
                                   INTEGRATION_PARMS *parms) ;
int MRISupsampleIco(MRI_SURFACE *mris, MRI_SURFACE *mris_new) ;
int MRIScopyVolGeomFromMRI(MRI_SURFACE *mris, MRI *mri) ;

// the next function assumes the surface vertex positions are 
//    set in the voxel coordinates of srcMri
// and it constructs a valid RAS for the surface
void MRISsetVolumeForSurface(MRI_SURFACE *mris, MRI *srcMri);

MRI *MRISremoveRippedFromMask(MRIS *surf, MRI *mask, MRI *outmask);
int MRISremoveIntersections(MRI_SURFACE *mris) ;
int MRIScopyMarkedToMarked2(MRI_SURFACE *mris) ;
int MRIScopyMarked2ToMarked(MRI_SURFACE *mris) ;
int MRIScopyMarkedToMarked3(MRI_SURFACE *mris) ;
int MRIScopyMarked3ToMarked(MRI_SURFACE *mris) ;
int MRISexpandMarked(MRI_SURFACE *mris) ;
double MRISsmoothingArea(MRIS *mris, int vtxno, int niters);

int MRISopenMarked(MRI_SURFACE *mris, int order) ;
int MRIScloseMarked(MRI_SURFACE *mris, int order) ;
int MRISerodeMarked(MRI_SURFACE *mris, int ndil) ;
int MRISdilateMarked(MRI_SURFACE *mris, int ndil) ;
int MRISerodeRipped(MRI_SURFACE *mris, int ndil) ;
int MRISdilateRipped(MRI_SURFACE *mris, int ndil) ;
MRI *MRISdilateConfined(MRIS *surf,
                        MRI *mask,
                        int annotidmask,
                        int niters,
                        int newid);
MRI *MRISfbirnMask_SFG_Cing(MRIS *surf);
MRI *MRISfbirnMask_MOF_RACing(MRIS *surf);

int   MRISvalidVertices(MRI_SURFACE *mris) ;
int MRISmarkedVertices(MRI_SURFACE *mris) ;

int MRIScomputeClassStatistics(MRI_SURFACE *mris,
                               MRI *mri,
                               float *pwhite_mean,
                               float *pwhite_std,
                               float *pgray_mean,
                               float *pgray_std) ;
int MRIScomputeClassModes(MRI_SURFACE *mris,
                          MRI *mri,
                          float *pwhite_mode,
                          float *pgray_mode,
                          float *pcsf_mode);
int MRISrasToVoxel(MRI_SURFACE *mris,
                   MRI *mri,
                   double xs, double ys, double zs,
                   double *pxv, double *pyv, double *pzv) ;
int MRISrestoreRipFlags(MRI_SURFACE *mris) ;
int MRISstoreRipFlags(MRI_SURFACE *mris) ;
int MRISripMedialWall(MRI_SURFACE *mris) ;
int MRISripMarked(MRI_SURFACE *mris) ;
int MRISripUnmarked(MRI_SURFACE *mris) ;
int MRISzeroMedialWallCurvature(MRI_SURFACE *mris) ;
int MRISvertexNormalInVoxelCoords(MRI_SURFACE *mris,
                                  MRI *mri,
                                  int vno,
                                  double *pnx, double *pny, double *pnz) ;

#define MRISgetCoords(v,c,vx,vy,vz) \
 switch(c) { \
   case ORIGINAL_VERTICES:  (*vx) = (v)->origx;  (*vy) = (v)->origy;  (*vz) = (v)->origz; break; \
   case TMP_VERTICES:       (*vx) = (v)->tx;     (*vy) = (v)->ty;     (*vz) = (v)->tz; break; \
   case CANONICAL_VERTICES: (*vx) = (v)->cx;     (*vy) = (v)->cy;     (*vz) = (v)->cz; break; \
   case CURRENT_VERTICES:   (*vx) = (v)->x;      (*vy) = (v)->y;      (*vz) = (v)->z; break; \
   case TARGET_VERTICES:   (*vx) = (v)->targx;   (*vy) = (v)->targy;  (*vz) = (v)->targz; break; \
   case INFLATED_VERTICES:  (*vx) = (v)->infx;   (*vy) = (v)->infy;   (*vz) = (v)->infz; break; \
   case FLATTENED_VERTICES: (*vx) = (v)->fx;     (*vy) = (v)->fy;     (*vz) = (v)->fz; break; \
   case PIAL_VERTICES:      (*vx) = (v)->pialx;  (*vy) = (v)->pialy;  (*vz) = (v)->pialz; break; \
   case TMP2_VERTICES:      (*vx) = (v)->tx2;    (*vy) = (v)->ty2;    (*vz) = (v)->tz2; break; \
   case WHITE_VERTICES:     (*vx) = (v)->whitex; (*vy) = (v)->whitey; (*vz) = (v)->whitez; break; \
   default: break; \
 }

#include "label.h" // LABEL
int MRISlogOdds(MRI_SURFACE *mris, LABEL *area, double slope)  ;
MRI_SP  *MRISPorLabel(MRI_SP *mrisp, MRI_SURFACE *mris, LABEL *area) ;
MRI_SP  *MRISPandLabel(MRI_SP *mrisp, MRI_SURFACE *mris, LABEL *area) ;
MRI *MRISlabel2Mask(MRIS *surf, LABEL *lb, MRI *mask);
int   MRIScomputeAverageCircularPhaseGradient(MRI_SURFACE *mris,
    LABEL *area,
    float *pdx,
    float *pdy,
    float *pdz);

int MRISmaskLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISmaskNotLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISripLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISripNotLabel(MRI_SURFACE *mris, LABEL *area) ;
int MRISsegmentMarked(MRI_SURFACE *mris,
                      LABEL ***plabel_array,
                      int *pnlabels,
                      float min_label_area) ;
int MRISsegmentAnnotated(MRI_SURFACE *mris,
                         LABEL ***plabel_array,
                         int *pnlabels,
                         float min_label_area) ;

int MRISaverageGradients(MRI_SURFACE *mris, int num_avgs) ;
int MRISaverageGradientsFast(MRI_SURFACE *mris, int num_avgs);
int MRISaverageGradientsFastCheck(int num_avgs);

int MRISnormalTermWithGaussianCurvature(MRI_SURFACE *mris,double l_lambda) ;
int MRISnormalSpringTermWithGaussianCurvature(MRI_SURFACE *mris,
                                              double gaussian_norm,
                                              double l_spring) ;
int MRISmakeDensityMap(MRI_SURFACE *mris,
                       double resolution,
                       double radius,
                       int diag_no,
                       MRI **pmri);
double MRIScomputeWhiteVolume(MRI_SURFACE *mris,
                              MRI *mri_aseg,
                              double resolution);
int MRIShistoThresholdCurvature(MRI_SURFACE *mris, float thresh_pct);
int MRISthresholdCurvature(MRI_SURFACE *mris, float thresh, int use_abs);
int MRISbinarizeCurvature(MRI_SURFACE *mris,
                          float thresh,
                          float low,
                          float high,
                          int use_abs);
int MRISsetVal2(MRI_SURFACE *mris, float val);
MRI *MRIScomputeDistanceToSurface(MRI_SURFACE *mris,
                                  MRI *mri,
                                  float resolution) ;
int MRISdistanceTransform(MRI_SURFACE *mris,LABEL *area, int mode) ;
int MRISinvertMarks(MRI_SURFACE *mris) ;


/* value field labels for swap_vertex_fields(), also used in tcl script */
#define FIELD_CURV      0
#define FIELD_CURVBAK   1
#define FIELD_VAL       2
#define FIELD_IMAG_VAL  3
#define FIELD_VAL2      4
#define FIELD_VALBAK    5
#define FIELD_VAL2BAK   6
#define FIELD_STAT      7

#define SCLV_VAL          0
#define SCLV_VAL2         1
#define SCLV_VALBAK       2
#define SCLV_VAL2BAK      3
#define SCLV_VALSTAT      4
#define SCLV_IMAG_VAL     5
#define SCLV_MEAN         6
#define SCLV_MEAN_IMAG    7
#define SCLV_STD_ERROR    8
#define NUM_SCALAR_VALUES 9
HISTOGRAM *MRISgetHistogram(MRI_SURFACE *mris, int nbins, int field);

// Discrete Principal Curvature and Related  vvvvvvvvvvvvvvvvvv
// Column / width for formatted output
#define		G_LC				50
#define		G_RC				30

int		slprints(
    char*		apch_txt
);

void	cprintf(
	char*		apch_left,
	float		af_right
);

void	cprints(
	char*		apch_left,
	char*		apch_right
);

void	cprintd(
	char*		apch_left,
	int		a_right
);

short	FACES_aroundVertex_reorder(
	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			pv_geometricOrder
);

short	FACE_vertexIndex_find(
    	FACE*			pFace,
    	int 			avertex 
);

short	VECTOR_elementIndex_findNotEqual(
	VECTOR*			apV,
	float			af_searchTerm
);

short	VECTOR_elementIndex_find(
	VECTOR*			apV,
	float			af_searchTerm
);

short	MRIS_vertexProgress_print(
    	MRIS*			apmris,
    	int			avertex,
    	char*			apch_message
);

int	FACE_vertexIndexAtMask_find(
	FACE*			apFACE_I,
	VECTOR*			apv_verticesCommon
);

short	VERTICES_commonInFaces_find(
	FACE*			apFACE_I,
	FACE*			apFACE_J,
	VECTOR*			apv_verticesCommon
);

short	FACES_Hcurvature_determineSign(
    	MRIS*			apmris,
    	FACE*			apFACE_O,
    	FACE*			apFACE_I
);

int	VERTEX_faceAngles_determine(
    	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			apv_angle
);

int	VERTEX_faceMinMaxAngles_determine(
    	MRIS*			apmris,
    	int			avertex,
    	int*			ap_minIndex,
    	float*			apf_minAngle,
    	int*			ap_maxIndex,
    	float*			apf_maxAngle
);

int	MRIS_facesAtVertices_reorder(
    	MRIS*			apmris
);

int	MRIScomputeGeometricProperties(
    	MRIS*			apmris
);

short	FACES_aroundVertex_reorder(
    	MRIS*			apmris,
    	int			avertex,
    	VECTOR*			pv_geometricOrder
);

float	FACES_angleNormal_find(
    	MRIS*			apmris,
    	FACE*			apFACE_I,
    	FACE*			apFACE_J
);

float	FACES_commonEdgeLength_find(
    	MRIS*			apmris,
    	FACE*			apFACE_I,
    	FACE*			apFACE_J
);

short	MRIS_discreteKH_compute(
	MRIS*			apmris
);

short	MRIS_discretek1k2_compute(
	MRIS*			apmris,
	short			ab_signedPrincipals
);

// The discrete curvature calculations are based on the Gauss-Bonnet Scheme.
// 
// For 'K' and 'H', see
// 
// @article{1280456,
//  author = {Evgeni Magid and Octavian Soldea and Ehud Rivlin},
//  title = {A comparison of Gaussian and mean curvature estimation 
// 	  methods on triangular meshes of range image data},
//  journal = {Comput. Vis. Image Underst.},
//  volume = {107},
//  number = {3},
//  year = {2007},
//  issn = {1077-3142},
//  pages = {139--159},
//  doi = {http://dx.doi.org/10.1016/j.cviu.2006.09.007},
//  publisher = {Elsevier Science Inc.},
//  address = {New York, NY, USA},
//  }
// 
// and for k1 and k2 see:
// 
// @misc{ meyer02discrete,
//   author = "M. MEYER and M. DESBRUN and P. SCHR and A. BARR",
//   title = "Discrete DifferentialGeometry Operators for Triangulated 2-Manifolds",
//   text = "MEYER, M., DESBRUN, M., SCHR ODER, P., AND BARR, 
// 	  A. H. Discrete DifferentialGeometry
//     	  Operators for Triangulated 2-Manifolds, 2002. VisMath.",
//   year = "2002",
//   url = "citeseer.ist.psu.edu/meyer02discrete.html" }
// 

short	MRIScomputeSecondFundamentalFormDiscrete(
	MRIS*			apmris,
	short			ab_signedPrincipals
);

int  	MRISminMaxCurvatureIndicesLookup(
  	MRI_SURFACE*  		apmris,
  	int*   			ap_vertexMin,
  	int*   			ap_vertexMax
);

int  	MRISvertexCurvature_set(
  	MRI_SURFACE*  		apmris,
  	int   			aindex,
  	float   		af_val
);

int	MRISzeroCurvature(
  	MRI_SURFACE*  		apmris
);

int	MRISuseK1Curvature(
  	MRI_SURFACE*  		mris
);

int	MRISuseK2Curvature(
  	MRI_SURFACE*  		mris
);

// end Discrete Principal Curvature and Related ^^^^^^^^^^^^^^^^^^

void UpdateMRIS(MRI_SURFACE *mris,const char *fname);
int  MRISreadTransform(MRIS *mris,const char *fname);
int MRISaddToValues(MRI_SURFACE *mris, float val) ;
int MRISsetValues(MRI_SURFACE *mris, float val) ;
LABEL *MRISannotation_to_label(MRI_SURFACE *mris, int annot_index) ;


// thickness stuff
int MRISminimizeThicknessFunctional(MRI_SURFACE *mris,
                                    INTEGRATION_PARMS *parms,
                                    float max_thick) ;
int MRISmeasureDistanceBetweenSurfaces(MRI_SURFACE *mris,
                                       MRI_SURFACE *mris2,
                                       int signed_dist) ;
int MRISwriteCoordsToIco(MRI_SURFACE *mris,
                         MRI_SURFACE *mris_ico,
                         int which_vertices);
int MRISvertexCoord2XYZ_float (VERTEX * v,
                               int which,
                               float  *x, float  *y, float  *z);
int MRISsampleFaceNormal(MRI_SURFACE *mris, int fno, double x, double y, double z, 
                         float *px, float *py, float *pz) ;
int
MRISsampleFaceCoordsCanonical(MHT *mht, MRI_SURFACE *mris, float x, float y, float z, int which, 
                              float *px, float *py, float *pz) ;
int MRISsampleFaceCoords(MRI_SURFACE *mris,
                         int fno,
                         double x, double y, double z,
                         int which_coords,
                         int which_barcyentric,
                         float *px, float *py, float *pz);
MRI *MRISlaplacian(MRI_SURFACE *mris,
                   MRI *mri_cmatrix,
                   double inner_width,
                   double outer_width);
double MRISsampleValue(MRI_SURFACE *mris,
                       FACE *f,
                       double xp, double yp, double zp, 
                       int which, MRI *mri_vals) ;
int MRIScopyAnnotationsToMarkedIndex(MRI_SURFACE *mris) ;
int MRISmaxMarked(MRI_SURFACE *mris) ;
int MRISscaleVertexCoordinates(MRI_SURFACE *mris, double scale) ;
int MRIScurvToMarked(MRI_SURFACE *mris) ;
int MRISreadMarked(MRI_SURFACE *mris, const char *sname) ;
int MRISstoreTangentPlanes(MRI_SURFACE *mris, int which_vertices) ;
double MRISsampleFace(MRI_SURFACE *mris, int fno, int which, double x, double y, double z, double val0, double val1, double val2);
int MRISrepositionSurface(MRI_SURFACE *mris, MRI *mri, int *target_vnos, float *target_vals, 
                          int nv, int nsize, double sigma, int flags)  ;
int MRISrepositionSurfaceToCoordinate(MRI_SURFACE *mris, MRI *mri, int target_vno, 
                                      float tx, 
                                      float ty, 
                                      float tz, 
                                      int nsize, double sigma, int flags)  ;
int face_barycentric_coords(MRI_SURFACE *mris, int fno, int which_vertices,
                            double cx, double cy, double cz, double *pl1, double *pl2, double *pl3) ;

MRI *MRIScomputeFlattenedVolume(MRI_SURFACE *mris,
                                MRI *mri,
                                double res,
                                int nsamples,
                                int normalize,
                                MRI **pmri_vertices,
                                int smooth_iters,
                                double wm_dist,
                                double outside_dist);
int MRIStrinarizeCurvature(MRI_SURFACE *mris, float binarize_thresh) ;
int MRISthresholdValIntoMarked(MRI_SURFACE *mris, float thresh) ;
int MRISremoveCompressedRegions(MRI_SURFACE *mris, double min_dist) ;
int  MRISweightedSoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs) ;
int MRIStaubinSmooth(MRI_SURFACE *mris, int niters, double lambda, double mu, int which) ;
MRI_SURFACE *MRISconcat(MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris) ;

#define TAUBIN_UNIFORM_WEIGHTS   0
#define TAUBIN_INVERSE_WEIGHTS   1
#define TAUBIN_EDGE_WEIGHTS      2

MRI *MRIcomputeLaminarVolumeFractions(MRI_SURFACE *mris, double res, MRI *mri_src, MRI *mri_vfracs) ;


#endif // MRISURF_H
