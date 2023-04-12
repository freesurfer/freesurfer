#pragma once
/**
 * @brief MRI_SURFACE utilities.
 *
 * Utilities, constants and structure definitions for manipulation
 * and i/o of surfaces derived from MRI volumes.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

/* To find the VERTEX and FACE classes, see mrisurf_FACE_VERTEX_MRIS_generated.h*/

#include "mrisurf_aaa.h"

#include "minc_volume_io.h"
#include "label.h"
#include "mrishash.h"
#include "json.h"
using json = nlohmann::json;

struct OverlayInfoStruct;

#define CONTRAST_T1    0
#define CONTRAST_T2    1
#define CONTRAST_FLAIR 2

typedef const MRIS MRIS_const;
    // Ideally the MRIS and all the things it points to would be unchangeable via this object but C can't express this concept easily.



void MRISctr(MRIS *mris, int max_vertices, int max_faces, int nvertices, int nfaces);
void MRISdtr(MRIS *mris);
    //
    // These are only way to create a surface and destroy a surface.
    // The destruction frees any child structures.
    //
    // There are functions below for editing the surface by adding vertex positions, edges, face information, etc.
    // There is even one for creating one similar to a subset of another's vertices and faces - MRIScreateWithSimilarTopologyAsSubset

MRI_SURFACE* MRISoverAlloc              (                   int max_vertices, int max_faces, int nvertices, int nfaces) ;
MRI_SURFACE* MRISalloc                  (                                                    int nvertices, int nfaces) ;
    //
    // Allocates an MRIS then calls MRISctr

void MRISfree(MRIS **pmris) ;
    //
    // The only way to delete a surface.  All the substructures are also freed, and the *pmris set to nullptr
    
void MRISreallocVerticesAndFaces(MRIS *mris, int nvertices, int nfaces) ;
    //
    // Used by code that is deforming the surface

void MRISacquireNverticesFrozen(MRIS *mris);
void MRISreleaseNverticesFrozen(MRIS *mris);
    //
    // Used by some code that depends on this number not changing
    // but lots of such dependent code have not been changed to use this
    
// Make and free temp properties for all the vertices
// whilest making sure they are the right size by stopping the nvertices changing when any such exist
//
float* MRISmakeFloatPerVertex(MRIS *mris);
void   MRISfreeFloatPerVertex(MRIS *mris, float** pp);

// The following create a copy of a surface, whilest deleteing some vertices and some faces
// The faces that are kept must not reference any vertices which are not kept.
// There is a function to help generate such a mapping...
//
MRIS* MRIScreateWithSimilarTopologyAsSubset(
    MRIS_const * src,
    size_t       nvertices,     // the mapToDstVno entries must each be less than this
    int const*   mapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t       nfaces,        // the mapToDstFno entries must each be less than this
    int const*   mapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to
    //
    // Used to extract some of the vertices, some of the faces, and to renumber them.
    // mrisurf_deform uses this for several purposes.

MRIS* MRIScreateWithSimilarXYZAsSubset(
    MRIS_const * src,
    size_t       nvertices,     // the mapToDstVno entries must each be less than this
    int const*   mapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t       nfaces,        // the mapToDstFno entries must each be less than this
    int const*   mapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to

MRIS* MRIScreateWithSimilarPropertiesAsSubset(
    MRIS_const * src,
    size_t       nvertices,     // the mapToDstVno entries must each be less than this
    int const*   mapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t       nfaces,        // the mapToDstFno entries must each be less than this
    int const*   mapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to

void MRIScreateSimilarTopologyMapsForNonripped(
    MRIS_const * src,
    size_t     * pnvertices,     // the mapToDstVno entries must each be less than this
    int const* * pmapToDstVno,   // src->nvertices entries, with the entries being -1 (vertex should be ignored) or the vno within the dst surface the face maps to
    size_t     * pnfaces,        // the mapToDstFno entries must each be less than this
    int const* * pmapToDstFno);  // src->nfaces entries, with the entries being -1 (face should be ignored) or the fno within the dst surface the face maps to
    // The caller should free the * pmapToDstVno and * pmapToDstFno

// There are various fields in the VERTEX and FACE and others that are used for many purposes at different
// times, and clashing uses could cause a big problem.  Start working towards a reservation system.
//
typedef enum MRIS_TempAssigned {
    MRIS_TempAssigned_Vertex_marked,
    MRIS_TempAssigned_Vertex_marked2,
    MRIS_TempAssigned__end
} MRIS_TempAssigned;

int  MRIS_acquireTemp      (MRIS* mris, MRIS_TempAssigned temp);                               // save the result to use later to ...
void MRIS_checkAcquiredTemp(MRIS* mris, MRIS_TempAssigned temp, int MRIS_acquireTemp_result);  // ... check that you own it
void MRIS_releaseTemp      (MRIS* mris, MRIS_TempAssigned temp, int MRIS_acquireTemp_result);  // ... be allowed to release it


// Support for writing traces that can be compared across test runs to help find where differences got introduced  
//
void mrisVertexHash(MRIS_HASH* hash, MRIS const * mris, int vno);

void mris_hash_init (MRIS_HASH* hash, MRIS const * mris);
void mris_hash_add  (MRIS_HASH* hash, MRIS const * mris);
void mris_hash_print(MRIS_HASH const* hash, FILE* file);
void mris_print_hash(FILE* file, MRIS const * mris, const char* prefix, const char* suffix);
void mris_print_diff(FILE* file, MRIS const * lhs, MRIS const * rhs);

// This structs are used with the TESS functions
typedef struct tface_type_
{
  int imnr,i,j,f;
  int num;
  int v[4];
}
tface_type;
typedef struct tvertex_type_
{
  int imnr,i,j;
  int num;
  int f[9];
}
tvertex_type;

typedef struct {
  // Structure for computing info at an equal number of hops from
  // a given vertex (not exactly the same as equal mm distance)
  int cvtx; // center vertex
  int nhops;
  char *hit; // vector to indicate whether vertex has been hit
  int *nperhop; // number of vertices per hop
  int **vtxlist; // list of vertices for each hop
  int *nperhop_alloced; // number of vertices alloced per hop
} SURFHOPLIST;


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
  float        sigma ;                    /* blurring scale */
  float        radius ;
  float        scale ;
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


#include "dtk.fs.h"
class SurfacePointSet {
public:
  MRIS *surf;
  json *pPointSet;
  int m_nhops=2;
  int m_fill_holes = 1;
  int m_prune_by_angle = 1;
  double AngleDegThresh = 60; // prune points where vector relative to normal is more than this (abs)
  int m_prune_by_dist = 1;
  double DistMmThresh = 1; // prune points when vertex is within this dist of the target
  int m_debug = 0;
  FILE *m_debug_fp = stdout;
  MRI *mri = NULL; // template from surf
  std::vector<SURFHOPLIST *> m_shl;
  std::vector<int> m_shl_vtxlist;
  std::vector<int> psvtxlist; // vertices that map from the point set
  std::vector<int> vtxlist; // all vertices
  std::vector<int> npervtxlist; // number of target points that map to each vertex
  std::vector<std::vector<double>> txyzlist; // target coords for all vertices
  std::vector<std::vector<double>> dxyz; // delta from vertex to target
  std::vector<double> angle; // angle between surface normal and vector to target
  std::vector<double> dist; // dist from vertex to its target
  int MapPointSet(void);
  fsPointSet ConvertToPointSet(void); // for debugging, per-vertex, tkreg
  fsPointSet VerticesToPointSet(void); // for debugging, per-vertex, tkreg
  int Print(FILE *fp);
  MRI *MakeMask(void);
  double CostAndGrad(double weight, int ComputeGradient);
  DTK_TRACK_SET *ConvertToTrack(int nsteps);
  int WriteAsPatch(const char *fname,int ndil);
  int PruneByAngle(void);
  int PruneByDistance(void);
  ~SurfacePointSet();
};

class INTEGRATION_PARMS {
 public:
  double  tol ;               /* tolerance for terminating a step */
  float   l_angle ;           /* coefficient of angle term */
  float   l_pangle ;          /* coefficient of "positive" angle term - only penalize narrower angles */
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
  float   l_norm ;            /* coefficient of surface self-intersection */
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
  float   l_spring_nzr;       /* coefficient of spring term with non-zero resting length*/
  float   l_spring_nzr_len;   /* resting length of spring term with non-zero resting length*/
  float   l_hinge;            /* coefficient of hinge angle term */
  float   l_repulse ;         /* repulsize force on tessellation */
  float   l_repulse_ratio ;   /* repulsize force on tessellation */
  float   l_boundary ;        /* coefficient of boundary term */
  float   l_dist ;            /* coefficient of distance term */
  float   l_location ;        // target location term
  float   l_targetpointset ;  // target pointset location term
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
  
  FILE    * const fp ;        /* for logging results, write by calling INTEGRATION_PARMS_<various> functions */
  
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
  int          nonmax ;  // apply nonmax suppression in reg
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
  double       target_intensity ;
  double       stressthresh ;
  int          explode_flag ;
  SurfacePointSet  *TargetPointSet;
  
  /*
    Introduce all initializers in an effort to avoid some memset() calls
    which gcc8 finds unsettling.

    I do find myself contemplating whether this class may benefit from
    some splitting up - whether by composition, inheritance, or both.
   */
  INTEGRATION_PARMS(FILE* file = nullptr)
    : tol(0), l_angle(0), l_pangle(0), l_area(0), l_parea(0),
      l_nlarea(0), l_nldist(0), l_thick_normal(0), l_thick_spring(0),
      l_ashburner_triangle(0), l_ashburner_lambda(0), l_corr(0),
      l_ocorr(0), l_pcorr(0), l_curv(0), l_norm(0), l_scurv(0), l_lap(0),
      l_link(0), l_spring(0), l_nlspring(0), l_max_spring(0),
      l_spring_norm(0), l_tspring(0), l_nltspring(0), l_nspring(0),
      l_spring_nzr(0), l_spring_nzr_len(0), l_hinge(0), l_repulse(0),
      l_repulse_ratio(0), l_boundary(0), l_dist(0), l_location(0), l_targetpointset(0),
      l_neg(0), l_intensity(0), l_sphere(0), l_expand(0), l_grad(0),
      l_convex(0), l_tsmooth(0), l_surf_repulse(0), l_osurf_repulse(0),
      l_external(0), l_thick_parallel(0), l_thick_min(0), l_shrinkwrap(0),
      l_expandwrap(0), l_unfold(0), l_dura(0), l_histo(0), l_map(0),
      l_map2d(0), dura_thresh(0), mri_dura(nullptr), n_averages(0),
      min_averages(0), first_pass_averages(0), nbhd_size(0), max_nbrs(0),
      write_iterations(0),
      base_name(), /* Should default initialize array to zero */
      projection(0), niterations(0), a(0), b(0), c(0), start_t(0), t(0),
      fp(file), // Highlighted as the sole configurable value
      Hdesired(0), integration_type(0), momentum(0), dt(0), base_dt(0),
      flags(0), dt_increase(0), dt_decrease(0), error_ratio(0), epsilon(0),
      desired_rms_height(0), starting_sse(0), ending_sse(0), scale(0),
      mrisp(nullptr), frame_no(0), mrisp_template(nullptr),
      mrisp_blurred_template(nullptr), area_coef_scale(0), sigma(0),
      nfields(0),
      fields(), /* Array should initialize to zero */
      mri_brain(nullptr), mri_smooth(nullptr), user_parms(nullptr),
      mri_dist(nullptr), target_radius(0), ignore_energy(0), check_tol(0),
      overlay_dir(nullptr), nsurfaces(0), mri_ll(nullptr), rmin(0), rmax(0),
      var_smoothness(0), vsmoothness(nullptr), dist_error(nullptr),
      area_error(nullptr), geometry_error(nullptr), which_norm(0),
      abs_norm(0), grad_dir(0), fill_interior(0), rms(0), complete_dist_mat(0),
      nsubjects(0), nlabels(0), mht_array(nullptr), mris_array(nullptr),
      mris_ico(nullptr), mht(nullptr), smooth_averages(0), ico_order(0),
      remove_neg(0), mri_hires(nullptr), mri_hires_smooth(nullptr),
      mri_vno(nullptr), mri_template(nullptr), which_surface(0),
      trinarize_thresh(0), nonmax(0), smooth_intersections(0), uncompress(0),
      min_dist(0), h_wm(nullptr), h_gm(nullptr), h_nonbrain(nullptr),
      mri_labels(nullptr), mri_white(nullptr), mri_aseg(nullptr), hwm(nullptr),
      hgm(nullptr), hout(nullptr), h2d_wm(nullptr), h2d_gm(nullptr), 
      h2d_out(nullptr), h2d(nullptr), mri_volume_fractions(nullptr), 
      mri_dtrans(nullptr), resolution(0), target_intensity(0), stressthresh(0),
      explode_flag(0) {}
  
};

void INTEGRATION_PARMS_copy   (INTEGRATION_PARMS* dst, INTEGRATION_PARMS const * src);

void INTEGRATION_PARMS_setFp  (INTEGRATION_PARMS* parms, FILE* file);
void INTEGRATION_PARMS_openFp (INTEGRATION_PARMS* parms, const char* name, const char* mode);
void INTEGRATION_PARMS_closeFp(INTEGRATION_PARMS* parms);
void INTEGRATION_PARMS_copyFp (INTEGRATION_PARMS* dst, INTEGRATION_PARMS const * src);


extern double (*gMRISexternalGradient)                   (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern double (*gMRISexternalSSE)                        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern double (*gMRISexternalRMS)                        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern int    (*gMRISexternalTimestep)                   (MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
extern int    (*gMRISexternalRipVertices)                (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
extern int    (*gMRISexternalClearSSEStatus)             (MRI_SURFACE *mris) ;
extern int    (*gMRISexternalReduceSSEIncreasedGradients)(MRI_SURFACE *mris, double pct) ;

// These are for backwards compatibility for when we don't want to fix
// the vertex area. Default is to fix, but this can be changed by setting to 0.
//
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
                          float *dmin, int which_vertices);
double MRISfindMinDistanceVertexWithDotCheckXYZ(MRI_SURFACE *mris, double xs, double ys, double zs, 
						MRI *mri, double dot_dir, int *pvtxno);
double MRISfindMinDistanceVertexWithDotCheck(MRI_SURFACE *mris, int c, int r, int s, 
					     MRI *mri, double dot_dir, int *pvtxno);                              


double       MRIScomputeCorrelationError(MRI_SURFACE *mris, MRI_SP *mrisp_template, int fno) ;
int          MRISallocExtraGradients(MRI_SURFACE *mris) ;
MRI_SURFACE  *MRISread(const char *fname, bool dotkrRasConvert=true) ;
MRI_SURFACE  *MRISreadOverAlloc(const char *fname, double nVFMultiplier) ;
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
int          MRISreadPatch(MRI_SURFACE *mris,const  char *pname, bool dotkrRasConvert=true) ;
int          MRISreadPatchNoRemove(MRI_SURFACE *mris,const  char *pname, bool dotkrRasConvert=true) ;
int          MRISreadTriangleProperties(MRI_SURFACE *mris,
                                        const  char *mris_fname) ;
int          MRISreadBinaryCurvature(MRI_SURFACE *mris,
                                     const  char *mris_fname) ;
int          MRISreadCurvatureFile(MRI_SURFACE *mris,const char *fname, MRI *curvmri=NULL, std::vector<OverlayInfoStruct> *poverlayinfo=NULL) ;
float        *MRISreadNewCurvatureVector(MRI_SURFACE *mris,
                                         const  char *sname) ;
float        *MRISreadCurvatureVector(MRI_SURFACE *mris,const  char *sname) ;
int          MRISreadFloatFile(MRI_SURFACE *mris,const char *fname) ;
#define MRISreadCurvature MRISreadCurvatureFile

int mrisReadAsciiCurvatureFile(MRI_SURFACE *mris, const char *fname, MRI *curvmri=NULL);
int mrisWriteAsciiCurvatureFile(MRI_SURFACE *mris, char *fname);
MRI_SURFACE *MRISreadVTK(MRI_SURFACE *mris, const char *fname, MRI *curvmri=NULL);

MRI *MRISloadSurfVals(const char *srcvalfile,
                      const char *typestring,
                      MRI_SURFACE *Surf,
                      const char *subject,
                      const char *hemi,
                      const char *subjectsdir);
int          MRISreadValues(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadAnnotation(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteVertexLocations(MRI_SURFACE *mris, char *fname, int which_vertices) ;
int          MRISimportVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices);
int          MRISwriteAnnotation(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadCTABFromAnnotationIfPresent(const char *fname,
                                                 COLOR_TABLE** out_table);
int          MRISisCTABPresentInAnnotation(const char *fname, int* present);
int          MRISreadValuesBak(MRI_SURFACE *mris,const  char *fname) ;
int          MRISreadImagValues(MRI_SURFACE *mris,const  char *fname) ;
int          MRIScopyImagValuesToValues(MRI_SURFACE *mris) ;
int          MRIScopyMarksToAnnotation(MRI_SURFACE *mris) ;
int          MRIScopyValsToAnnotations(MRI_SURFACE *mris) ;
int          MRIScopyValuesToImagValues(MRI_SURFACE *mris) ;
int          MRIScopyStatsToValues(MRI_SURFACE *mris) ;
int          MRIScopyStatsFromValues(MRI_SURFACE *mris) ;


int MRISsetCroppedToZero(MRI_SURFACE *mris) ;
int MRIScopyFromCropped(MRI_SURFACE *mris, int which) ;
int MRIScopyToCropped(MRI_SURFACE *mris, int which) ;
int MRISwriteCropped(MRI_SURFACE *mris, const char *fname) ;

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
int mrisNormalFace(MRIS *mris, int fac, int n, float norm[]);
int          MRISwritePrincipalDirection(MRI_SURFACE *mris, int dir_index, const  char *fname) ;
int          MRISwriteVTK(MRI_SURFACE *mris,const  char *fname);
int          MRISwriteCurvVTK(MRI_SURFACE *mris, const char *fname);
int          MRISwriteGeo(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteICO(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteSTL(MRI_SURFACE *mris, const char *fname) ;
int          MRISwritePatchAscii(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteDists(MRI_SURFACE *mris,const  char *fname) ;
int          MRISwriteCurvature(MRI_SURFACE *mris, const  char *fname, const char *curv_name=NULL) ;
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

int   MRISintegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_avgs);
int   mrisLogIntegrationParms(FILE *fp, MRI_SURFACE *mris,
			      INTEGRATION_PARMS *parms) ;

int          MRISsampleDistances(MRI_SURFACE *mris, int *nbr_count,int n_nbrs);
int          MRISsampleAtEachDistance(MRI_SURFACE *mris, int nbhd_size,
                                      int nbrs_per_distance) ;
int          MRISscaleDistances(MRI_SURFACE *mris, float scale) ;
MRI_SURFACE  *MRISradialProjectOntoEllipsoid(MRI_SURFACE *mris_src,
    MRI_SURFACE *mris_dst,
    float a, float b, float c);
    
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

int          MRISsurfaceRASToVoxel(MRI_SURFACE *mris, MRI *mri, double r, 
                                   double a, double s, 
                                   double *px, double *py, double *pz) ;
                                   
// THE FOLLOWING IS NOT THREAD SAFE!
int          MRISsurfaceRASToVoxelCached(MRI_SURFACE *mris,
                                         MRI *mri,
                                         double r, double a, double s, 
                                         double *px, double *py, double *pz) ;


typedef struct MRIS_SurfRAS2VoxelMap {
    MRI*    mri;                            // was held in mris->mri_sras2vox
    MATRIX* sras2vox;                       // was held in mris->m_sras2vox
    VECTOR * volatile v1[_MAX_FS_THREADS];  // used to avoid repeated allocations
    VECTOR * volatile v2[_MAX_FS_THREADS];  // used to avoid repeated allocations
} MRIS_SurfRAS2VoxelMap;

void MRIS_useRAS2VoxelMap(MRIS_SurfRAS2VoxelMap * cache_nonconst,   // accesses cache thread safely
        MRI const * const mri,
        double r, double a, double s, double *px, double *py, double *pz);
    
void MRIS_loadRAS2VoxelMap(MRIS_SurfRAS2VoxelMap* cache,            // not thread safe
        MRI const * const mri, MRI_SURFACE const * const mris);
void MRIS_unloadRAS2VoxelMap(MRIS_SurfRAS2VoxelMap* cache);         // not thread safe

MRIS_SurfRAS2VoxelMap* MRIS_makeRAS2VoxelMap(                       // not thread safe
        MRI const * const mri, MRI_SURFACE const * const mris);
void MRIS_freeRAS2VoxelMap(MRIS_SurfRAS2VoxelMap** const cachePtr); // not thread safe

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

MRI          *MRISwriteIntoVolume(MRI_SURFACE *mris, MRI *mri, int which) ;
MRI_SURFACE  *MRISreadFromVolume(MRI *mri, MRI_SURFACE *mris, int which) ;



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
void         MRISclearCurvAndVal2(MRIS *mris) ;
void         MRISclearD(MRIS *mris) ;
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
int MRIScountRipped(MRIS *mris);
int MRIScountAllMarked(MRIS *mris);
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

typedef struct MRISPfunctionValResultForAlpha {
    double curr;
    double next;
} MRISPfunctionValResultForAlpha;

void MRISPfunctionVal_radiusR(                                                      // returns the value that would be stored in resultsForEachFno[0] for fnoLo
                              MRI_SURFACE_PARAMETERIZATION *mrisp,  
                              MRISPfunctionValResultForAlpha* resultsForEachAlpha,  // must be numAlphas elements
                              float r, float x, float y, float z, 
                              int fno, bool getNextAlso,                            // always fills in resultsForEachAlpha.curr for fno, optionally fills in .next for fno+1
                              const float* alphas, float numAlphas,                 // rotate x,y,z by these alphas (radians) and get the values
                              bool trace) ;                                         // note: this rotation is around the z axis, hence z does not change
                             
double       MRISPfunctionValTraceable(MRI_SURFACE_PARAMETERIZATION *mrisp,
                              float desired_radius,
                              float x, float y, float z, int fno, bool trace) ;
double       MRISPfunctionVal(MRI_SURFACE_PARAMETERIZATION *mrisp,
                              float desired_radius,
                              float x, float y, float z, int fno) ;
                              
MRI_SP       *MRIStoParameterizationBarycentric(MRI_SURFACE *mris, MRI_SP *mrisp, float scale, int fno) ;
MRI_SURFACE  *MRISfromParameterizationBarycentric(MRI_SP *mrisp, MRI_SURFACE *mris, int fno) ;

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
int MRIScoordsFromParameterizationBarycentric(MRIS *mris, MRI_SP *mrisp, int which_vertices);


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
int          MRISPwrite(MRI_SP *mrisp, const char *fname) ;

int          MRISwriteArea(MRI_SURFACE *mris,const  char *sname) ;
int          MRISwriteMarked(MRI_SURFACE *mris,const  char *sname) ;
int          MRISmarkedToCurv(MRI_SURFACE *mris) ;
int          MRISdToCurv(MRI_SURFACE *mris) ;

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
#define VERTEX_NORMALS      12
#define PIAL_NORMALS        13
#define WHITE_NORMALS       14

// Two ways of saving and restoring the VERTEX:xyz values
//
// The push/pop is preferable to using VERTEX members because it uses all the entries in a cache line
//
typedef struct MRISsavedXYZ MRISsavedXYZ;
MRISsavedXYZ* MRISsaveXYZ(MRIS *mris);
void          MRISloadXYZ(MRIS *mris, MRISsavedXYZ const *  pMRISsavedXYZ);  // does set
void          MRISfreeXYZ(MRIS *mris, MRISsavedXYZ      ** ppMRISsavedXYZ);  // does not set
void          MRISpopXYZ (MRIS *mris, MRISsavedXYZ      ** ppMRISsavedXYZ);  // set then free

int MRISsaveVertexPositions   (MRIS *mris, int which) ;
int MRISrestoreVertexPositions(MRIS *mris, int which) ;

int MRISrestoreNormals(MRIS *mris, int which) ;
int MRISsaveNormals   (MRIS *mris, int which) ;

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

#define DEFAULT_A  122.0f
#define DEFAULT_B  122.0f
#define DEFAULT_C  122.0f
#define DEFAULT_RADIUS  100.0f

#define MAX_DIM    DEFAULT_B

#define DT_INCREASE  1.1 /* 1.03*/
#define DT_DECREASE  0.5
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
int   MRISsetValBaks(MRI_SURFACE *mris, float val) ;
int   MRISaverageD(MRI_SURFACE *mris, int navgs) ;
int   MRISgaussianFilterD(MRI_SURFACE *mris, double wt) ;
int   MRISaverageVals(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageVal2baks(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageVal2s(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageMarkedVals(MRI_SURFACE *mris, int navgs) ;
int   MRISaverageMarkedStats(MRI_SURFACE *mris, int navgs) ;
int   MRISsoapBubbleVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISsoapBubbleOrigVertexPositions(MRI_SURFACE *mris, int navgs) ;
int   MRISsoapBubbleTargetVertexPositions(MRI_SURFACE *mris, int navgs) ;
MRI   *MRISwriteSurfaceIntoVolume(MRI_SURFACE *mris, MRI *mri_template,
                                  MRI *mri) ;
int   MRISmeasureCorticalThickness(MRI_SURFACE *mris, int nbhd_size,
                                   float max_thickness) ;

int  MRISmeasureThicknessFromCorrespondence(MRI_SURFACE *mris, MHT *mht, float max_thick) ;
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
int   MRISreverseCoords(MRI_SURFACE *mris, int which_direction, int reverse_face_order, int which_coords) ;
int   MRISreverse(MRI_SURFACE *mris, int which, int reverse_face_order) ;
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
MRIScomputeBorderValuesV6(MRI_SURFACE *mris,MRI *mri_brain,
			  MRI *mri_smooth, double inside_hi, double border_hi,
			  double border_low, double outside_low, double outside_hi,
			  double sigma, float max_thickness, FILE *log_fp,
			  int which, MRI *mri_mask, double thresh,
			  int flags,  MRI *mri_aseg,int junk1, int junk2);
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
                              MRI *mri_mask, double thresh, int flags, MRI *mri_aseg,
                              int vno_start, int vno_stop);
int MRIScomputePialTargetLocationsMultiModal(MRI_SURFACE *mris,
                              MRI *mri_T2,
                              LABEL **labels,
                              int nlabels,
                              int contrast_type, MRI *mri_aseg, double T2_min_inside, double T2_max_inside, 
			      double T2_min_outside, double T2_max_outside, double max_outward_dist,
			      double left_inside_peak_pct,
			      double right_inside_peak_pct,
			      double left_outside_peak_pct,
			      double right_outside_peak_pct,
			      double wm_weight,
 			      double pial_sigma,
 			      MRI *mri_T1);
int  MRIScomputeWhiteSurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                                   MRI *mri_smooth);
int  MRIScomputeGraySurfaceValues(MRI_SURFACE *mris, MRI *mri_brain,
                                  MRI *mri_smooth, float gray_surface);
int  MRIScomputeDistanceErrors(MRI_SURFACE *mris, int nbhd_size,int max_nbrs);
int  MRISsetAllMarks(MRI_SURFACE *mris, int mark) ;
int  MRISscaleCurvature(MRI_SURFACE *mris, float scale) ;
int  MRISwriteTriangularSurface(MRI_SURFACE *mris,const char *fname) ;

int  MRISbuildFileName_read(MRI_SURFACE *mris, const char *sname, char *fname) ;
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
int MRISltaMultiply(MRIS *surf, const LTA *lta);

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
int MRISprintVertexInfo(FILE *fp, MRIS *surf, int vertexno);
int MRISprintSurfQualityStats(FILE *fp, MRIS *surf);
int MRISprettyPrintSurfQualityStats(FILE *fp, MRIS *surf);

int MRISmergeIcosahedrons(MRI_SURFACE *mri_src, MRI_SURFACE *mri_dst) ;
//int MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico) ;
MRI *MRIScoverSeg(MRIS *mris, MRI *mri_bin, MRI *mri_cover_seg, int surftype);

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

int MRISmarkOrientationChanges(MRI_SURFACE *mris);
MRIS* MRISextractMainComponent(MRI_SURFACE *mris,
                               int do_not_extract,
                               int verbose,
                               int *ncpts);
MRIS* MRISextractMarkedVertices(MRIS *mris);
MRIS* MRISremoveRippedSurfaceElements(MRIS *mris);

MRI_SURFACE *MRIScorrectTopology(MRI_SURFACE *mris, MRI *mri, 
   MRI *mri_wm, int nsmooth, TOPOLOGY_PARMS *parms, const char *defectbasename);

int MRISsmoothOnSphere(MRIS* mris, int niters);
int mrisCountIntersectingFaces(MRIS *mris, int*flist , int nfaces);
int MRIScountNegativeFaces(MRI_SURFACE *mris) ;
int MRISevertSurface(MRI_SURFACE *mris) ;
int MRISripDefectiveFaces(MRI_SURFACE *mris) ;
int MRISdefects2Seg(MRIS *surf, MRI *defects, int offset, MRI *vol);
int MRISunrip(MRI_SURFACE *mris) ;
int MRISdivideLongEdges(MRI_SURFACE *mris, double thresh) ;
int MRISdivideEdges(MRI_SURFACE *mris, int npoints) ;
int MRISremoveTriangleLinks(MRI_SURFACE *mris) ;
int MRISsetOriginalFileName(const char *orig_name) ;
int MRISsetSulcFileName(const char *sulc_name) ;
int MRISsetCurvatureName(int nth, const char *name);
int MRISprintCurvatureNames(FILE *fp);
int MRISsetInflatedFileName(char *inflated_name) ;
int MRISsetRegistrationSigmas(float *sigmas, int nsigmas) ;

int MRISextractVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices) ;
int MRISimporttVertexCoords(MRI_SURFACE *mris, float *locations[3], int which_vertices) ;

float* MRISexportCurv(MRIS* mris);
int MRISextractCurvatureVector(MRI_SURFACE *mris, float *curvs) ;
int MRISextractCurvatureDoubleVector(MRI_SURFACE *mris, double *curvs) ;
int MRISextractCurvatureVectorDouble(MRI_SURFACE *mris, double *curvs, int offset) ;
#define MRISexportCurvatureVector  MRISextractCurvatureVector

int MRISimportCurvatureVector(MRI_SURFACE *mris, float *curvs) ;
int MRISimportValVector(MRI_SURFACE *mris, float *vals) ;
int MRISexportValVector(MRI_SURFACE *mris, float *vals) ;
int MRISexportValVectorDouble(MRI_SURFACE *mris, double *vals, int offset) ;
int MRISimportValFromMatrixColumn(MRI_SURFACE *mris, MATRIX *m, int col) ;
int MRISimportValFromMRI(MRI_SURFACE *mris, MRI *mri, int frame) ;
MRI *MRISreadParameterizationToSurface(MRI_SURFACE *mris, char *fname) ;

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
float  MRISdistanceToSurface(MRI_SURFACE *mris, MHT *mht,
                             float x0, float y0, float z0,
                             float nx, float ny, float nz) ;
int    MRISexpandSurface(MRI_SURFACE *mris,
                         float distance,
                         INTEGRATION_PARMS *parms, int use_thickness, int nsurfaces) ;
int MRISripZeroThicknessRegions(MRI_SURFACE *mris) ;


/* cortical ribbon */
MRI   *MRISribbon(MRI_SURFACE *inner_mris,
                  MRI_SURFACE *outer_mris,
                  MRI *mri_src,
                  MRI *mri_dst);
MRI   *MRISaccentuate(MRI *mri_src,
                      MRI *mri_dst,
                      int lo_thresh,
                      int hi_thresh);
MRI *MRISmakeBoundingVolume(MRIS *mris, double resolution);
MRI *MRISfillInterior(MRI_SURFACE *mris,
                      double resolution,
                      MRI *mri_interior) ;
int MRISfillInteriorRibbonTest(char *subject, int UseNew, FILE *fp);
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
#define MRIS_LABEL_FILE               10

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
#define VERTEX_MARKS       11


// consolidated from two identical copies in mris_flatten and mris_sphere
int MRISscaleUp(MRI_SURFACE *mris) ;

int MRISclearOrigArea(MRIS *mris) ;
void MRISclearOrigAreaAndVal2(MRIS *mris);
int MRISclearOrigDistances(MRIS *mris) ;
int MRIScombine(MRIS *mris_src, MRIS *mris_total,
                MRIS_HASH_TABLE *mht, int which) ;
int MRISsphericalCopy(MRI_SURFACE *mris_src, MRI_SURFACE *mris_total,
                      MRIS_HASH_TABLE *mht, int which) ;
int   MRISorigAreaToCurv(MRI_SURFACE *mris) ;
int   MRISareaToCurv(MRI_SURFACE *mris) ;
int   MRISnormalize(MRI_SURFACE *mris, int dof, int which) ;

int  MRIScopyMRI(MRIS *Surf, MRI *Src, int Frame, const char *Field);
MRI *MRIcopyMRIS(MRI *mri, MRIS *surf, int Frame, const char *Field);

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
int MRISaddCommandLine(MRI_SURFACE *mris, const std::string& cmdline);
int MRISgetReadFrame(void);
int MRISabsCurvature(MRI_SURFACE *mris) ;
int MRISabsVals(MRI_SURFACE *mris) ;
int MRISsmoothFrames(MRI_SURFACE *mris, MRI *mri, int navgs) ;
int MRISwriteFrameToValues(MRI_SURFACE *mris, MRI *mri, int frame) ;
int MRISreadFrameFromValues(MRI_SURFACE *mris, MRI *mri, int frame) ;
MRI *MRISar1(MRIS *surf, MRI *src, MRI *mask, MRI *ar1);
MRI *MRISfwhmFromAR1Map(MRIS *surf, MRI *mask, MRI *ar1map);
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
int MRISremoveIntersections(MRI_SURFACE *mris, int FillHoles) ;
int MRIScopyMarkedToMarked2(MRI_SURFACE *mris) ;
int MRIScopyMarked2ToMarked(MRI_SURFACE *mris) ;
int MRIScopyMarkedToMarked3(MRI_SURFACE *mris) ;
int MRIScopyMarked3ToMarked(MRI_SURFACE *mris) ;
int MRISexpandMarked(MRI_SURFACE *mris) ;
int MRISnotMarked(MRI_SURFACE *mri) ;
double MRISsmoothingArea(MRIS *mris, int vtxno, int niters);

int MRISopenMarked(MRI_SURFACE *mris, int order) ;
int MRIScloseMarked(MRI_SURFACE *mris, int order) ;
int MRISerodeMarked(MRI_SURFACE *mris, int ndil) ;
int MRISdilateMarked(MRI_SURFACE *mris, int ndil) ;
MRI *MRISdilateVertexToSum(int vno, MRIS *surf, MRI *measure, double targetSum);
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
int MRISmarkVerticesWithValOverThresh(MRI_SURFACE *mris, float thresh) ;

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
                          float *pcsf_mode,
			  float *pwhite_std,
			  float *pgray_std,
			  float *pcsf_std);
int MRISrasToVoxel(MRI_SURFACE *mris,
                   MRI *mri,
                   double xs, double ys, double zs,
                   double *pxv, double *pyv, double *pzv) ;

int MRISripMedialWall(MRI_SURFACE *mris) ;
int MRISripMarked(MRI_SURFACE *mris) ;
int MRISripUnmarked(MRI_SURFACE *mris) ;
int MRISzeroMedialWallCurvature(MRI_SURFACE *mris) ;
int MRISvertexNormalInVoxelCoords(MRI_SURFACE *mris,
                                  MRI *mri,
                                  int vno,
                                  double *pnx, double *pny, double *pnz) ;

float* MRISgetVertexArray(MRIS *mris);
float* MRISgetVertexNormalArray(MRIS *mris);
int*   MRISgetFaceArray(MRIS *mris);
float* MRISgetFaceNormalArray(MRIS *mris);
MRIS*  MRISfromVerticesAndFaces(const float *vertices, int nvertices, const int *faces, int nfaces);

#define MRISgetCoords(v,c,vx,vy,vz) \
 MRISvertexCoord2XYZ_float(v,c,vx,vy,vz)

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
	const char*		apch_left,
	float		af_right
);

void	cprints(
	const char*		apch_left,
	const char*		apch_right
);

void	cprintd(
	const char*		apch_left,
	int		a_right
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
    	const char*			apch_message
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
    	int			apFACE_O_fno,
    	int			apFACE_I_fno
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

float	FACES_angleNormal_find(
    	MRIS*			apmris,
    	int			apFACE_I_fno,
    	int			apFACE_J_fno
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
                  
#define MRISvertexCoord2XYZ_float  MRISvertexCoord2XYZ
#define MRISvertexCoord2XYZ_double MRISvertexCoord2XYZ
int MRISvertexCoord2XYZ (VERTEX const * v, int which, float  *x, float  *y, float  *z);
int MRISvertexCoord2XYZ (VERTEX const * v, int which, double *x, double *y, double *z);

int face_barycentric_coords(
  double V0[3], double V1[3], double V2[3],
  double  cx,
  double  cy,
  double  cz,
  double *pl1,
  double *pl2,
  double *pl3);

int face_barycentric_coords(MRIS const * mris, int fno, int which_vertices,
                            double cx, double cy, double cz, double *pl1, double *pl2, double *pl3) ;

                               
int MRISsampleFaceNormal(MRI_SURFACE *mris, int fno, double x, double y, double z, 
                         float *px, float *py, float *pz) ;
int
MRISsampleFaceCoordsCanonical(MHT *mht, MRI_SURFACE *mris, float x, float y, float z, int which, 
                              float *px, float *py, float *pz) ;
MRI *MRISmapToSurface(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst, MRI *mri_src_features, MRI *mri_dst_features) ;
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
int MRIScurvToD(MRI_SURFACE *mris) ;
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

int MRISvertexNormalToVoxelScaled(MRI_SURFACE *mris,
				  VERTEX *v,
				  MRI *mri,
				  double *pnx, double *pny, double *pnz) ;
int
MRISvertexNormalToVoxel(MRI_SURFACE *mris,
			VERTEX *v,
			MRI *mri,
			double *pnx, double *pny, double *pnz) ;
MRI *MRIcomputeLaminarVolumeFractions(MRI_SURFACE *mris, double res, MRI *mri_src, MRI *mri_vfracs) ;
MRIS *MRIStessellate(MRI *mri,  int value, int all_flag);
void TESSaddFace(MRI *mri, int imnr, int i, int j, int f, int prev_flag, int *pface_index, 
		 tface_type *face, int *face_index_table0, int *face_index_table1);
void TESScheckFace(MRI *mri, int im0, int i0, int j0, int im1, int i1,int j1,
		int f, int n, int v_ind, int prev_flag, int all_flag, int value,
		tface_type *face, int *pface_index, int *face_index_table0, 
		   int *face_index_table1,	tvertex_type *vertex);
int TESSaddVertex(MRI *mri, int imnr, int i, int j, int *pvertex_index,  int *vertex_index_table, tvertex_type *vertex);
int TESSfacep(MRI *mri, int im0, int i0, int j0, int im1, int i1, int j1, int value, int all_flag);

#define SURFACE_SMOOTH_STEPS_TO_SIGMA(iter)   (sqrt((double)iter) * M_PI / 2.0)
#define SIGMA_TO_SURFACE_SMOOTH_STEPS(sigma)  SQR(2.0*sigma/M_PI)

SURFHOPLIST *SetSurfHopListAlloc(MRI_SURFACE *Surf, int nHops);
SURFHOPLIST *SetSurfHopList(int CenterVtx, MRI_SURFACE *Surf, int nHops);
int SurfHopListFree(SURFHOPLIST **shl0);
MRI *MRISarN(MRIS *surf, MRI *src, MRI *mask, MRI *arN, int N);
MRI *MRISsmoothKernel(MRIS *surf, MRI *src, MRI *mask, MRI *mrikern, MATRIX *globkern, SURFHOPLIST ***pshl, MRI *out);
int MRISmeasureLaplaceStreamlines(MRI_SURFACE *mris, MRI *mri_laplace, MRI *mri_intensity, MRI *mri_profiles) ;
MRI *MRISsolveLaplaceEquation(MRI_SURFACE *mris, MRI *mri, double res) ;
int MRIScountEdges(MRIS *surf);
int MRISedges(MRIS *surf);
int MRIScorners(MRIS *surf);
MRIS *MRIScopyMetadata(MRIS const * source, MRIS * target);
int MRISfixAverageSurf7(MRIS *surf7);
double mrisRmsValError(MRI_SURFACE *mris, MRI *mri);
void MRIScalculateCenterCOG2(MRIS *mris, double *xCOG, double *yCOG, double *zCOG);

// for sorting vertices using qsort
typedef struct{
  int vtxno;
  float x,y,z;
} VERTEX_SORT; 
int CompareVertexCoords(const void *v1, const void *v2);
// for sorting faces using qsort
typedef struct{
  int faceno;
  int v0, v1, v2;
} FACE_SORT; 
int CompareFaceVertices(const void *vf1, const void *vf2);
// for making the surface deterministic after decimation
MRIS *MRISsortVertices(MRIS *mris0);


// Create a surface
//
MRIS *MRISclone(MRIS const * mris_src) ;
MRIS* MRISunion(MRIS const * mris, MRIS const * mris2);

// Various basic export and import
//
char* MRISexportVertexRipflags(MRIS* mris) ;
void  MRISimportVertexRipflags(MRIS* mris, const char*) ;

char* MRISexportFaceRipflags(MRIS* mris) ;
void  MRISimportFaceRipflags(MRIS* mris, const char*) ;

int  MRISrestoreRipFlags(MRIS *mris) ;
int  MRISstoreRipFlags  (MRIS *mris) ;

char* MRISexportVertexRipflags(MRIS* mris) ;

// mrisurf_topology needed by more
//
//  Vertices, like Faces, come into existence when the surface is created with a vertex and face count.
//  Edges are implicit (MRI_EDGE is more than just an edge), and are created by telling each of the end vertices that they are neighbors.
//  Faces get associated with three edges associated with three vertices (VERTICES_PER_FACE is 3)
//
bool mrisCheckVertexVertexTopologyWkr(const char* file, int line, MRIS const * mris, bool always);
bool mrisCheckVertexFaceTopologyWkr  (const char* file, int line, MRIS const * mris, bool always);
inline static bool returnTrue() { return true; };
#define mrisCheckVertexFaceTopology(_MRIS)   returnTrue() // mrisCheckVertexFaceTopologyWkr  (__FILE__,__LINE__,_MRIS,false)
#define mrisCheckVertexVertexTopology(_MRIS) returnTrue() // mrisCheckVertexVertexTopologyWkr(__FILE__,__LINE__,_MRIS,false)

//  Vertices
//
int mrisVertexVSize                 (MRIS const * mris, int vno);
static int  mrisVertexNeighborIndex (MRIS const * mris, int vno1, int vno2);
static bool mrisVerticesAreNeighbors(MRIS const * mris, int vno1, int vno2);

void mrisAddEdge   (MRIS* mris, int vno1, int vno2);
void mrisRemoveEdge(MRIS *mris, int vno1, int vno2);


// Neighbourhoods
//
void mrisVertexReplacingNeighbors(MRIS * mris, int vno, int vnum);
void mrisForgetNeighborhoods     (MRIS * mris);

int  MRISresetNeighborhoodSize      (MRIS *mris, int nsize) ;
void MRISsetNeighborhoodSize        (MRIS *mris, int nsize) ;   // Doesn't compute distances if not already present
void MRISsetNeighborhoodSizeAndDist (MRIS *mris, int nsize) ;   // Always computes distances

int  MRISfindNeighborsAtVertex      (MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops);
    // sets vlist[*] to the neighboring vno
    // sets hops [*] to -1 for [vno] and the number of hops for all entries returned in the vlist
    // returns the number of neighbors


// Inputs to the metric properties
//
//      It is a good idea to call MRISfreeDistsButNotOrig(MRIS *mris);
//      before calling these, to invalidate the now-wrong distances
//
//      However there are places where the old consequences of xyz are used to compute the new values of xyz
//      so this tactic isn't always usable - especially when the previous normals are being used to move the xyz
//
void MRISsetXYZwkr(MRIS *mris, int vno, float x, float y, float z, const char * file, int line, bool* laterTime);
#define MRISsetXYZ(_MRIS,_VNO, _X,_Y,_Z) { \
    static bool _laterTime; \
    MRISsetXYZwkr((_MRIS),(_VNO),(_X),(_Y),(_Z), __FILE__, __LINE__, &_laterTime); \
  }
    //
    // VERTEX:xyz can be set directly, one at a time, or via one of the following operations
    // that iterate across many vertices.
    //
    // This results in all the derived properties (distances, face areas, normals, angles, ...) being invalid
    // until recomputed.  However the use of invalid properties is not yet detected.

void MRIScopyXYZ(MRIS *mris, MRIS* mris_from);

void MRISmemalignNFloats(size_t n, float** ppx, float** ppy, float** ppz);
void MRISexportXYZ(MRIS *mris,       float*       * ppx,       float*       * ppy,       float*       * ppz);
void MRISimportXYZ(MRIS *mris, const float* const    px, const float* const    py, const float* const   ppz);
    //
    // By importing and exporting into three arrays, it is possible to use h/w vector instructions more effectively.
    // 
    // the three vectors are malloced and filled in with the xyz values
    //       the vectors are cache-block aligned so loops processing them can be very efficient
    //
    // the xyz values are set from the vectors, the vectors are NOT freed by this call
    // the mris [xyz] lo,hi,ctr are set during import


// Deforming the MRIS
//
void mrisFindMiddleOfGray(MRIS *mris);


void MRISscaleThenTranslate (MRIS *mris, double sx, double sy, double sz, double dx, double dy, double dz);   // new = old * s + d


int  MRIStranslate (MRIS *mris, float dx, float dy, float dz);
void MRISmoveOrigin(MRIS *mris, float x0, float y0, float z0);
int  MRISscale     (MRIS *mris, double scale);

void MRISblendXYZandTXYZ(MRIS* mris, float xyzScale, float txyzScale);          // x = x*xyzScale + tx*txyzScale  etc.
void MRISblendXYZandNXYZ(MRIS* mris,                 float nxyzScale);          // x = x          + nx*nxyzScale  etc.


void MRIStranslate_along_vertex_dxdydz(MRIS*    dst, MRIS*    src, double dt);  // dst.x = src.x + src.dx*dt, etc.
void MRIStranslate_along_vertex_dxdydz(MRIS_MP* dst, MRIS_MP* src, double dt);

MRIS* MRISrotate(MRIS* mris_src, MRIS* mris_dst, float alpha, float beta, float gamma) ;

void mrisDisturbVertices(MRIS *mris, double amount);

MRIS* MRIScenter(MRIS *mris_src, MRIS *mris_dst) ;
void  MRIScenterSphere(MRIS *mris);
void  MRISrecenter(MRIS *mris, int which_move, int which_target) ;

MRIS* MRISprojectOntoTranslatedSphere(MRIS* mris_src, MRIS* mris_dst, 
    double r,
    double x0, double y0, double z0);


// xyz's immediate consequences
//
int mrisComputeSurfaceDimensions(MRIS *mris);
    // xyz lo/hi  and the cxyz

// dist and dist_orig
//      can be freed at any time
// dist is created by calls to MRISsetNeighborhoodSizeAndDist
// dist_orig must be explicitly created
//
void MRISfreeDistsButNotOrig(MRIS*    mris);
void MRISfreeDistsButNotOrig(MRISPV*  mris);
void MRISfreeDistsButNotOrig(MRIS_MP* mris);

void MRISmakeDistOrig (MRIS *mris, int vno);                        // makes it the same size as the current VERTEX.dist
void MRISgrowDistOrig (MRIS *mris, int vno, int minimumCapacity);   // same size as current or bigger
void MRISfreeDistOrigs(MRIS *mris);

int mrisComputeVertexDistances(MRIS *mris);

//  Faces
//
void mrisSetVertexFaceIndex(MRIS *mris, int vno, int fno);
    // is being used outside mrissurf_topology but shouldn't be
    
void setFaceAttachmentDeferred(MRIS* mris, bool to);                                // for performance reasons, defer adding them; or do all the deferred ones
int mrisCountAttachedFaces(MRIS* mris, int vno0, int vno1);
bool mrisCanAttachFaceToVertices(MRIS* mris, int vno1, int vno2, int vno3);         // returns whether such a face would be legal
void mrisAttachFaceToEdges   (MRIS* mris, int fno, int vno1, int vno2, int vno3);   // edges must already exist
void mrisAttachFaceToVertices(MRIS* mris, int fno, int vno1, int vno2, int vno3);   // adds any needed edges

void  MRISflipFaceAroundV1(MRIS *mris, int fno);
void  MRISreverseFaceOrder(MRIS *mris);


// Adding missing parts of the topology representation
//
void mrisCompleteTopology(MRIS *mris);
    //
    // The best way to build a tesselation is to simply create all the vertices then ...
    //      setFaceAttachmentDeferred(true)
    //      lots of calls to mrisAttachFaceToVertices()  
    //      setFaceAttachmentDeferred(false)
    //
    // However lots of existing code just adds a misc set of edges and faces by mechanisms
    // then they call this function to add any missing edges and to calculate the mris->avg_nbrs


// Vertices and Faces can have their ripflag set
//
void MRISsetRipInFacesWithRippedVertices(MRIS* mris);


// Removing 
//
void MRISrenumberRemovingRippedFacesAndVertices(MRIS *mris);

void MRISremoveRipped(MRIS *mris);
    // cleans up the neighbours etc. but does not renumber


// Marked
//
bool mrisAnyVertexOfFaceMarked(MRIS *mris, int fno);


// Vals
//
void MRISsetOriginalXYZfromXYZ(MRIS *mris);
    //
    // This includes copying the MRIS::status to the MRIS::origxyz_status
    
void MRISsetOriginalXYZwkr(MRIS *mris, int vno, float x, float y, float z, const char* file, int line, bool* laterTime);
#define MRISsetOriginalXYZ(_MRIS,_VNO,_X,_Y,_Z) \
    { static bool laterTime; MRISsetOriginalXYZwkr((_MRIS),(_VNO),(_X),(_Y),(_Z),__FILE__,__LINE__, &laterTime); }
    //
    // The values being set need to match the MRIS::origxyz_status

int mrisComputeOriginalVertexDistances(MRIS *mris);
void mrisComputeOriginalVertexDistancesIfNecessaryWkr(MRIS *mris, bool* laterTime, const char* file, int line);
#define mrisComputeOriginalVertexDistancesIfNecessary(_MRIS) \
  { static bool laterTime;  \
    mrisComputeOriginalVertexDistancesIfNecessaryWkr((_MRIS), &laterTime, __FILE__, __LINE__); \
  }

void MRIScheckForNans(MRIS *mris);
void MRIScheckIsPolyhedron(MRIS *mris, const char* file, int line);

int StuffVertexCoords(MRIS *surf, int vertexno, double p[3]);
int StuffFaceCoords(MRIS *surf, int faceno, int cornerno, double p[3]);
double MinDistToTriangleBF(double p1[3], double p2[3], double p3[3], double ptest[3], double pmin[3], double dL);
int MRISdistanceBetweenSurfacesExact(MRIS *surf1, MRIS *surf2);
int MRISnorm2Pointset(MRIS *mris, int vno, double dstart, double dend, double dstep, FILE *fp);
MRI *MRISextractNormalMask(MRIS *surf, int vno, double dstart, double dend, double dstep, double UpsampleFactor);
MRI *MRISsampleProfile(MRIS *mris, MRI *mri, double dstart, double dend, double dstep, double sigma, int interptype, MRI *profile);
int MatlabPlotFace(FILE *fp, MRIS *surf, int faceno, char color, double NormLen);
int MatlabPlotVertex(FILE *fp, MRIS *surf, int vno, char color, double NormLen);
int MatlabPlotVertexNbhd(FILE *fp, MRIS *surf, int cvno, int nhops, char color, double NormLen);
  
int MRISshrinkFaceCorner(MRIS *surf, int faceno, int nthv, double dist, double *vx, double *vy, double *vz);
double MRISshrinkFace(MRIS *surf, int faceno, double newareafraction);
int MRISshrinkFaces(MRIS *surf, double zthresh, int nmax);


#if defined(COMPILING_MRISURF_TOPOLOGY) || defined(COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED)
#define CONST_EXCEPT_MRISURF_TOPOLOGY 
#else
#define CONST_EXCEPT_MRISURF_TOPOLOGY const
#endif

#if defined(COMPILING_MRISURF_METRIC_PROPERTIES) || defined(COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND)
#define CONST_EXCEPT_MRISURF_METRIC_PROPERTIES 
#else
#define CONST_EXCEPT_MRISURF_METRIC_PROPERTIES const
#endif
    //
    // Used to find and control where various fields are written


struct face_topology_type_ {    // not used much yet
  vertices_per_face_t v;        // and being overtaken by events...
};

//  The types needed to store and describe the various properties of the mesh
//  that are independent of the representation of the mesh are found in
//  mrisurf_aaa.h
//
//  We have several ways to code against our triangular mesh, so that we can
//  support several representations and restricted access to the properties.
//
//  Abstractly, each of these representations defines a Surface as a set of 
//  Faces and Vertices, each with a set of properties that can be read and
//  perhaps written.
//
//  The traditional way is to have the code fully understand the representation
//  so the code says
//      MRIS*   mris;
//      VERTEX* v = &mris->vertices[vno];
//      if (v->ripflag) ...
//  but this style has three problems
//      1) The data is not stored in a memory-traffic efficient manner
//      2) The code has to know the data structures
//      3) The code has access to all the data members at all times,
//              rather than having read and write access only to those
//              fields that it should currently have access to
//
//  This is mrisurf_FACE_VERTEX_MRIS*.h
//
//  Most of the code uses this representation directly

#include "mrisurf_FACE_VERTEX_MRIS_generated.h"

//  To provide restricted access to these structures, using code that
//  can be redirected to other representations, we have a set of namespaces
//  which each define Surface, Vertex, and Face classes, and which provide
//  member functions to get and set the various properties of those
//  entities regardless of how they are represented.  These classes are
//  effectively pointers rather than structs - Vertex replaces VERTEX*.
//
//  After the includes in the .cpp file, specify the namespace that describes
//  what it is allowed to access
//
//      using namespace SurfaceFromMRIS::Topology;
//  
//  The Surface Face Vertex classes replace the traditional pointers,
//  so you write
//
//      Vertex v2 = v1.v(3);
//      if (v2.ripflag()) continue;
//      v2.set_marked(true);
//
//  instead of the traditional
//
//      VERTEX* v2 = &mris->vertices[v1->v[3]];
//      if (v2->ripflag) continue;
//      v2->marked = true;
//
#include "mrisurf_SurfaceFromMRIS_generated.h"


// An alternative representation has the various fields of the VERTEX and FACE above
// spread out into one array per field.
//
// This representation exploits the two facts that (a) we almost never want to
// deal with all the fields of a VERTEX as a whole, and (b) that we tend to use
// the vertices and faces with nearby indexs are about the same time, so that
// the cachelines brought into and written out of the cache will now be fully
// used.
//
// Almost none of the code uses this representation directly,
// but obviously the constructors and destructors must.
//
//  class MRISPV is the MRIS with its Properties stored in Vectors
//
#include "mrisurf_MRIS_PropertiesInVectors.h"

// Similar Surface Face Vertex classes and namespaces are generated to access
// this representation, so that template classes and functions can be created
// which can be instantiated to get high speed access to this representation
// also, yet having common source.
//
#include "mrisurf_SurfaceFromMRISPV_generated.h"


// Static function implementations
//

short        modVnum  (MRIS const *mris, int vno, short add, bool clearFirst = false);
static short setVnum  (MRIS const *mris, int vno, short to)     { return modVnum(mris,vno, to,true );     }
static short clearVnum(MRIS const *mris, int vno)               { return modVnum(mris,vno,  0,true );     }
static short vnumAdd  (MRIS const *mris, int vno, short add=+1) { return modVnum(mris,vno,add,false)-add; }
static short addVnum  (MRIS const *mris, int vno, short add=-1) { return modVnum(mris,vno,add,false);     }

static short VERTEXvnum(VERTEX_TOPOLOGY const * v, int i) {
  switch (i) {
  case 1: return v->vnum;
  case 2: return v->v2num;
  case 3: return v->v3num;
  default: cheapAssert(false); return 0;
  }    
}

static int mrisVertexNeighborIndex(MRIS const *mris, int vno1, int vno2) {
  cheapAssert(0 <= vno1 && vno1 < mris->nvertices);
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno1];
  if (!vt->v) return -1;    // happens in invalid surfaces 
  int n;
  for (n = 0; n < vt->vnum; n++) {
    if (vt->v[n] == vno2) return n;
  }
  return -1;
}


static bool mrisVerticesAreNeighbors(MRIS const * const mris, int const vno1, int const vno2)
{
  return 0 <= mrisVertexNeighborIndex(mris, vno1, vno2);
}

int MRISripMidline(MRI_SURFACE *mris, MRI *mri_aseg, MRI *mri_brain, const char *hemi, int which, int fix_mtl);
int MRIcomputeLabelNormal(MRI *mri_aseg, int x0, int y0, int z0,int label, int whalf, double *pnx, double *pny,
			  double *pnz, int use_abs);
int MRIScopyCoords(MRIS *surf, MRIS *surfcoords);
int MRISfindExpansionRegions(MRI_SURFACE *mris);
int MRISwriteField(MRIS *surf, const char **fields, int nfields, const char *outname);
MRI *MRISflatMap2MRI(MRIS *flatmap, MRI *overlay, double res, int DoAverage, MRI *out);

/**
  class AutoDetGWStats. This class houses functions used to compute
  intensity limits used in MRIScomputeBorderValues() when placing both
  the white and pial surfaces on a T1w image. This is functionality
  that used to be in mris_make_surfaces.cpp. It has mostly been copied
  over, which is why it is not very well organized. 
 */
class AutoDetGWStats
{
public:
  MRIS *mrisAD, *mrisADlh, *mrisADrh; // surface used to autodetect stats
  MRI *mri_T1, *mri_wm;
  const char *wm_name = "wm" ;
  const char *orig_name = "orig";
  //In mris_make_surfaces, "brain" is the default, but brain.finalsurfs is always used in recon-all
  const char *T1_name = "brain.finalsurfs"; 
  int hemicode = 0; //1=left, 2=right
  int use_mode = 1;
  float variablesigma = 3.0;
  double std_scale = 1.0;
  float adWHITE_MATTER_MEAN = 110;
  float MAX_WHITE = 120;
  float MAX_BORDER_WHITE = 105;
  float MIN_BORDER_WHITE = 85;
  float MIN_GRAY_AT_WHITE_BORDER = 70;
  float MAX_GRAY = 95;
  float MID_GRAY;
  float MIN_GRAY_AT_CSF_BORDER = 40;
  float MAX_GRAY_AT_CSF_BORDER = 75;
  float MIN_CSF = 10;
  float adMAX_CSF = 40;
  float white_mean, white_std, gray_mean, gray_std ;
  float white_mode, gray_mode ;
  float lh_white_mode, lh_gray_mode, rh_white_mode, rh_gray_mode ;
  float max_border_white = MAX_BORDER_WHITE;
  float min_border_white = MIN_BORDER_WHITE;
  float min_gray_at_white_border = MIN_GRAY_AT_WHITE_BORDER;
  float max_gray = MAX_GRAY;
  float max_gray_at_csf_border = MAX_GRAY_AT_CSF_BORDER;
  float min_gray_at_csf_border = MIN_GRAY_AT_CSF_BORDER;
  float min_csf = MIN_CSF;
  float max_csf = adMAX_CSF ;
  double max_gray_scale = 0.0 ;  
  double max_scale_down = .2;
  double white_inside_hi;
  double white_border_hi;
  double white_border_low;
  double white_border_low_factor=1;
  double white_outside_low;
  double white_outside_hi;
  double pial_inside_hi;
  double pial_border_hi;
  double pial_border_low;
  double pial_outside_low;
  double pial_outside_hi;
  // These indicate whether there was a manual override.
  int  max_border_white_set = 0, min_border_white_set = 0, min_gray_at_white_border_set = 0,
    max_gray_set = 0,max_gray_at_csf_border_set = 0, min_gray_at_csf_border_set = 0,
    min_csf_set = 0, max_csf_set = 0 ;
  int AutoDetectStats(const char *subject, const char *hemistr);
  int AutoDetectStats(void);
  int Write(char *fname);
  int Print(FILE *fp);
  int Read(char *fname); // from file name
  int ReadStream(FILE *fp); // read from stream
};

/*!
  \fn class ClosestVertex 
  \brief Manages finding closest vertex
 */
class ClosestVertex 
{
public:
  MRIS *surf=NULL;
  LTA *lta=NULL; // can be in either direction, and it will figure it out
  MHT *hash=NULL;
  int hashres=16;
  int which_vertices = CURRENT_VERTICES;
  int m_debug=0;
  MATRIX *scanner2tkreg=NULL;
  MATRIX *tkreg2tkreg=NULL;
  int InitMatrices(void);
  int PrintMatrices(FILE *fp);
  int ClosestTkReg(double x, double y, double z, double *dist);
  int ClosestScanner(double x, double y, double z, double *dist);
  int Closest(double x, double y, double z, double *dist, MATRIX *M);
  int InitHash(void){ 
    if(m_debug) printf("Using hash %d\n",hashres);
    hash = MHTcreateVertexTable_Resolution(surf, which_vertices, hashres);
    return(0);
  }
  int ReadLTA(char *ltafile){ // can pass it "nofile" to ignore, good for cmdargs
    if(strcmp(ltafile,"nofile")!=0) {
      lta = LTAread(ltafile);
      if(lta == NULL) return(1);
      else            return(0);
    }
    return(0);
  }
  int WriteVertexDist(char *outfile, int vno, double dist)
  { // can pass it "nofile" to ignore, good for cmdargs
    double vx=0,vy=0,vz=0;
    if(vno >=0 && vno < surf->nvertices){
      // should take into account which_vertices
      vx = surf->vertices[vno].x;
      vy = surf->vertices[vno].y;
      vz = surf->vertices[vno].z;
    }
    printf("%8d  %8.4f  %8.4f %8.4f %8.4f\n",vno,dist,vx,vy,vz);fflush(stdout);
    if(outfile != NULL && strcmp(outfile,"nofile")!=0) {
      FILE *fp = stdout;
      fp = fopen(outfile,"w");
      if(fp==NULL){
	printf("ERROR: could not open %s for writing\n",outfile);
	return(1);
      }
      fprintf(fp,"%8d  %8.4f  %8.4f %8.4f %8.4f\n",vno,dist,vx,vy,vz);
      fclose(fp);
    }
    return(0);
  }
};
