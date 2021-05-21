/*
 *
 */
/*
 * surfaces Author: Bruce Fischl
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrisurf_base.h"

int UnitizeNormalFace = 1;
int RmsValErrorRecord = 0;


#if (!SPHERE_INTERSECTION)
static int mrisComputeCanonicalEdgeBasis(
    MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2, double origin[3], double e0[3], double e1[3]);
#endif

#if AVERAGE_AREAS
static int mrisAverageAreas(MRI_SURFACE *mris, int num_avgs, int which);
#endif

/*--------------------- CONSTANTS AND MACROS -------------------------*/

/* 16777215 = 0xFFFFFF */

/* 77557 is horizontal with positive area */
/* 126906 has dy = 0 with negative area */
/* 102961 has a = 0 */
/* 77115 is > pi/2, face 1 */
/* v 82875, f 0 has negative area and horizontal */
/* v 115365, f 1 has negative area and is vertical */
/* v 75530, f 4 has negative area and is vertical */
#if 0
#define DEBUG_FACE(vno, fno) (((fno) == 2) && (Gdiag & DIAG_SURFACE) && (vno == 79881))
#endif
#define DEBUG_FACE(vno, fno) (((fno) == 4) && (Gdiag & DIAG_SURFACE) && (DEBUG_VERTEX(vno)))
#define VDEBUG_FACE(fno) (DEBUG_FACE(fno) && 0)
#define DEBUG_VERTEX(v) (((v) == 75530) && (Gdiag & DIAG_SURFACE) && 1)
#define VDEBUG_VERTEX(v) (((v) == 77115) && (Gdiag & DIAG_SURFACE) && 0)

/*--------------------------------------------------------------------*/
/*------------------------ STATIC DATA -------------------------------*/

/*-------------------------- FUNCTIONS -------------------------------*/
double (*gMRISexternalGradient)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL;
double (*gMRISexternalSSE)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL;
double (*gMRISexternalRMS)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL;
int (*gMRISexternalRipVertices)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL;
int (*gMRISexternalClearSSEStatus)(MRI_SURFACE *mris) = NULL;
int (*gMRISexternalReduceSSEIncreasedGradients)(MRI_SURFACE *mris, double pct) = NULL;


/*-----------------------------------------------------
  ------------------------------------------------------*/
#ifdef BEVIN_REPRODUCIBLES_CHECK
static void reproducible_check(double cell, double val, int line, int* count) 
{
    (*count)++;
    if (cell == val) return;
    fprintf(stderr, "reproducible_check %s:%d diff %g:%g %g using %d threads, count:%d\n",
        __FILE__, line, cell, val, cell-val, omp_get_max_threads(), *count);
    exit(1);
}
#endif


// ------------------- Structures -------------------------//
// now in topology/topo_parms.h


// -------------------- Declaration of Macros --------------------- //


// ----------------- Declaration of Static Functions ---------------------- //
int MRISSfree(SMALL_SURFACE **pmriss)
{
  SMALL_SURFACE *mriss;

  mriss = *pmriss;
  *pmriss = NULL;
  free(mriss->vertices);
  free(mriss);
  return (NO_ERROR);
}

int MRISaddCommandLine(MRI_SURFACE *mris, const std::string& cmdline)
{
  if (mris->ncmds >= MAX_CMDS)
    fs::error() << "can't add cmd to surface since max cmds (" << mris->ncmds <<  ") has been reached";

  int i = mris->ncmds++;
  mris->cmdlines[i] = (char *)calloc(cmdline.size() + 1, sizeof(char));
  strcpy(mris->cmdlines[i], cmdline.c_str());
  return NO_ERROR;
}

// Support for writing traces that can be compared across test runs to help find where differences got introduced  
//
static size_t showHashCalc;

static bool vertix_n_hash_add(size_t vectorSize, MRIS_HASH* hashVector, MRIS const ** mrisPVector, FILE* showDiff, int vno)
{
    unsigned int i;
    #define SEP
    #define ELTP(TARGET, MBR) // don't hash pointers.   Sometime may implement hashing their target
    #define ELTX(TYPE,   MBR) // don't hash excluded elements
#ifdef SEPARATE_VERTEX_TOPOLOGY
    #define ELTT(TYPE,   MBR) \
        for (i = 0; i < vectorSize; i++) {                                                              \
            MRIS_HASH  * hash = &hashVector[i];                                                         \
            MRIS const * mris = mrisPVector[i];                                                         \
            VERTEX_TOPOLOGY const * vt = &mris->vertices_topology[vno];                                 \
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(&vt->MBR), sizeof(vt->MBR));        \
            if (showHashCalc) {                                                                         \
                fprintf(stdout, "After %s hash is %ld\n", #MBR, hash->hash);                            \
            }                                                                                           \
            if (showDiff && i > 0 && hash->hash != hashVector[0].hash) {                                \
                fprintf(showDiff, "Differ at vertices_topology:%d field %s\n", vno, #MBR);              \
                return false;                                                                           \
            }                                                                                           \
        }                                                                                               \
        // end of macro
    LIST_OF_VERTEX_TOPOLOGY_ELTS
    #undef ELTT
#endif
    #define ELTT(TYPE,   MBR) \
        for (i = 0; i < vectorSize; i++) {                                                              \
            MRIS_HASH  * hash = &hashVector[i];                                                         \
            MRIS const * mris = mrisPVector[i];                                                         \
            VERTEX const * v = &mris->vertices[vno];                                                    \
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(&v->MBR), sizeof(v->MBR));          \
            if (showHashCalc) {                                                                         \
                fprintf(stdout, "After %s hash is %ld\n", #MBR, hash->hash);                            \
            }                                                                                           \
            if (showDiff && i > 0 && hash->hash != hashVector[0].hash) {                                \
                fprintf(showDiff, "Differ at vertices:%d field %s\n", vno, #MBR);                       \
                return false;                                                                           \
            }                                                                                           \
        }                                                                                               \
        // end of macro
    #undef ELTT
    #undef ELTX
    #undef ELTP
    #undef SEP
    
    // Now include some of the pointer targets
    //
    for (i = 0; i < vectorSize; i++) {
        MRIS_HASH  * hash = &hashVector[i];
        MRIS const * mris = mrisPVector[i];
        VERTEX_TOPOLOGY const * vt = &mris->vertices_topology[vno];
        VERTEX          const * v  = &mris->vertices         [vno];
        int vsize = mrisVertexVSize(mris, vno);
        if (vt->v) {
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(vt->v),        vsize * sizeof(vt->v[0]));
            if (showHashCalc) {
                fprintf(stdout, "After v hash is %ld\n", hash->hash);
            }
        }
        if (vt->f) {
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(vt->f),        vt->num * sizeof(vt->f[0]));
            if (showHashCalc) {
                fprintf(stdout, "After f hash is %ld\n", hash->hash);
            }
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(vt->n),        vt->num * sizeof(vt->n[0]));
            if (showHashCalc) {
                fprintf(stdout, "After n hash is %ld\n", hash->hash);
            }
        }
        if (v->dist) {
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(v->dist),      vsize * sizeof(v->dist[0]));
            if (showHashCalc) {
                fprintf(stdout, "After dist hash is %ld\n", hash->hash);
            }
        }
        if (v->dist_orig) {
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(v->dist_orig), vsize * sizeof(v->dist_orig[0]));
            if (showHashCalc) {
                fprintf(stdout, "After dist_orig hash is %ld\n", hash->hash);
            }
        }
        if (showDiff && i > 0 && hash->hash != hashVector[0].hash) {
            fprintf(showDiff, "Differ at vertices:%d field v f n dist and dist_orig\n", vno);
            return false;
        }
    }

    return true;
}

static bool face_n_hash_add(size_t vectorSize, MRIS_HASH* hashVector, MRIS const ** mrisPVector, FILE* showDiff, int fno)
{
    unsigned int i;
    #define SEP
    #define ELTP(TARGET,NAME) // don't hash pointers
    #define ELTX(TARGET,NAME) // don't hash these fiellds
    #define ELTT(TYPE,       MBR) \
        for (i = 0; i < vectorSize; i++) {                                                              \
            MRIS_HASH  * hash = &hashVector[i];                                                         \
            MRIS const * mris = mrisPVector[i];                                                         \
            FACE* face = &mris->faces[fno];                                                             \
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(&face->MBR), sizeof(face->MBR));    \
            if (showDiff && i > 0 && hash->hash != hashVector[0].hash) {                                \
                fprintf(showDiff, "Differ at face:%d field %s\n", fno, #MBR);                           \
                return false;                                                                           \
            }                                                                                           \
        }                                                                                               \
        // end of macro
    LIST_OF_FACE_ELTS
    #undef ELTP
    #undef ELTX
    #undef ELTT
    #undef SEP
    return true;
}

static bool mris_n_hash_add(size_t vectorSize, MRIS_HASH* hashVector, MRIS const ** mrisPVector, FILE* showDiff)
{
    size_t i;
    #define SEP
    #define ELTP(TARGET, MBR) // don't hash pointers.   Sometime may implement hashing their target
    #define ELTX(TYPE,   MBR) 
    #define ELTT(TYPE,   MBR)                                                                           \
        for (i = 0; i < vectorSize; i++) {                                                              \
            MRIS_HASH  * hash = &hashVector[i];                                                         \
            MRIS const * mris = mrisPVector[i];                                                         \
            hash->hash = fnv_add(hash->hash, (const unsigned char*)(&mris->MBR), sizeof(mris->MBR));    \
            if (showHashCalc) {                                                                         \
                fprintf(stdout, "After mris.%s hash is %ld\n", #MBR, hash->hash);                       \
            }                                                                                           \
            if (showDiff && i > 0 && hash->hash != hashVector[0].hash) {                                \
                fprintf(showDiff, "Differ at field %s\n", #MBR);                                        \
                return false;                                                                           \
            }                                                                                           \
        }                                                                                               \
        // end of macro
    LIST_OF_MRIS_ELTS
    #undef ELTT
    #undef ELTX
    #undef ELTP
    #undef SEP

    // Now include some of the pointer targets
    //
    int vno;
    for (vno = 0; vno < mrisPVector[0]->nvertices; vno++) {
        if (!vertix_n_hash_add(vectorSize, hashVector, mrisPVector, showDiff, vno)) return false;
    }
    
    int fno;
    for (fno = 0; fno < mrisPVector[0]->nfaces; fno++) {
        if (!face_n_hash_add(vectorSize, hashVector, mrisPVector, showDiff, fno)) return false;
    }
    
    return true;
}

void mrisVertexHash(MRIS_HASH* hash, MRIS const * mris, int vno) {
    hash->hash = fnv_init();
    vertix_n_hash_add(1, hash, &mris, nullptr, vno);
}

void mris_hash_add(MRIS_HASH* hash, MRIS const * mris)
{
    mris_n_hash_add(1, hash, &mris, NULL);
}

void mris_hash_init (MRIS_HASH* hash, MRIS const * mris)
{
    hash->hash = fnv_init();
    if (mris) mris_hash_add(hash, mris);
}

void mris_hash_print(MRIS_HASH const* hash, FILE* file)
{
    fprintf(file, "%ld", hash->hash);
}

void mris_print_hash(FILE* file, MRIS const * mris, const char* prefix, const char* suffix) {
    MRIS_HASH hash;
    double const * pd = &mris->avg_vertex_dist;
    void*  const * pp = (void**)pd;
    fprintf(stdout, "mris.nsize:%d mris.avg_vertex_dist:%f %p\n", mris->nsize, *pd, *pp);
    
    static size_t 
        showHashCount = 0, 
        showHashLimit = 0;  // 0 means never shows details

    bool showHash = (++showHashCount == showHashLimit);
    
    if (showHash) { showHashCalc++; showHashLimit *= 2; fprintf(stdout, "showHashCount:%ld\n", showHashCount); }
    mris_hash_init(&hash, mris);
    if (showHash) --showHashCalc;

    fprintf(file, "%sMRIS_HASH{",prefix);
    mris_hash_print(&hash, file);
    fprintf(file, "}%s",suffix);
}


void mris_print_diff(FILE* file, MRIS const * lhs, MRIS const * rhs) {
    MRIS_HASH hashPair[2]; hashPair[0].hash = hashPair[1].hash = fnv_init();
    MRIS const * mrisPair[2]; mrisPair[0] = lhs; mrisPair[1] = rhs; 
    mris_n_hash_add(2, hashPair, mrisPair, file);
}

/*!
  \fn int MRISripMidline()

  \brief Was "fix_midline() in mris_make_surfaces.cpp. Note: "fix"
  here means to "hold in place", not to "repair something that is
  broken".  It also does more than midline.

  Sets the v->ripped flag at certain vertices (which then holds them
  in place later on). Does not unrip vertices; if a vertex is ripped,
  then it is not procsssed. All unripped v->marked are set to 0 at the
  beginning; all v->marked (ripped or otherwise) are set to 0 at end
  of the function. All unripped v->marked2 are set to 0 at the
  begining of the function. The ripped vertices will have
  v->val=intensity at the vertex and v->d=0.

  It finds the midline mostly by finding subcortical GM, ventricular
  CSF, or CC in the aseg.presurf within a few mm of the vertex. These
  vertices are marked (v->marked=1).  It will also mark some vertices
  that are near putamen or a lesion. All ripped vertices are marked.

  The marks are spatially filtered by eliminating clusters that have
  fewer than 5 vertices or, if annotation is present, that are less
  than 60% "unknown". This last step unmarks the putamen marks done
  before (bug).

  If there is an annotation, it will only make a difference in two
  places: (1) for entorhinal (where vertices are never marked/ripped),
  and (2) in some areas of "unknown" cortex (see above). #1 probably
  does not have much of an impact because because there is usually
  some WM between entorhinal and hippo. #2 probably removes a few
  vertices (and causes the putamen bug). 

  v->marked2 - affects MRIScortexLabel(), forces labeling as cortex. Also
  influences this program the second time around in that vertices with
  v->marked2=1 are not reprocessed. 

  v->val2=1 for Lesions when fitting GRAY_CSF
  #FIX
*/
int MRISripMidline(MRI_SURFACE *mris, MRI *mri_aseg, MRI *mri_brain, const char *hemi, int which, int fix_mtl)
{
  int vno, label, contra_wm_label, nvox=0, total_vox=0, adjacent=0;
  int wm_label, gm_label, nlabels, n, index, annotation, entorhinal_index ;
  VERTEX   *v ;
  double   xv, yv, zv, val, xs, ys, zs, d, nx, ny, nz ;
  LABEL    **labels ;
  int nmarked, nmarked2, nripped;

  printf("Entering: MRISripMidline(): inhibiting deformation at non-cortical midline structures...\n") ;
  printf("  which=%d, fix_mtl=%d, using annot = %d\n",which, fix_mtl, mris->ct != NULL);
  printf("#FML0# MRISripMidline(): nripped=%d\n",MRIScountRipped(mris));
  fflush(stdout);

  if (stricmp(hemi, "lh") == 0){
    contra_wm_label = Right_Cerebral_White_Matter ;
    wm_label = Left_Cerebral_White_Matter ;
    gm_label = Left_Cerebral_Cortex ;
  }
  else{
    contra_wm_label = Left_Cerebral_White_Matter ;
    wm_label = Right_Cerebral_White_Matter ;
    gm_label = Right_Cerebral_Cortex ;
  }

  // Clear the deck
  MRISclearMarks(mris) ;  // Sets v->marked=0  for all unripped vertices
  MRISclearMark2s(mris) ; // Sets v->marked2=0 for all unripped vertices

  if (mris->ct)
    CTABfindName(mris->ct, "entorhinal", &entorhinal_index);
  else
    entorhinal_index = -1 ;

  // Loop over vertices ==================================
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

    if(v->ripflag || v->marked2 > 0){
      // It was ripped previously - it should still be excluded, but
      // still mark it.  Note: all unripped marked2 are set to zero
      // above, so if(v->marked2>0) is meaningless.
      v->marked = 1 ;
      continue ;
    }

    if (mris->ct)
      CTABfindAnnotation(mris->ct, v->annotation, &index);

    if (vno == Gdiag_no ) {
      printf("vno %d: annotation %x, index %d, EC %d\n", vno, v->annotation, index, entorhinal_index) ;
      DiagBreak() ;
    }

    // don't freeze vertices that are in EC and very close to hippocampus
    // How do you know from this that the vertex is close to hippo?
    if (mris->ct && index == entorhinal_index)  
      continue ;

    // search outwards. Hidden parameters 2mm and 0.5mm
    for (d = 0 ; d <= 2 ; d += 0.5) {
      xs = v->x + d*v->nx ;
      ys = v->y + d*v->ny ;
      zs = v->z + d*v->nz ;

      // Sample the aseg at this distance
      MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      label = nint(val) ;
      if (vno == Gdiag_no){
	printf("vno %d: dist %2.2f - %s (%d)\n", vno, d, cma_label_to_name(label), label) ;
	DiagBreak() ;
      }

      if (d > 0 && label == gm_label)
	break ;  // cortical gray matter of this hemisphere - doesn't matter what is outside it

      // If the vertex lands on a certain label, then v->marked=1 and
      // v->d=0 and v->val=val at this location. The val may be
      // changed at later distances but marked and d will not even if
      // the point lands in an acceptable label. marked2 may be set to
      // 1 if the label is a lesion or WMSA.  Labels excluded from
      // thist list include ipsilateral cortex and WM.
      if(label == contra_wm_label ||
          label == Left_Lateral_Ventricle ||
          label == Left_vessel ||
          label == Right_vessel ||
          label == Optic_Chiasm ||
          label == Left_choroid_plexus ||
          label == Right_choroid_plexus ||
          label == Third_Ventricle ||
          label == Right_Lateral_Ventricle ||
          ((label == Left_Accumbens_area ||
            label == Right_Accumbens_area) &&  // only for gray/white
           which == GRAY_WHITE) ||
          ((label == Left_Lesion ||
            label == Right_Lesion ||
            label == WM_hypointensities ||
            label == Left_WM_hypointensities ||
            label == Right_non_WM_hypointensities ||
            label == Left_non_WM_hypointensities ||
            label == Right_WM_hypointensities) &&  // only for gray/white
           which == GRAY_WHITE) ||
          label == Left_Caudate ||
          label == Right_Caudate ||
          label == Left_Pallidum ||
          label == Right_Pallidum ||
          IS_CC(label) ||
          ((IS_HIPPO(label)  || IS_AMYGDALA(label)) && fix_mtl)  ||
          label == Right_Thalamus_Proper ||
          label == Left_Thalamus_Proper ||
          label == Brain_Stem ||
          label == Left_VentralDC ||
          label == Right_VentralDC)
      {
	// below are labels where the intensities aren't useful, so just freeze surface there
	if (label == Left_Lesion || label == Right_Lesion || IS_WMSA(label))
	  v->marked2 = 1 ; // afects the cortex.label
        if (label == Left_Putamen || label == Right_Putamen)
          DiagBreak() ;
        if (vno == Gdiag_no)
          DiagBreak() ;
	// Sample brain intensity at vertex (not at distance from vertex)
        MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ; // not redundant
        MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
        v->val = val ;
        v->d = 0 ;
        v->marked = 1 ;
      }
    } // end loop over distance

    if (vno == Gdiag_no)
      DiagBreak() ;

    // Compute the normal to the edge of the label if this is putamen
    MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
    MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
    label = nint(val) ;
    if (IS_PUTAMEN(label)) {
      // 3=whalf in voxels, hidden parameter
      // 1=use_abs, hidden parameter
      // nx, ny, nz are in vox coords making orientation a hidden parameter
      // Are nx,ny,nz ever used?
      if(vno == Gdiag_no) printf("vno=%d, putamen label: running compute_label_normal\n",vno);
      MRIcomputeLabelNormal(mri_aseg, xv, yv, zv, label, 3, &nx, &ny, &nz, 1) ;
    }
    else
    {
      nx = ny = nz = 0 ;
    }

    /*
      for gray/white surface, if we are in insula, don't want to let the
      surface diverge into the putamen.
    */
    if (which == GRAY_WHITE) {

      // search inwards. Hidden parameters 2mm and 0.5mm
      for (d = 0 ; d <= 2 ; d += 0.5) {

	// Sample the aseg at this distance
        xs = v->x - d*v->nx ;
        ys = v->y - d*v->ny ;
        zs = v->z - d*v->nz ;
        MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
        MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
        label = nint(val) ;
	if (vno == Gdiag_no){
	  printf("vno %d: dist %2.2f - %s (%d)\n", vno, -d, cma_label_to_name(label), label) ;
	  DiagBreak() ;
	}

	// these are labels where the intensities aren't useful, so just freeze surface there
	// hidden parameter d<1mm. v->marked2 affects the cortical label
	if (d < 1 && (label == Left_Lesion || label == Right_Lesion || IS_WMSA(label)))
	  v->marked = v->marked2 = 1 ;

        if (IS_PUTAMEN(label) || IS_ACCUMBENS(label) || IS_CLAUSTRUM(label)){
	  // 3=whalf in voxels, hidden parameter. Voxel resolution is
	  //   a hidden parameter since whalf is in voxels
	  // 1=use_abs, hidden parameter, but abs used regardless below
	  // nx, ny, nz are in vox coords making orientation a hidden parameter
          MRIcomputeLabelNormal(mri_aseg, xv, yv, zv, label, 3,&nx, &ny, &nz, 1) ;
          if(fabs(nx) > fabs(ny) && fabs(nx) > fabs(nz))  {
	    // edge is mostly oriented in column direction
            if(vno == Gdiag_no)
              DiagBreak() ;
	    if(IS_ACCUMBENS(label)){
	      // Same as put and claust but val not set
	      v->d = 0 ;
	      v->marked = 1 ;
	      if (Gdiag & DIAG_SHOW && vno == Gdiag_no)
		printf("marking vertex %d as adjacent to accumbens on midline\n", vno);
	    }
	    else {
	      // Sample brain intensity at the vertex (not the distance)
	      MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ; // not redundant
	      MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
	      v->val = val ;
	      v->d = 0 ;
	      v->marked = 1 ;
	      if (Gdiag & DIAG_SHOW && vno == Gdiag_no)
		printf("marking vertex %d as adjacent to putamen/claustrum in insula\n", vno);
	    }
	    break ;
          }
        } // if putamen, accumbens, or claustrum

      } // loop over distance
    }// if GRAY_WHITE


    // search inwards 2mm step 0.5mm (hidden parameters)
    for (d = 0 ; d <= 2 ; d += 0.5) {

      // Sample the aseg at this distance
      xs = v->x - d*v->nx ;
      ys = v->y - d*v->ny ;
      zs = v->z - d*v->nz ;
      MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      label = nint(val) ;
 
      // Found ispilateral GM or WM next to surface (1.1mm hidden
      // parameter), so don't go any deeper.  This is probably why the
      // surface near entorhinal is not affected (or not much) even
      // when no annotation is passed.
      if(d < 1.1 && (label == wm_label || label == gm_label))
        break ;  

      if ((which == GRAY_CSF) && (d < 1) && (label == Left_Lesion || label == Right_Lesion || IS_WMSA(label)))
	v->val2 = 1 ;

      if ((label == contra_wm_label ||
           label == Left_vessel ||
           label == Right_vessel ||
           label == Optic_Chiasm ||
           label == Left_choroid_plexus ||
           label == Right_choroid_plexus ||
           label == Left_Lateral_Ventricle ||
           label == Third_Ventricle ||
           label == Right_Lateral_Ventricle ||
           ((label == Left_Accumbens_area ||
             label == Right_Accumbens_area) &&
            which == GRAY_WHITE)||
           label == Left_Caudate ||
           label == Right_Caudate ||
           label == Left_Pallidum ||
           IS_CC(label) ||
           label == Right_Thalamus_Proper ||
           label == Left_Thalamus_Proper ||
           label == Right_Pallidum ||
           label == Brain_Stem ||
           label == Left_VentralDC ||
           label == Right_VentralDC) ||
          // putamen can be adjacent to insula in aseg for pial
          (which == GRAY_WHITE && (d < 1.1) &&
           (label == Left_Putamen || label == Right_Putamen)))

      {
        if (label == Left_Putamen || label == Right_Putamen)
          DiagBreak() ;

        if((label == Left_Lateral_Ventricle || label == Right_Lateral_Ventricle) && d > 1)  
        {
	  // In calcarine ventricle can be pretty close to wm surface
	  // but could affect any vertex that is near L/R Lat Vent.
	  // d>1mm hidden parameter
          break ;
        }
        if (vno == Gdiag_no)
          DiagBreak() ;

        MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
        MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
        v->val = val ;
        v->d = 0 ;
        v->marked = 1 ;
      }

    } // end loop over inward distance

    /* Now check for putamen superior to this vertex. If there's a lot
       of it there, then we are in basal forebrain and not cortex. */
    if (which == GRAY_WHITE) {
      // Project in the row direction (not the normal direction) 10mm step 0.5 mm (hidden par)
      // For a conformed volume, row direction = superior/inferior
      adjacent = total_vox = nvox = 0;
      for (d = 0 ; d <= 10 ; d += 0.5, total_vox++){
	// Sample superiorly. This leads to a dependence on how the
	// head is oriented with respect to the voxel coordinates
	// (hidden parameter) since the surface coordinates are
	// aligned with the voxel coords.  It also leads to a
	// dependence on voxel axis orientation since the row
	// direction is not always superior/inf.
        xs = v->x ;
        ys = v->y ;
        zs = v->z + d ;  // z is along the row direction
        MRISsurfaceRASToVoxelCached(mris, mri_aseg, xs, ys, zs, &xv, &yv, &zv);
        MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
        label = nint(val) ;
        if (label == Left_Putamen || label == Right_Putamen) {
          nvox++ ;
          if (d < 1.5) // d<1.5mm hidden parameter
            adjacent = 1 ;  // right next to putamen
        }
      } // loop over distance

      // if more than 50% of the samples are putamen and vertex is within 1.5mm of putamen
      // and edge of putamen is mostly in the row direction
      if (adjacent && (double)nvox/(double)total_vox > 0.5)
      {
        MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
        MRIsampleVolumeType(mri_aseg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
        label = nint(val) ;
	// 3=whalf in voxels, hidden parameter. Voxel resolution is
	//   a hidden parameter since whalf is in voxels
	// 1=use_abs, hidden parameter, but abs used regardless below
	// nx, ny, nz are in vox coords making orientation a hidden parameter
        MRIcomputeLabelNormal(mri_aseg, xv, yv, zv, label, 3, &nx, &ny, &nz, 1) ;
        if (ny > 0 && fabs(ny) > fabs(nx) &&  fabs(ny) > fabs(nz))  {
	  // Normal is mostly in the row direction. Note that ny will always be >=0 since use_abs=1
          if (vno == Gdiag_no)
            DiagBreak() ;
          MRISvertexToVoxel(mris, v, mri_aseg, &xv, &yv, &zv) ;
          MRIsampleVolume(mri_brain, xv, yv, zv, &val) ;
          v->val = val ;
          v->d = 0 ;
          v->marked = 1 ;
        }
      }

    } // end if GRAY_WHITE

  } // end loop over vertices =================================================

  if (Gdiag_no >= 0) {
    int index ;
    v = &mris->vertices[Gdiag_no] ;
    if (mris->ct)
      CTABfindAnnotation(mris->ct, v->annotation, &index);
    else
      index = -1 ;
    printf("v %d: ripflag = %d before connected components, annot %d (%d)\n",
           Gdiag_no, v->marked, v->annotation, index) ;
    if (v->marked == 0)
      DiagBreak() ;
    else
      DiagBreak() ;
  }

  // Dilate and erode the marked vertices
  MRISdilateMarked(mris, 3) ;
  MRISerodeMarked(mris, 3) ;

  // Get a list of marked clusters (ie, connected components) whose
  // area is greater than 1mm2 (hidden parameter)
  MRISsegmentMarked(mris, &labels, &nlabels, 1) ;
  
  if (Gdiag_no > 0)
    printf("v %d: ripflag = %d after morphology\n", Gdiag_no, mris->vertices[Gdiag_no].marked) ;

  // Go through all the clusters
  for (n = 0 ; n < nlabels ; n++)  {

    // If there are fewer than 5 points in the cluster, then discard it
    if(labels[n]->n_points < 5) {
      int i, threadno = 0 ;
      #ifdef HAVE_OPENMP 
      threadno = omp_get_thread_num();
      #endif
      printf("removing %d vertices from ripped group in thread:%d\n",labels[n]->n_points,threadno); 
      for (i = 0 ; i < labels[n]->n_points ; i++) 
        mris->vertices[labels[n]->lv[i].vno].marked = 0 ;
    }

    // If there is an annoation that includes "unknown", then ...
    if (mris->ct && CTABfindName(mris->ct, "unknown", &index) == NO_ERROR) {
      double pct_unknown;
      int    i ;
      VERTEX *v ;
      CTABannotationAtIndex(mris->ct, index, &annotation) ;

      // Compute the fraction of this cluster that is in the unknown annotation
      for (pct_unknown = 0.0, i = 0 ; i < labels[n]->n_points ; i++) {
	v = &mris->vertices[labels[n]->lv[i].vno] ;
        if ((v->annotation == annotation || v->annotation == 0) &&
	    (v->marked2 == 0))  // will be 1 if a lesion or WMSA (keep frozen if so)
        {
          pct_unknown = pct_unknown + 1 ;
        }
      }
      pct_unknown /= (double)labels[n]->n_points ;

      // If this fraction is < 0.6 (hidden parameter), unmark all vertices in the cluster
      // Won't this undo the marks near putamen or lesion? Answer: yes, this is a bug
      // that makes the surface worse when the annotation is loaded.
      if (pct_unknown < .6) {
        printf("deleting segment %d with %d points - only %2.2f%% unknown, v=%d\n",
               n,labels[n]->n_points,100*pct_unknown, labels[n]->lv[0].vno) ;
        for (i = 0 ; i < labels[n]->n_points ; i++) {
	  v = &mris->vertices[labels[n]->lv[i].vno] ;
	  if (v->marked2 == 0){
	    mris->vertices[labels[n]->lv[i].vno].marked = 0 ;
	    if (labels[n]->lv[i].vno  == Gdiag_no)
	      printf("removing mark from v %d due to non-unknown aparc\n",Gdiag_no) ;
	  }
        }
      }
    }

    LabelFree(&labels[n]) ;

  } // end loop over clusters
  free(labels) ;

  if (Gdiag_no > 0)
    printf("v %d: ripflag = %d after connected components\n",
           Gdiag_no, mris->vertices[Gdiag_no].marked) ;

  // This is where the effects of this routine are instantiated
  // Sets v->ripflag=1 if v->marked==1  Note: does not unrip any vertices.
  MRISripMarked(mris) ; 

  // Count up the final number of marked vertices
  nmarked=0, nmarked2=0, nripped=0;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    if(mris->vertices[vno].marked)  nmarked++;
    if(mris->vertices[vno].marked2) nmarked2++;
    if(mris->vertices[vno].ripflag) nripped++;
  }
  printf("#FML# MRISripMidline(): nmarked=%d, nmarked2=%d, nripped=%d\n",nmarked,nmarked2,nripped);
  fflush(stdout);

  // Sets all marked vertices to 0, even the ripped ones
  MRISsetAllMarks(mris, 0) ;

  return(NO_ERROR) ;
}

/*!
  \fn MRIcomputeLabelNormal()
  \brief This function searches whalf voxels around the given voxel
  for voxels that are on the edge of the given label. It then computes
  the average normal of the edge over these voxels. If use_abs=1, then
  the average of the absolute value of the normals is computed. The
  normal is in voxel units, and if the volume orientation changes, then the
  interpretation of the normal changes making orientation effectively
  a hidden parameter.
*/
int MRIcomputeLabelNormal(MRI *mri_aseg, int x0, int y0, int z0,
                     int label, int whalf, double *pnx, double *pny,
                     double *pnz, int use_abs)
{
  int xi, yi, zi, xk, yk, zk, nvox = 0, val, dx, dy, dz, xn, yn, zn ;
  double  nx, ny, nz, mag ;

  nx = ny = nz = 0.0 ;

  // Search an area of 3x3x3 around xyz0 (whalf usually=3)
  for (xk = -whalf ; xk <= whalf ; xk++)  {
    xi = mri_aseg->xi[x0+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)    {
      yi = mri_aseg->yi[y0+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)      {
        zi = mri_aseg->zi[z0+zk] ;

        val = (int)MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
        if (val != label) continue ;

	// If this point is within the label, then look at only the 6 face neighbors
        for (dx = -1 ; dx <= 1 ; dx++)  {
          for (dy = -1 ; dy <= 1 ; dy++)  {
            for (dz = -1 ; dz <= 1 ; dz++)   {
              if (fabs(dx) + fabs(dy) + fabs(dz) != 1) continue ;  // only 8-connected nbrs (??)
                
              xn = mri_aseg->xi[xi+dx] ;
              yn = mri_aseg->yi[yi+dy] ;
              zn = mri_aseg->zi[zi+dz] ;
              val = (int)MRIgetVoxVal(mri_aseg, xn, yn, zn, 0) ;
              if (val != label) {
		// This voxel not the target label but is at the edge of the target label
		// "surface" of label - interface between label and non-label
                nvox++ ;
                if(use_abs){
                  nx += fabs(dx) ;
                  ny += fabs(dy) ;
                  nz += fabs(dz) ;
                }
                else
                {
                  nx += dx ;
                  ny += dy ;
                  nz += dz ;
                }
              }
            }
          }
        }
      }
    }
  }

  if (nvox > 0)
  {
    nx /= nvox ;
    ny /= nvox ;
    nz /= nvox ;
  }

  mag = sqrt(nx*nx + ny*ny + nz*nz) ;
  if (mag > 0)
  {
    nx /= mag ;
    ny /= mag ;
    nz /= mag ;
  }

  *pnx = nx ;
  *pny = ny ;
  *pnz = nz ;
  return(NO_ERROR) ;
}

/*!
  \fn int MRIScopyCoords(MRIS *surf, MRIS *surfcoords)
  \brief Transfers the xyz coords from surfcoords to the given surface.
*/
int MRIScopyCoords(MRIS *surf, MRIS *surfcoords)
{
  int k;

  if(surf == NULL){
    printf("ERROR: MRIScopyCoordsCoords(): surf is null\n");
    return(1);
  }
  if(surfcoords == NULL){
    printf("ERROR: MRIScopyCoordsCoords(): surfcoords is null\n");
    return(1);
  }
  if(surf->nvertices != surfcoords->nvertices){
    printf("ERROR: MRIScopyCoordsCoords(): surf and surfcoords have diff no of vetices\n");
    return(1);
  }

  for(k=0; k < surf->nvertices; k++){
    VERTEX *v = &surf->vertices[k];
    //if(v->ripflag) continue;
    v->x = surfcoords->vertices[k].x;
    v->y = surfcoords->vertices[k].y;
    v->z = surfcoords->vertices[k].z;
  }
  return(0);
}
