/*
 * @file utilities operating on Original
 *
 */
/*
 * surfaces Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2017/02/16 19:43:03 $
 *    $Revision: 1.793 2011
 *
 * $ © copyright-2014 The General Hospital Corporation (Boston, MA) "MGH"
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
int BorderValsHiRes   = 0;
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

/*---------------------------------------------------------------
  MRISurfSrcVersion() - returns CVS version of this file.
  ---------------------------------------------------------------*/
const char *MRISurfSrcVersion(void) { return ("$Id$"); }


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

int MRISaddCommandLine(MRI_SURFACE *mris, char *cmdline)
{
  int i;
  if (mris->ncmds >= MAX_CMDS)
    ErrorExit(ERROR_NOMEMORY, "MRISaddCommandLine: can't add cmd %s (%d)", cmdline, mris->ncmds);

  i = mris->ncmds++;
  mris->cmdlines[i] = (char *)calloc(strlen(cmdline) + 1, sizeof(char));
  strcpy(mris->cmdlines[i], cmdline);
  return (NO_ERROR);
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
