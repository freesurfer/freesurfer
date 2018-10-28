/**
 * @file  mrishash.h
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which 
 * need to compare vertices with one another or to a point.  See: 
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2015/03/18 17:04:00 $
 *    $Revision: 1.26 $
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


#ifndef MRISHASH_ONCE_H
#define MRISHASH_ONCE_H

// define the following to get a single inclusion of non-renamed functions
#define MRISHASH_VANILLA_FUNCS


//--------------------------
typedef struct _mht MRIS_HASH_TABLE, MHT ;

// Ad hoc test functions
int MHT_gw_version(void);  // version of that unit
void MHTfindReportCounts(int * BucketsChecked, 
                         int * BucketsPresent, 
                         int * VtxNumByMHT);
int MHTtestIsMRISselfIntersecting(MRIS const *mris, float res);


void MHTfree(MRIS_HASH_TABLE**);

int MHTwhich(MRIS_HASH_TABLE*);

//------------------------------------------------
// Surface --> MHT, store Face Numbers
//------------------------------------------------

MRIS_HASH_TABLE* MHTcreateFaceTable(
    MRIS const   *mris) ;

MRIS_HASH_TABLE *MHTcreateFaceTable_Resolution(
    MRIS const *mris, 
    int   which, 
    float res) ;

// Add/remove the faces of which vertex V is a part
int  MHTaddAllFaces(   MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX_TOPOLOGY const *v) ;
int  MHTremoveAllFaces(MRIS_HASH_TABLE *mht, MRIS const *mris, VERTEX_TOPOLOGY const *v) ;

//------------------------------------------------
// Surface --> MHT, store Vertex Numbers
//------------------------------------------------
MRIS_HASH_TABLE *MHTcreateVertexTable(
    MRIS const *mris, 
    int which) ;
                                    
MRIS_HASH_TABLE *MHTcreateVertexTable_Resolution(
    MRIS const *mris,
    int which,
    float res) ;

//------------------------------------------------
// Surface self-intersection (Uses MHT initialized with FACES)
//------------------------------------------------
int MHTdoesFaceIntersect(MRIS_HASH_TABLE *mht, MRIS const *mris,int fno);


int MHTisVectorFilled(MRIS_HASH_TABLE const *mht,    MRIS const *mris,
                         int vno,  float dx, float dy, float dz) ;

//------------------------------------------------
// Find nearest vertex/vertices (Uses MHT initialized with VERTICES)
//------------------------------------------------
//------- new generic find function ------------
int MHTfindClosestVertexGeneric(MRIS_HASH_TABLE *mht, 
                                MRIS const *mris,
                                double probex, double probey, double probez,
                                double in_max_distance_mm, 
                                int in_max_halfmhts,
                                VERTEX **pvtx, 
                                int *vtxnum, 
                                double *vtx_distance);

//------- original mrishash find functions ------------

VERTEX *MHTfindClosestVertex  (MRIS_HASH_TABLE *mht, 
                               MRIS const *mris, 
                               VERTEX const *v) ;
int     MHTfindClosestVertexNo(MRIS_HASH_TABLE *mht, 
                               MRIS const *mris, 
                               VERTEX const *v, 
                               float *min_dist);
int     MHTfindClosestVertexNoXYZ(MRIS_HASH_TABLE *mht, 
                               MRIS const *mris, 
                               float x, float y, float z, 
                               float *min_dist);

                             
VERTEX *MHTfindClosestVertexSet(MRIS_HASH_TABLE *mht, 
                                MRIS const *mris, 
                                VERTEX const *v, 
                                int which) ;
VERTEX * MHTfindClosestVertexSetInDirection(MRIS_HASH_TABLE *mht, 
                                            MRIS const *mris, 
                                            VERTEX const *v, 
                                            int which,
                                            double nx, double ny, double nz);
int    *MHTgetAllVerticesWithinDistance(MRIS_HASH_TABLE *mht, 
                                        MRIS const *mris,
                                        int vno, 
                                        float max_dist, 
                                        int *pvnum);
VERTEX *MHTfindClosestVertexInTable(MRIS_HASH_TABLE *mht, 
                                    MRIS const *mris,
                                    float x, float y, float z, int do_global_search) ;

//------------------------------------------------
// Diagnostic
//------------------------------------------------
int MHTcheckFaces(MRIS const *mris,MRIS_HASH_TABLE *mht) ;
int MHTcheckSurface(MRIS const *mris,MRIS_HASH_TABLE *mht);


//------------------------------------------------
// utilities for finding closest face
//------------------------------------------------
int MHTfindClosestFaceGeneric(MRIS_HASH_TABLE *mht, 
                              MRIS const *mris,
                              //---------- inputs --------------
                              double probex, double probey, double probez,
                              // How far to search: set one or both
                              double in_max_distance_mm, /* Use large number 
                                                            to ignore */
                              int    in_max_mhts,  /* Use -1 to ignore */
                              // only faces that projection is interior to (Use -1 to ignore )
                              int    project_into_face, 
                              //---------- outputs -------------
                              FACE **pface, 
                              int *pfno, 
                              double *pface_distance);
                              
int mhtBruteForceClosestFace(MRIS const *mris, 
                             float x, float y, float z, 
                             int which,                  // which surface within mris to search
                             float *dmin);

//-----
// Support parallelism
//-----

void MHT_maybeParallel_begin();     // Note: Can be nested!
void MHT_maybeParallel_end();

#endif // END #ifndef MRISHASH_ONCE_H
