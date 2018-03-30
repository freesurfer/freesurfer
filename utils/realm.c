#include "realm.h"

#ifndef freeAndNULL
#define freeAndNULL(PTR) { free((PTR)); (PTR) = NULL; }
#endif

/**
 * @file  realm.c
 * @brief support quickly scanning all the vertices or faces on an MRI for
 *        for those that might intersect a brick
 *
 * 
 */
/*
 * Original Author: Bevin R Brett
 * CVS Revision Info:
 *    $Author: brbrett $
 *    $Date: 2018/02/21 15:00:00 $
 *    $Revision: 1.0 $
 *
 * Copyright © 2018 The General Hospital Corporation (Boston, MA) "MGH"
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

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <strings.h>

#include "fnv_hash.h"

static int chkBnd(int lo, int b, int hi) {
    costlyAssert(lo <= b);
    costlyAssert(b < hi);
    return b;
}

#ifdef REALM_UNIT_TEST
    //
    // (cd utils; rm -f a.out ; gcc -o a.out -I ../include realm.c -lm |& less ; ./a.out)
    //
    
    #define freeAndNULL(PTRVAR) { free((PTRVAR)); (PTRVAR) = NULL; }

    #define MAX_FACES_PER_VERTEX 3
    
    typedef struct VERTEX {
        float someX,someY,someZ;
        int num;                        // number of faces
        int f[MAX_FACES_PER_VERTEX];    // fno of the faces this vertex touches
    } VERTEX;

    void getSomeXYZ(VERTEX const * vertex, float* x, float* y, float* z) {
        *x = vertex->someX;
        *y = vertex->someY;
        *z = vertex->someZ;
    }
    
    #define VERTICES_PER_FACE 3
    typedef int vertices_per_face_t[VERTICES_PER_FACE];

    typedef struct FACE {
        vertices_per_face_t v;
    } FACE;

    struct MRIS {
      int     nvertices;
      VERTEX* vertices;
      int     nfaces;
      FACE*   faces;
    };
    
    static float MIN(float lhs, float rhs) { return (lhs < rhs) ? lhs : rhs; }
    static float MAX(float lhs, float rhs) { return (lhs > rhs) ? lhs : rhs; }
    
    static int int_compare(const void* lhs_ptr, const void* rhs_ptr) {
        int lhs = *(int*)lhs_ptr;
        int rhs = *(int*)rhs_ptr;
        return lhs - rhs;
    }

    void* qsort_ctx;
    static int vno_compare(const void* lhs_ptr, const void* rhs_ptr) {
        int    lhs = *(int*)lhs_ptr;
        int    rhs = *(int*)rhs_ptr;
        float* ctx = (float*)qsort_ctx;
        
        return ctx[lhs] - ctx[rhs];
    }

    typedef struct PossiblyIntersectingGreatArcs_callback_context {
      int  capacity;
      int  size;
      int* keys;
    } PossiblyIntersectingGreatArcs_callback_context;

    static bool possiblyIntersectingGreatArcs_callback (void* void_ctx, int key, bool* isHit) {
      PossiblyIntersectingGreatArcs_callback_context* ctx = (PossiblyIntersectingGreatArcs_callback_context*)void_ctx;
      if (ctx->size == ctx->capacity) {
        ctx->capacity *= 2; if (!ctx->capacity) ctx->capacity = 64; 
        ctx->keys = (int*)realloc(ctx->keys, ctx->capacity*sizeof(int));
      }
      ctx->keys[ctx->size++] = key;
      *isHit = false;
      return true;  // keep sending them to me
    }

    void test(int nvertices, int useDuplicates) {
        fprintf(stdout,"Test nvertices:%d useDuplicates:%d\n", nvertices, useDuplicates);
        
        int fBenefitCount = 0, fBenefitLimit = 1, fNoBenefitCount = 0, fHasBenefitCount = 0;

        MRIS mris;

        // add the vertices
        //
        mris.nvertices = nvertices;
        mris.vertices = (VERTEX*)calloc(mris.nvertices, sizeof(VERTEX));

        int vno;
        for (vno = 0; vno < mris.nvertices; vno++) {
            int key = vno > useDuplicates ? vno : 936; 
            VERTEX* v = &mris.vertices[vno];
            v->someX = (key*321)%51; 
            v->someY = (key*7321)%71; 
            v->someZ = (key*17321)%91;
        }
        vno = 0;
        float xMin = mris.vertices[vno].someX, xMax = xMin,
              yMin = mris.vertices[vno].someY, yMax = yMin, 
              zMin = mris.vertices[vno].someZ, zMax = zMin;
        for (vno = 1; vno < mris.nvertices; vno++) {
            VERTEX* v = &mris.vertices[vno];
            xMin = MIN(xMin, v->someX);
            yMin = MIN(yMin, v->someY);
            zMin = MIN(zMin, v->someZ);
            xMax = MAX(xMax, v->someX);
            yMax = MAX(yMax, v->someY);
            zMax = MAX(zMax, v->someZ);
        }
        

        // add mostly small faces
        //
        mris.nfaces = (nvertices > 2) ? nvertices - 2 : 0;    // see below
        mris.faces  = (FACE*)calloc(mris.nfaces, sizeof(FACE));

        const float delta_x = MAX(2, (xMax - xMin)/30 );
        const float delta_y = MAX(2, (yMax - yMin)/30 );
        const float delta_z = MAX(2, (zMax - zMin)/30 );
        
        int*   vnos  = (int*  )calloc(nvertices,sizeof(int));
        float* ctx_x = (float*)calloc(nvertices,sizeof(float));
        float* ctx_y = (float*)calloc(nvertices,sizeof(float));
        float* ctx_z = (float*)calloc(nvertices,sizeof(float));
        {
            int i; 
            for (i=0; i<nvertices; i++) {
                vnos [i] = i;
                ctx_x[i] = mris.vertices[i].someX;
                ctx_y[i] = mris.vertices[i].someY;
                ctx_z[i] = mris.vertices[i].someZ;
            }
        }
        
        //  choose a set of close x's
        qsort_ctx = ctx_x;
        qsort(vnos, nvertices, sizeof(int), vno_compare);

        int fno = 0;
        int i = 0;
        while (i+2 < nvertices) {
            int iLo = i; i++;
            float center_x = ctx_x[vnos[iLo]];
            while (i < nvertices && fabs(center_x - ctx_x[vnos[i]]) < delta_x) i++;
            // [iLo..i) will be emitted by the following loop
            
            //  sort the subrange by y
            qsort_ctx = ctx_y;
            qsort(vnos+iLo, i-iLo, sizeof(int), vno_compare);
            
            int j=iLo;
            while (j < i) {
                //  choose a set of close x's with close y's
                int jLo = j; j++;
                float center_y = ctx_y[vnos[jLo]];
                while (j < i && fabs(center_y - ctx_y[vnos[j]]) < delta_y) j++;
                // [jLo..j) will be emitted by the following loop
                
                // sort the sub-sub range by z
                qsort_ctx = ctx_z;
                qsort(vnos+jLo, j-jLo, sizeof(int), vno_compare);
                
                int k=jLo;
                while (k < j) {
                    //  choose a set of close x's with close y's with close z's
                    int kLo = k; k++;
                    float center_z = ctx_z[vnos[kLo]];
                    while (k < j && fabs(center_z - ctx_z[vnos[k]]) < delta_z) k++;
                    // [kLo..k) will be emitted by the following loop
                    
                    // make the faces
                    int m;
                    for (m = kLo; m+2 < k; m++) {
                        int v0 = m;
                        int v1 = m+1; if (v1 > j) v1 -= j-jLo;
                        int v2 = m+2; if (v2 > j) v2 -= j-jLo;
                        if (fno >= mris.nfaces) *(int*)-1 = 0;
                        FACE* face = &mris.faces[fno];
                        face->v[0] = vnos[v0];
                        face->v[1] = vnos[v1]; 
                        face->v[2] = vnos[v2];
                        int vi;
                        for (vi = 0; vi < 3; vi++) {
                            int vno = face->v[vi];
                            VERTEX* v = &mris.vertices[vno];
                            if (v->num >= MAX_FACES_PER_VERTEX) *(int*)-1 = 0;;
                            v->f[v->num++] = fno;
                        }
                        fno++;
                    }
                }
            }
        }
        mris.nfaces = fno;  // decrease to the correct number

        freeAndNULL(vnos );               
        freeAndNULL(ctx_x);
        freeAndNULL(ctx_y);
        freeAndNULL(ctx_z);
        
#if 1          
        if (1) {
            fprintf(stderr, "Testing GreatArcSet\n");
            
            int  size;
            int* keys;

            if (1) {
                // Make a GreatArcSet
                //
                GreatArcSet* gas = makeGreatArcSet(&mris);

                // Throw in some one edge
                //
                mris.vertices[0].someX = 1;         // vertex 0
                mris.vertices[0].someY = 0;
                mris.vertices[0].someZ = -1.1;
                mris.vertices[1].someX = 1;         // vertex 1
                mris.vertices[1].someY = 0;
                mris.vertices[1].someZ = -1.2;
                insertGreatArc(gas, 77, 0,1);       // edge connecting them
                
                // See if an intersection is detected
                //
                {
                    PossiblyIntersectingGreatArcs_callback_context context;
                    bzero(&context, sizeof(context));
                    possiblyIntersectingGreatArcs(gas, &context, possiblyIntersectingGreatArcs_callback, 1,0,-1.2, 1,0,-1.1, false);
                    int size  = context.size;
                    int* keys = context.keys; 
                    if (size    !=  1) fprintf(stderr, "same lines failed, size:%d\n", size); else
                    if (keys[0] != 77) fprintf(stderr, "same lines wrong key\n"); else
                                       fprintf(stderr, "same lines correct key - good!\n");
                    freeAndNULL(keys);
                }
                
                // Done
                //
                freeGreatArcSet(&gas);
            }
            
            if (1) {
            
                // Make a GreatArcSet
                //
                GreatArcSet* gas = makeGreatArcSet(&mris);

                // Throw in edges to make squared graph paper
                //
                int key = 0;
                int vno = 0;
                int x,y;
                for (x = 0; x < 100; x++)             
                for (y = 0; y < 100; y++) {
                    if (vno + 3 >= mris.nvertices) break;
                               
                    mris.vertices[vno].someX = x;
                    mris.vertices[vno].someY = y;
                    mris.vertices[vno].someZ = 1000.0;
                    vno++;
                    mris.vertices[vno].someX = x+1;
                    mris.vertices[vno].someY = y;
                    mris.vertices[vno].someZ = 1000.0;
                    vno++;
                    mris.vertices[vno].someX = x;
                    mris.vertices[vno].someY = y+1;
                    mris.vertices[vno].someZ = 1000.0;
                    vno++;

                    insertGreatArc(gas, key++, vno-3,vno-2);
                    insertGreatArc(gas, key++, vno-3,vno-1);
                }

                // Lookup a variety of lines and make sure that at least the right ones are found
                // and not too many others
                //                
                for (x = 5; x < 100; x+=5)             
                for (y = 5; y < 100; y+=5) {
                    PossiblyIntersectingGreatArcs_callback_context context;
                    bzero(&context, sizeof(context));

                    possiblyIntersectingGreatArcs(gas, &context, possiblyIntersectingGreatArcs_callback,  x-0.5,y-0.5,1000.0, x+x%5,y+y%3,1000.0, true);

                    size = context.size;
                    keys = context.keys;
                    fprintf(stderr, "%d possible intersections found near %d\n", size, x*10000+y*100+0);
                    int i;
                    for (i = 0; i < size; i++) {
                        if (i == 40) {
                            fprintf(stderr, " ...");
                            break;
                        }
                        fprintf(stderr, "  %d", keys[i]);
                    }
                    fprintf(stderr, "\n");
                    freeAndNULL(keys);
                }
            
                // Done
                //
                freeGreatArcSet(&gas);
            }
        }
#endif
        
        if (0) {
            fprintf(stdout, "Testing RealmTree\n");
        
            // Make a RealmTree
            //
            RealmTree* realmTree = makeRealmTree(&mris, getSomeXYZ);
            if (0) summarizeRealmTree(realmTree);

            // Move some vertices
            //
            {
                int vno;
                for (vno = 0; vno < mris.nvertices; vno++) {
                    if (vno % 77 >= 3) continue;
                    VERTEX* v = &mris.vertices[vno];
                    v->someX += 0.1 * (xMax - v->someX);
                    v->someY += 0.1 * (yMax - v->someY);
                    v->someZ += 0.1 * (zMax - v->someZ);
                    noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                }
                for (vno = 0; vno < mris.nvertices; vno++) {
                    if ((vno * 123) % 31 >= 2) continue;
                    VERTEX* v = &mris.vertices[vno];
                    v->someX -= 0.1 * (v->someX - xMin);
                    v->someY -= 0.1 * (v->someY - yMin);
                    v->someZ -= 0.1 * (v->someZ - zMin);
                    noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                }
                fprintf(stdout,"Moved vertices, now checking\n");
                checkRealmTree(realmTree, &mris, getSomeXYZ);
                fprintf(stdout,"Checked realmTree now updating\n");
                updateRealmTree(realmTree, &mris, getSomeXYZ);
                fprintf(stdout,"Updated realmTree, now checking\n");
                checkRealmTree(realmTree, &mris, getSomeXYZ);
                fprintf(stdout,"Checked realmTree\n");
            }

            // Be nasty, deliberately go outside in all the different directions
            //
            if (mris.nvertices >= 6) {
                int vno; VERTEX* v;
                vno = 0; v = &mris.vertices[vno]; v->someX = xMin - 0.1; noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                                                                         updateRealmTree          (realmTree, &mris, getSomeXYZ);
                vno = 1; v = &mris.vertices[vno]; v->someY = yMin - 0.1; noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                                                                         updateRealmTree          (realmTree, &mris, getSomeXYZ);
                vno = 2; v = &mris.vertices[vno]; v->someZ = zMin - 0.1; noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                                                                         updateRealmTree          (realmTree, &mris, getSomeXYZ);
                vno = 3; v = &mris.vertices[vno]; v->someX = xMax + 0.1; noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                                                                         updateRealmTree          (realmTree, &mris, getSomeXYZ);
                vno = 4; v = &mris.vertices[vno]; v->someY = yMax + 0.1; noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                                                                         updateRealmTree          (realmTree, &mris, getSomeXYZ);
                vno = 5; v = &mris.vertices[vno]; v->someZ = zMax + 0.1; noteIfXYZChangedRealmTree(realmTree, &mris, getSomeXYZ, vno);
                                                                         updateRealmTree          (realmTree, &mris, getSomeXYZ);
            }

            // Check varous realms
            //
            int fLimit = 1;
            int fCount = 0;
            float xfLo, xfHi;
            float yfLo, yfHi;
            float zfLo, zfHi;
            for (xfLo = -0.1; xfLo <= 1.2; xfLo += 0.1)     // check also when the realm exceeds the original bounds
            for (xfHi = -0.1; xfHi <= 1.2; xfHi += 0.1)     // because this can happen...
            for (yfLo = -0.1; yfLo <= 1.2; yfLo += 0.1)
            for (yfHi = -0.1; yfHi <= 1.2; yfHi += 0.1)
            for (zfLo = -0.1; zfLo <= 1.2; zfLo += 0.1)
            for (zfHi = -0.1; zfHi <= 1.2; zfHi += 0.1)
            {
                float xLo = xMin +    xfLo *(xMax-xMin);
                float xHi = xMax - (1-xfHi)*(xMax-xLo);
                float yLo = yMin +    yfLo *(yMax-yMin);
                float yHi = yMax - (1-yfHi)*(yMax-yLo);
                float zLo = zMin +    zfLo *(zMax-zMin);
                float zHi = zMax - (1-zfHi)*(zMax-zLo);

                fCount++;
                if (fCount == fLimit) {
                    fLimit *= 2;
                    fprintf(stdout,"fCount:%d x:%f..%f y:%f.%f z:%f..%f\n", fCount, xLo, xHi, yLo, yHi, zLo, zHi);
                }

                Realm* realm = 
                    makeRealm(realmTree, 
                        xLo, xHi, 
                        yLo, yHi,
                        zLo, zHi);

                RealmIterator realmIterator;
                initRealmIterator(&realmIterator, realm);

                int* states = (int*)calloc(mris.nvertices, sizeof(int));

                int counter = 1;
                int vno;
                for (;;) {
    #ifdef REALM_UNIT_TEST
                    if (false && (counter == 1 || counter == 122)) {
                        fprintf(stdout,"counter:%d ri.i:%ld ri.p:%p\n", counter, realmIterator.i, realmIterator.p); 
                    }
    #endif
                    vno = realmNextMightTouchVno(realm, &realmIterator);
    #ifdef REALM_UNIT_TEST
                    if (vno < -1 || vno >= mris.nvertices) {
                        fprintf(stdout,"ERROR, vno:%d is illegal\n", vno); 
                        exit(1);
                    }
    #endif
                    if (0 > vno) break;
                    if (counter == 0 || states[vno]) {
                        fprintf(stdout,"ERROR, vno:%d reported again when counter:%d, was reported counter:%d\n", vno, counter, states[vno]); 
                        exit(1);
                    }
                    states[vno] = counter++;
                }

                // No vno should have been visited more than once
                // No unreported vno should be in the region
                for (vno = 0; vno < mris.nvertices; vno++) {
                    if (states[vno] > 1 ) 
                    if (states[vno] == 0) {
                       VERTEX* v = &mris.vertices[vno];
                       if (xLo <= v->someX && v->someX < xHi 
                       &&  yLo <= v->someY && v->someY < yHi
                       &&  zLo <= v->someZ && v->someZ < zHi) fprintf(stdout,"ERROR, vno:%d was not reported\n", vno);
                    }
                }

                // Check that at least the needed fno's are reported and that none is reported twice
                //            
                int  fnosCapacity = realmNumberOfMightTouchFno(realm);
                int* fnos         = (int*)calloc(fnosCapacity, sizeof(int));
                int  fnosSize     = realmMightTouchFno(realm, fnos, fnosCapacity);

                qsort(fnos, fnosSize, sizeof(int), int_compare);

                int fnosI;
                for (fnosI = 0; fnosI < fnosSize-1; fnosI++) {
                    if (fnos[fnosI] >= fnos[fnosI + 1]) {
                        fprintf(stdout,"ERROR, fnos[fnosI]:%d fnos[fnosI+1]:%d\n", fnos[fnosI], fnos[fnosI + 1]);
                    }
                }

                fnosI = 0;
                int fno;
                for (fno = 0; fno < mris.nfaces; fno++) {
                    FACE const * face = &mris.faces[fno];
                    int vi = 0;
                    VERTEX const * vertex = &mris.vertices[face->v[vi]];
                    float fxLo = vertex->someX, fxHi = fxLo,
                          fyLo = vertex->someY, fyHi = fyLo,
                          fzLo = vertex->someZ, fzHi = fzLo;
                    for (vi = 0; vi < VERTICES_PER_FACE; vi++) {
                        fxLo = MIN(fxLo, vertex->someX); fxHi = MAX(fxHi, vertex->someX);
                        fyLo = MIN(fyLo, vertex->someY); fyHi = MAX(fyHi, vertex->someY);
                        fzLo = MIN(fzLo, vertex->someZ); fzHi = MAX(fzHi, vertex->someZ);
                    }
                    bool wontIntersect =  
                        fxHi < xLo || xHi <= fxLo ||
                        fyHi < yLo || yHi <= fyLo ||
                        fzHi < zLo || zHi <= fzLo;
                    if (wontIntersect) continue;                            // might or might not be in the list
                    while (fnosI < fnosSize && fnos[fnosI] < fno) fnosI++;  // skip the ones that were reported but need not be
                    if (fnosI == fnosSize || fnos[fnosI] != fno) {
                        fprintf(stdout,"ERROR, fno:%d was not reported\n", fno);
                    }
                }

                // We are only interested in the benefits when the realm is much smaller than the volume
                //
                if (mris.nfaces > 0 &&
                    (xHi - xLo) < (xMax - xMin)/4 &&
                    (yHi - yLo) < (yMax - yMin)/4 &&
                    (zHi - zLo) < (zMax - zMin)/4
                    ) {

                    if (fnosSize*3 > mris.nfaces*2) fNoBenefitCount++; else fHasBenefitCount++;

                    if (++fBenefitCount == fBenefitLimit) {
                        if (fBenefitLimit < 1000) fBenefitLimit *= 2; else fBenefitLimit += 1000;
                        fprintf(stdout,"fnosSize:%d mris.nfaces:%d fNoBenefitCount:%d fHasBenefitCount:%d\n", 
                            fnosSize, mris.nfaces, fNoBenefitCount, fHasBenefitCount);
                    }
                }

                // Done
                //
                freeAndNULL(fnos);
                freeAndNULL(states);
                freeRealm(&realm);
            }

            freeRealmTree(&realmTree);
        }
                
        freeAndNULL(mris.vertices);
        freeAndNULL(mris.faces); 
    }
    
    int main() {
        int useDuplicates;
        for (useDuplicates = 0; useDuplicates <= 1000; useDuplicates += 200) {
            if (0) test(1000, useDuplicates);
            if (0) test(1, useDuplicates);
            if (0) test(2, useDuplicates);
            if (0) test(3, useDuplicates);
            if (0) test(4, useDuplicates);
            if (0) test(100, useDuplicates);
            if (0) test(10000, useDuplicates);
            if (1) test(100000, useDuplicates);
            if (0) test(0, useDuplicates);
        }
        return 0;
    }
#endif



// RealmTree
//
typedef struct RealmTreeNode RealmTreeNode;
struct RealmTreeNode {
    float xLo, xHi, yLo, yHi, zLo, zHi;
    RealmTreeNode*  parent;
    int             depth;    
#define childrenSizeLog2 3                              // 2x 2y 2z
#define childrenSize     (1<<childrenSizeLog2)      
#define maxVnosSizeLog2  20                             // only support 1M vno's
#define vnosBuffSize     ((sizeof(RealmTreeNode*)*childrenSize/sizeof(int)) - 2)    // 2 for vnosSize and vnosCapacity
    int* vnos;                                          // NULL for non-leaf nodes, either &vnosBuff or 
    union {
        RealmTreeNode*  children[childrenSize];
        struct {
            int         vnosSize;                       // above 2 assumes these are the same as the vnosBuff elements
            int         vnosCapacity;
            int         vnosBuff[vnosBuffSize];
        };
    };
    int nFaces;                                         // the number of faces in the following list
    int firstFnoPlus1;                                  // the first of a list of faces for which this node is the deepest node they fully fit within
};
static const unsigned long childIndexBits =      childrenSizeLog2;
static const unsigned long childIndexMask = ((1<<childrenSizeLog2) - 1);
static const unsigned long leafIndexBits  =      maxVnosSizeLog2;
static const unsigned long leafIndexMask  = ((1<<maxVnosSizeLog2 ) - 1);

typedef struct Captured_VERTEX_xyz {
    float x,y,z;
} Captured_VERTEX_xyz;

struct RealmTree {
    MRIS const  *           mris;
    int                     saved_nvertices;    // detect if these mris change
    int                     saved_nfaces;
    
    Captured_VERTEX_xyz*    captured_VERTEX_xyz;
    RealmTreeNode**         vnoToRealmTreeNode;
    RealmTreeNode**         fnoToRealmTreeNode;
    
    // links in chains of fno off each RealmTreeNode
    //      nextFnoPlus1[fno]         == 0   means end of chain
    //  
    int*                    nextFnoPlus1;

    // links in chains of vno's whose update is pending
    //      nextVnoToUpdatePlus1[vno] ==  0  means vno not in chain
    //                                   -1  means it is the last in the chain
    int                     firstVnoToUpdatePlus1;
    int*                    nextVnoToUpdatePlus1;

    // links in chains of fno's whose update are pending
    // this is only used during updateRealmTree
    //      kept here to avoid need to reallocate
    //      full of zero's between uses
    //
    //      nextFnoToUpdatePlus1[fno] ==  0  means fno not in chain
    //                                   -1  means it is the last in the chain
    int*                    nextFnoToUpdatePlus1;

    // the root node of the tree of nodes
    // the tree has 8 children below each parent, being a 2x 2y 2z
    //  
    RealmTreeNode           root;

    // Some debugging support
    RealmTreeNode* interestingRealmTreeNode;
};


static void constructRealmTreeNode(RealmTreeNode *child, RealmTreeNode *parent) {
    child->parent = parent;
    child->depth  = parent ? parent->depth+1 : 0;
    child->vnos   = child->vnosBuff;
    child->vnosCapacity = vnosBuffSize;
}

static void destroyRealmTreeNode(RealmTreeNode *n) {
    if (!n->vnos) {
        int c;
        for (c = 0; c < childrenSize; c++) {
            RealmTreeNode * child = n->children[c]; n->children[c] = NULL;
            if (!child) continue;
            destroyRealmTreeNode(child);
            freeAndNULL(child);
        }
    } else {
        if (n->vnos != n->vnosBuff) freeAndNULL(n->vnos);
    }
}

static const int maxDepth = 
    //
    // The maxDepth that can be reached when there are many (x,y,z) VERY close to each other and a few a long way away
    //
    (   sizeof(((RealmIterator*)NULL)->i)*8     // available                        64
      - maxVnosSizeLog2                         // needed to index the leaf nodes   20, leaving 44
      - 1                                       // needed to mark top of search      1, leaving 43
    )
    / childrenSizeLog2;                         // bits needed per non-leaf level   43/3 = 14 - more than enough

static bool nodeContains(
    RealmTreeNode const *  const n, 
    float const x, float const y, float const z) {
    return  n->xLo <= x && x < n->xHi &&
            n->yLo <= y && y < n->yHi &&
            n->zLo <= z && z < n->zHi;
}

static RealmTreeNode const * upUntilContainsNode(RealmTreeNode const * n, 
    float const x, float const y, float const z) {
    while (n && !nodeContains(n,x,y,z)) n = n->parent;
    return n;
}

static RealmTreeNode* deepestCommonNode(RealmTreeNode* n1, RealmTreeNode* n2) {
    while (n1->depth > n2->depth) n1 = n1->parent;
    while (n1->depth < n2->depth) n2 = n2->parent;
    while (n1 != n2) { n1 = n1->parent; n2 = n2->parent; }
    return n1;
}

static int chooseChild(
    RealmTreeNode const * n,
    float x, float y,float z) {

#ifdef REALM_UNIT_TEST
    if (!nodeContains(n, x,y,z)) 
        *(int*)-1 = 0;
#endif

    float xMid = n->children[1]->xLo;
    float yMid = n->children[2]->yLo;
    float zMid = n->children[4]->zLo;
    
    int c = ((x < xMid) ? 0 : 1) + ((y < yMid) ? 0 : 2) + ((z < zMid) ? 0 : 4);
    
#ifdef REALM_UNIT_TEST
    if (!nodeContains(n->children[c], x,y,z)) 
        *(int*)-1 = 0;
#endif
  
    return c;
}

typedef enum Widen_Mask {
    Widen_xLo =  1, Widen_xHi = 2,
    Widen_yLo =  4, Widen_yHi = 8,
    Widen_zLo = 16, Widen_zHi = 32 } Widen_Mask;
     
static void widenSubtree_wkr(RealmTreeNode* n, float xLo, float xHi, float yLo, float yHi, float zLo, float zHi, Widen_Mask widen_Mask) {
    if (widen_Mask == 0) return;
    //
    if (widen_Mask & Widen_xLo) n->xLo = xLo; if (widen_Mask & Widen_xHi) n->xHi = xHi;
    if (widen_Mask & Widen_yLo) n->yLo = yLo; if (widen_Mask & Widen_yHi) n->yHi = yHi;
    if (widen_Mask & Widen_zLo) n->zLo = zLo; if (widen_Mask & Widen_zHi) n->zHi = zHi;
    
    if (n->vnos) return;    // leaf nodes

    // there are eight children
    //      zyx     
    //    0 000
    //    1 001
    //    2 010
    //    3 011
    //    4 100
    //    5 101
    //    6 110
    //    7 111
    // each should only have some of its bounds widened
    //
    widenSubtree_wkr(n->children[0], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xHi | Widen_yHi | Widen_zHi));    
    widenSubtree_wkr(n->children[1], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xLo | Widen_yHi | Widen_zHi));    
    widenSubtree_wkr(n->children[2], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xHi | Widen_yLo | Widen_zHi));    
    widenSubtree_wkr(n->children[3], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xLo | Widen_yLo | Widen_zHi));    
    widenSubtree_wkr(n->children[4], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xHi | Widen_yHi | Widen_zLo));    
    widenSubtree_wkr(n->children[5], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xLo | Widen_yHi | Widen_zLo));    
    widenSubtree_wkr(n->children[6], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xHi | Widen_yLo | Widen_zLo));    
    widenSubtree_wkr(n->children[7], xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask & ~(Widen_xLo | Widen_yLo | Widen_zLo));    
} 

static void widenSubtree(RealmTree* realmTree, float xLo, float xHi, float yLo, float yHi, float zLo, float zHi) {
    RealmTreeNode* n = &realmTree->root;
    Widen_Mask widen_Mask = 0;
    if (xLo < n->xLo) widen_Mask |= Widen_xLo;
    if (xHi > n->xHi) widen_Mask |= Widen_xHi;
    if (yLo < n->yLo) widen_Mask |= Widen_yLo;
    if (yHi > n->yHi) widen_Mask |= Widen_yHi;
    if (zLo < n->zLo) widen_Mask |= Widen_zLo;
    if (zHi > n->zHi) widen_Mask |= Widen_zHi;
    widenSubtree_wkr(n, xLo, xHi, yLo, yHi, zLo, zHi, widen_Mask);
}

static RealmTreeNode const * deepestContainingNode(RealmTreeNode const * n, float const x, float const y, float const z) {
    n = upUntilContainsNode(n, x, y, z);
    while (n && !n->vnos) {
        int c = chooseChild(n, x, y, z);
        n = n->children[c];
    }
    return n;
}

static RealmTreeNode* insertVnoIntoNode(
    RealmTree*      const realmTree,
    RealmTreeNode*  const n, 
    int const vno);

static RealmTreeNode* insertIntoChild(
    RealmTree*     realmTree,
    RealmTreeNode* n,
    int vno) {
    chkBnd(0, vno, realmTree->saved_nvertices);

    Captured_VERTEX_xyz const * const captured_xyz = &realmTree->captured_VERTEX_xyz[vno];
    float const x = captured_xyz->x, y = captured_xyz->y, z = captured_xyz->z;
    int c = chooseChild(n, x, y, z);
    return insertVnoIntoNode(realmTree, n->children[c], vno);
}

static RealmTreeNode* insertVnoIntoNode(
    RealmTree*      const realmTree,
    RealmTreeNode*  const n, 
    int const vno)
{
    chkBnd(0, vno, realmTree->saved_nvertices);

#ifdef REALM_UNIT_TEST
    Captured_VERTEX_xyz const * const captured_xyz = &realmTree->captured_VERTEX_xyz[vno];
    float const x = captured_xyz->x, y = captured_xyz->y, z = captured_xyz->z;

    MRIS const* mris = realmTree->mris;
    VERTEX const* v = &mris->vertices[vno];
    if (x != v->someX || y != v->someY || z != v->someZ) 
        fprintf(stderr, "vertex moved\n");
#endif
    
    // If this is a leaf
    //
    if (n->vnos) {
        
        // Must extend if full and can't split
        //
        if (n->vnosSize == n->vnosCapacity && n->depth+1 == maxDepth) {
            n->vnosCapacity *= 2;
            int* p = (int*)calloc(n->vnosCapacity, sizeof(int));
            int i;
            for (i = 0; i < n->vnosSize; i++) p[i] = n->vnos[i];
            if (n->vnos != n->vnosBuff) freeAndNULL(n->vnos);
            n->vnos = p;
        }
        
        // Can insert 
        //
        if (n->vnosSize < n->vnosCapacity) {
#ifdef REALM_UNIT_TEST    
            if (!nodeContains(n, x,y,z)) 
                *(int*)-1 = 0;
#endif  
            n->vnos[n->vnosSize++] = vno;
            realmTree->vnoToRealmTreeNode[vno] = n;
            return n;
        }
        
        // Must split

        // Chose the splitting values
        float xMid = (n->xLo + n->xHi)/2;
        float yMid = (n->yLo + n->yHi)/2;
        float zMid = (n->zLo + n->zHi)/2;

        // Save the vnos, since the next step overwrites them
        //
        int const vnosSize = n->vnosSize;   // n->vnosSize and n->vnosBuf will get overwritten by children
        int vnos[vnosBuffSize];
#ifdef REALM_UNIT_TEST    
        if (vnosSize > vnosBuffSize || n->vnos != n->vnosBuff) {
            *(int*)-1 = 0;
        }
#endif  
        int vi;
        for (vi = 0; vi < vnosSize; vi++) {
            vnos[vi] = n->vnos[vi];
        }
        n->vnos = NULL;

        // Make the children
        int c;
        for (c = 0; c < childrenSize; c++) { 
            RealmTreeNode* child = n->children[c] = 
                (RealmTreeNode*)calloc(1, sizeof(RealmTreeNode));
            constructRealmTreeNode(child, n);
#ifdef REALM_UNIT_TEST    
            if (child->depth >= maxDepth) *(int*)-1 = 0;
#endif
            if (c&1) child->xLo = xMid, child->xHi = n->xHi; else child->xLo = n->xLo, child->xHi = xMid; 
            if (c&2) child->yLo = yMid, child->yHi = n->yHi; else child->yLo = n->yLo, child->yHi = yMid;
            if (c&4) child->zLo = zMid, child->zHi = n->zHi; else child->zLo = n->zLo, child->zHi = zMid;
        }
        
        // Insert the saved vno into their child
        for (vi = 0; vi < vnosSize; vi++) {
            insertIntoChild(realmTree, n, vnos[vi]);
        }
    }

    // Insert this vno into the right child
    //
    return insertIntoChild(realmTree, n, vno);
}

static void removeVnoFromRealmTree(RealmTree* realmTree, int vno) {

    RealmTreeNode* n = realmTree->vnoToRealmTreeNode[chkBnd(0, vno, realmTree->saved_nvertices)];
    realmTree->vnoToRealmTreeNode[vno] = NULL;
    {
        int vi;
        for (vi = n->vnosSize - 1; n->vnos[vi] != vno; vi--);                   // find it, backwards since the active ones are at end
        do { n->vnos[vi] = n->vnos[vi+1]; } while (++vi < n->vnosSize - 1);     // remove it
    }
}


static RealmTreeNode* insertVnoNear(RealmTree* realmTree, RealmTreeNode* n, int vno) {
    chkBnd(0, vno, realmTree->saved_nvertices);
    Captured_VERTEX_xyz* captured_xyz = &realmTree->captured_VERTEX_xyz[vno];
    // Find the right subtree
    while (!nodeContains(n, captured_xyz->x,captured_xyz->y,captured_xyz->z)) {
        n = n->parent;
    }
    // Insert here, or deeper
    return insertVnoIntoNode(realmTree, n, vno);
}

static int countXYZChanges(RealmTree const * realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    int count = 0;
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        float x,y,z;
        getXYZ(vertex, &x, &y, &z);
        Captured_VERTEX_xyz* c = &realmTree->captured_VERTEX_xyz[vno];
        if (x != c->x || y != c->y || z != c->z) count++;
    }
    return count;
}


void freeRealmTree(RealmTree** realmTreePtr) {
    RealmTree* rt = *realmTreePtr; *realmTreePtr = NULL;
    if (!rt) return;
    
    freeAndNULL(rt->nextFnoToUpdatePlus1);
    freeAndNULL(rt->nextVnoToUpdatePlus1);
    freeAndNULL(rt->nextFnoPlus1);
    freeAndNULL(rt->fnoToRealmTreeNode);
    freeAndNULL(rt->vnoToRealmTreeNode);
    destroyRealmTreeNode(&rt->root);
    
    freeAndNULL(rt);
}

static float widenHi(float hi) {
    float step = FLT_MIN;
    while (hi + step == hi) {
        step *= 2.0;
    }
    return hi + step;
}

static RealmTreeNode * chooseRealmTreeNodeForFno(MRIS const * const mris, RealmTree const * const rt, int const fno) {
    chkBnd(0, fno, rt->saved_nfaces);

    FACE const * face = &mris->faces[fno];

    RealmTreeNode * n;
    RealmTreeNode * vertexNode;
    int vi,vno;

    vi = 0; 
        vno = face->v[vi]; 
        vertexNode = rt->vnoToRealmTreeNode[vno];
        n = vertexNode;

    for (vi = 1; vi < VERTICES_PER_FACE; vi++) { 
        vno = face->v[vi]; 
        vertexNode = rt->vnoToRealmTreeNode[vno];
        n = deepestCommonNode(n, vertexNode);
    }
    
    return n;
}

static const int interestingFno = -1; // 301539;

static bool isFnoInRealmTreeNode(RealmTree* realmTree, int fno) {
    RealmTreeNode* n = realmTree->fnoToRealmTreeNode[chkBnd(0,fno,realmTree->saved_nfaces)];

    int prevFno = -1;                                           // this list is usually very small
    int entryFno = n->firstFnoPlus1 - 1;                        // should this search be a performance problem
    while (entryFno >= 0 && entryFno != fno) {                  //      change to a per-node btree
        prevFno = entryFno;
        entryFno = realmTree->nextFnoPlus1[chkBnd(0,prevFno,realmTree->saved_nfaces)] - 1;
    }
    return (entryFno == fno);
}

static void insertFnoIntoRealmTreeNode(RealmTree* realmTree, RealmTreeNode* n, int fno) {
    realmTree->fnoToRealmTreeNode[chkBnd(0,fno,realmTree->saved_nfaces)] = n;
    realmTree->nextFnoPlus1[fno]       = n->firstFnoPlus1;
    n->firstFnoPlus1                   = fno + 1;
    // adjust the count
    n->nFaces++;

    if (false && n == realmTree->interestingRealmTreeNode) { 
        fprintf(stdout, "%s:%d interestingRealmTreeNode inserted into\n", __FILE__, __LINE__);
    }
    
    if (fno == interestingFno) { 
        fprintf(stdout, "%s:%d interestingFno inserted\n", __FILE__, __LINE__);
        realmTree->interestingRealmTreeNode = n;
    }
    
    if (realmTree->interestingRealmTreeNode) costlyAssert(isFnoInRealmTreeNode(realmTree, interestingFno));
}

static void removeFnoFromRealmTree(RealmTree* realmTree, int fno) {
    RealmTreeNode* n = realmTree->fnoToRealmTreeNode[chkBnd(0,fno,realmTree->saved_nfaces)];

    if (false && n == realmTree->interestingRealmTreeNode) { 
        fprintf(stdout, "%s:%d interestingRealmTreeNode removed from\n", __FILE__, __LINE__);
    }

    // find in the list
    int prevFno = -1;                                           // this list is usually very small
    int entryFno = n->firstFnoPlus1 - 1;                        // should this search be a performance problem
    while (entryFno != fno) {                                   //      change to a per-node btree
        prevFno = entryFno;
        entryFno = realmTree->nextFnoPlus1[chkBnd(0,prevFno,realmTree->saved_nfaces)] - 1;
    }

    // remove from the list
    realmTree->fnoToRealmTreeNode[fno] = NULL;
    if (prevFno < 0) {
        n->firstFnoPlus1 = realmTree->nextFnoPlus1[fno];
    } else {
        realmTree->nextFnoPlus1[prevFno] = realmTree->nextFnoPlus1[fno];
    }

    // be tidy...
    realmTree->nextFnoPlus1[fno] = 0;

    // adjust the count
    n->nFaces--;

    if (fno == interestingFno) { 
        static long count;
        count++;
        fprintf(stdout, "%s:%d interestingFno removed, count:%ld\n", __FILE__, __LINE__, count);
        realmTree->interestingRealmTreeNode = NULL;
        
        if (count == 248) {
            fprintf(stdout, "%s:%d breakpoint here\n", __FILE__, __LINE__);
        }
    }

    if (realmTree->interestingRealmTreeNode) costlyAssert(isFnoInRealmTreeNode(realmTree, interestingFno));
}


static void resizeRealmTree(RealmTree* rt, MRIS const * mris) {
    
    static int count;
    count++;

    {
        int change = mris->nvertices - rt->saved_nvertices;
        if (change > 0) {   
            rt->captured_VERTEX_xyz = (Captured_VERTEX_xyz*)realloc(rt->captured_VERTEX_xyz,  mris->nvertices*sizeof(Captured_VERTEX_xyz));
            rt->vnoToRealmTreeNode  = (RealmTreeNode**     )realloc(rt->vnoToRealmTreeNode,   mris->nvertices*sizeof(RealmTreeNode*     ));
            rt->nextVnoToUpdatePlus1 =  (int*              )realloc(rt->nextVnoToUpdatePlus1, mris->nvertices*sizeof(int                ));
            
            bzero(rt->captured_VERTEX_xyz + rt->saved_nvertices, change*sizeof(Captured_VERTEX_xyz));
            bzero(rt->vnoToRealmTreeNode  + rt->saved_nvertices, change*sizeof(RealmTreeNode*));
            bzero(rt->nextVnoToUpdatePlus1+ rt->saved_nvertices, change*sizeof(int));
        } else {
            // remove these from the tree
            int vno;
            for (vno = mris->nvertices; vno < rt->saved_nvertices; vno++) {
                RealmTreeNode* const n = rt->vnoToRealmTreeNode[vno]; if (!n) continue;
                removeVnoFromRealmTree(rt, vno);
            }
            // remove any in the vnoToUpdate list
            int* prevLink = &rt->firstVnoToUpdatePlus1;
            while ((vno = *prevLink - 1) >= 0) {
                if (vno >= mris->nvertices) {                           // should remain?
                    *prevLink = rt->nextVnoToUpdatePlus1[vno];          // no - remove from chain
                    rt->nextVnoToUpdatePlus1[vno] = 0;
                } else {
                    prevLink = &rt->nextVnoToUpdatePlus1[vno];          // yes - keep it
                }
            }
        }
        rt->saved_nvertices = mris->nvertices;
    }
    
    {
        int change =  mris->nfaces - rt->saved_nfaces;
        if (change > 0) {   
            rt->fnoToRealmTreeNode   = (RealmTreeNode**)realloc(rt->fnoToRealmTreeNode,   mris->nfaces*sizeof(RealmTreeNode*));
            rt->nextFnoPlus1         = (int*           )realloc(rt->nextFnoPlus1,         mris->nfaces*sizeof(int           ));
            rt->nextFnoToUpdatePlus1 = (int*           )realloc(rt->nextFnoToUpdatePlus1, mris->nfaces*sizeof(int           ));
             
            bzero(rt->fnoToRealmTreeNode   + rt->saved_nfaces, change*sizeof(RealmTreeNode*));
            bzero(rt->nextFnoPlus1         + rt->saved_nfaces, change*sizeof(int));
            bzero(rt->nextFnoToUpdatePlus1 + rt->saved_nfaces, change*sizeof(int));
        } else {
            // remove these from the tree
            int fno;
            for (fno = mris->nfaces; fno < rt->saved_nfaces; fno++) {
                RealmTreeNode* const n = rt->fnoToRealmTreeNode[fno]; if (!n) continue;
                removeFnoFromRealmTree(rt, fno);
            }
        }
        rt->saved_nfaces = mris->nfaces;
    }
}


RealmTree* makeRealmTree(MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    // Fills in the tree using the existing position of 
    // the vertices and faces
    RealmTree* rt = (RealmTree*)calloc(1, sizeof(RealmTree));
    constructRealmTreeNode(&rt->root, NULL);
    
    rt->mris = mris;
    
    resizeRealmTree(rt, mris);

    rt->firstVnoToUpdatePlus1 = -1; // end of list marker, since none pending

    if (mris->nvertices == 0) return rt;
    
    // Capture the xyz and calculate the outer box
    //
    int vno = 0;
    VERTEX       const * vertex0      = &mris->vertices[vno];
    Captured_VERTEX_xyz* captured_xyz = &rt->captured_VERTEX_xyz[vno];
    getXYZ(vertex0, &captured_xyz->x, &captured_xyz->y, &captured_xyz->z);
    float xLo = captured_xyz->x, yLo = captured_xyz->y, zLo = captured_xyz->z; 
    float xHi = xLo, yHi = yLo, zHi = zLo;
    for (vno = 1; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        captured_xyz = &rt->captured_VERTEX_xyz[vno];
        getXYZ(vertex, &captured_xyz->x, &captured_xyz->y, &captured_xyz->z);
        float x = captured_xyz->x, y = captured_xyz->y, z = captured_xyz->z; 
        xLo = MIN(xLo, x); yLo = MIN(yLo, y); zLo = MIN(zLo, z); 
        xHi = MAX(xHi, x); yHi = MAX(yHi, y); zHi = MAX(zHi, z); 
    }
    
    // Initialise the root node, and make it the recentNode
    //
    // Since contains is xLo <= x < xHi etc. the bounds need to be slightly wider than Hi
    // so that the Hi is in a Node
    //
    xHi = widenHi(xHi);
    yHi = widenHi(yHi);
    zHi = widenHi(zHi);

    RealmTreeNode* recentNode  = &rt->root;
    recentNode->xLo = xLo; recentNode->yLo = yLo; recentNode->zLo = zLo;
    recentNode->xHi = xHi; recentNode->yHi = yHi; recentNode->zHi = zHi;

    // Place all the vertices into nodes.  recentNode tries to speed up by assuming some locality.
    // 
    for (vno = 0; vno < mris->nvertices; vno++) {
        recentNode = insertVnoNear(rt, recentNode, vno);
    }

    // Place all the faces into nodes
    //
    int fno;
    for (fno = 0; fno < mris->nfaces; fno++) {
        insertFnoIntoRealmTreeNode(rt, chooseRealmTreeNodeForFno(mris, rt, fno), fno);
    }
        
    if (0) {
        fprintf(stdout,"%s:%d summarizeRealmTree after made\n", __FILE__, __LINE__);
        summarizeRealmTree(rt);
    }
    
    return rt;
}

void noteIfXYZChangedRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ, int vno) {
    VERTEX const * vertex = &mris->vertices[vno];

    chkBnd(0, vno, realmTree->saved_nvertices);

    float x,y,z;
    getXYZ(vertex, &x, &y, &z);
        // Get the new values
        
    Captured_VERTEX_xyz* c = &realmTree->captured_VERTEX_xyz[vno];
        // Get the old values

    if (x == c->x && y == c->y && z == c->z) {
        // ignore if has not moved
        // fprintf(stdout,"noteIfXYZChangedRealmTree vno:%d has not moved\n", vno);       // this happens a lot
        return;
    }
    
    if (x < realmTree->root.xLo || realmTree->root.xHi <= x || 
        y < realmTree->root.yLo || realmTree->root.yHi <= y ||
        z < realmTree->root.zLo || realmTree->root.zHi <= z 
    ) {
        // fprintf(stderr,"noteIfXYZChangedRealmTree vno:%d has ", vno);
        // fprintf(stderr,"moved outside root\n");                                 
            // this almost never happens
            // when it does, must widen the sides that are too tight
        
        widenSubtree(realmTree, x, widenHi(x), y, widenHi(y), z, widenHi(z)); 
        
    } else {
        // fprintf(stderr,"stayed inside root\n");                                  // this happens a few tens of times
    }
    
    // rather than updating now, batch them 
    //  so faces only get moved once when several of their vertexs move
    //
    if (realmTree->nextVnoToUpdatePlus1[vno] == 0) {                                // make sure not already in list
        realmTree->nextVnoToUpdatePlus1[vno] = realmTree->firstVnoToUpdatePlus1;    // append existing list to this node
        realmTree->firstVnoToUpdatePlus1 = vno + 1;                                 // make this vno the start of the list
    }
}

static int addFnoFaceSet(int firstFnoToUpdatePlus1, RealmTree* realmTree, int fno) {
    if (realmTree->nextFnoToUpdatePlus1[fno] == 0) {
        realmTree->nextFnoToUpdatePlus1[fno] = firstFnoToUpdatePlus1;   // add to list
        firstFnoToUpdatePlus1 = fno + 1;
    }
    return firstFnoToUpdatePlus1;
}

static int addFacesToFaceSet(int firstFnoToUpdatePlus1, RealmTree* realmTree, VERTEX const * const vertex) {
    int fi; 
    for (fi = 0; fi < vertex->num; fi++) {
        firstFnoToUpdatePlus1 = addFnoFaceSet(firstFnoToUpdatePlus1, realmTree, vertex->f[fi]);
    }
    return firstFnoToUpdatePlus1;
}

void updateRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    
    int previous_saved_nvertices = realmTree->saved_nvertices;
    int previous_saved_nfaces    = realmTree->saved_nfaces;
    
    resizeRealmTree(realmTree, mris);

    if (previous_saved_nfaces <= interestingFno && interestingFno < realmTree->saved_nfaces) {
        fprintf(stdout, "%s:%d interestingFno:%d should be added soon\n", __FILE__, __LINE__, interestingFno);
    }

    // Pending faces list
    //    
    int firstFnoToUpdatePlus1 = -1;

    // add the new vertices.  recentNode tries to speed up by assuming some locality.
    // 
    RealmTreeNode* recentNode = &realmTree->root;
    int vno;
    for (vno = previous_saved_nvertices; vno < mris->nvertices; vno++) {
        recentNode = insertVnoNear(realmTree, recentNode, vno);
        firstFnoToUpdatePlus1 = addFacesToFaceSet(firstFnoToUpdatePlus1, realmTree, &mris->vertices[vno]);
    }
    
    // add the new faces
    // 
    int fno;
    for (fno = previous_saved_nfaces; fno < mris->nfaces; fno++) {
        firstFnoToUpdatePlus1 = addFnoFaceSet(firstFnoToUpdatePlus1, realmTree, fno);
    }
    
    // process all the pending vno and build the pending face list
    //
    vno = realmTree->firstVnoToUpdatePlus1 - 1;
    while (vno >= 0) {
        VERTEX const * const vertex = &mris->vertices[chkBnd(0,vno,mris->nvertices)];

        // Get its new xyz
        //    
        Captured_VERTEX_xyz* captured_xyz = &realmTree->captured_VERTEX_xyz[vno];
        getXYZ(vertex, &captured_xyz->x, &captured_xyz->y, &captured_xyz->z);

        // Has it changed nodes?
        //
        RealmTreeNode*                n = realmTree->vnoToRealmTreeNode[vno];
        RealmTreeNode const * shouldBeN = deepestContainingNode(n, captured_xyz->x, captured_xyz->y, captured_xyz->z);
        
        if (n != shouldBeN) {
        
            // add the faces to the face set
            //
            firstFnoToUpdatePlus1 = addFacesToFaceSet(firstFnoToUpdatePlus1, realmTree, vertex);

            // remove it from the existing realmTreeNode
            //
            removeVnoFromRealmTree(realmTree, vno);

            // insert it
            //
            insertVnoNear(realmTree, n, vno);
        }
                
        // next
        //
        int next_vno = realmTree->nextVnoToUpdatePlus1[vno] - 1;        // get the next vno, if any
        realmTree->nextVnoToUpdatePlus1[vno] = 0;                       // clear the link for next time
        vno = next_vno;
    }
    realmTree->firstVnoToUpdatePlus1 = -1;				// mark the list as empty

    // process all the pending fno, resetting the links to zero
    //    
    fno = firstFnoToUpdatePlus1 - 1;
    while (fno >= 0) {

        RealmTreeNode * chosenForFno  = chooseRealmTreeNodeForFno(mris, realmTree, fno);
        RealmTreeNode * currentForFno = realmTree->fnoToRealmTreeNode[fno];
        if (chosenForFno != currentForFno) {
            if (currentForFno) removeFnoFromRealmTree(realmTree, fno);
            insertFnoIntoRealmTreeNode(realmTree, chosenForFno, fno);
        }
                
        int next_fno = realmTree->nextFnoToUpdatePlus1[fno] - 1;        // get the next fno, if any
        realmTree->nextFnoToUpdatePlus1[fno] = 0;                       // clear the link for next time
        fno = next_fno;
    }
}

int checkRealmTree(RealmTree const * realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    int count = countXYZChanges(realmTree, mris, getXYZ);
    if (count > 0) {
        fprintf(stdout, "%s:%d mris %d vertex xyz have changed\n", __FILE__, __LINE__, count);
    }
    return count;
}

void getRealmTreeBnds(
    RealmTree* realmTree, float* xLo, float* xHi, float* yLo, float* yHi, float* zLo, float* zHi) {
    *xLo = realmTree->root.xLo;
    *yLo = realmTree->root.yLo;
    *zLo = realmTree->root.zLo;
    *xHi = realmTree->root.xHi;
    *yHi = realmTree->root.yHi;
    *zHi = realmTree->root.zHi;
}


// Realm construction and destruction
//
struct Realm {
    RealmTree     const * realmTree;
    RealmTreeNode const * deepestContainingNode;
    float                 xLo, xHi, yLo, yHi, zLo, zHi;
};

void freeRealm(Realm** realmPtr) {
    freeAndNULL(*realmPtr);
    *realmPtr = NULL;
}

Realm* makeRealm(
    RealmTree const * realmTree, 
    float xLo, float xHi, 
    float yLo, float yHi,
    float zLo, float zHi) {
    //
    // Creates a realm that can be used to quickly find vertices and faces
    // that MIGHT intersect this brick.
    //
    // There are no vertices ON the high wall or outside any of the walls
    // so a realm - which represents the vertices and faces and not some unoccupied space
    // - can be shrunk to these bounds, and needs to be if the deepestContainingNode
    // is to be set correctly, otherwise a realm whose inhabitants might be but aren't outside
    // those bounds would not have any containingNode!

    Realm* r = (Realm*)calloc(1, sizeof(Realm));
    r->realmTree = realmTree;
    r->xLo = xLo = MAX(realmTree->root.xLo,xLo); 
    r->yLo = yLo = MAX(realmTree->root.yLo,yLo); 
    r->zLo = zLo = MAX(realmTree->root.zLo,zLo);
    r->xHi = xHi = MIN(realmTree->root.xHi,xHi); 
    r->yHi = yHi = MIN(realmTree->root.yHi,yHi); 
    r->zHi = zHi = MIN(realmTree->root.zHi,zHi);
    
    // An empty realm contains no points hence can not be intersected by a face
    //
    if (xHi <= xLo || yHi <= yLo || zHi <= zLo) return r;     // r->deepestContainingNode left NULL
    
    // == on the high bound is permissible, since points on that border are not in the realm
    //
    RealmTreeNode const * n = deepestContainingNode(&realmTree->root, xLo, yLo, zLo);
    while (xHi > n->xHi || yHi > n->yHi || zHi > n->zHi) n = n->parent;
    r->deepestContainingNode = n;

    return r;
}


// Quick tests
//
static bool nodeIntersectsRealm(RealmTreeNode const * c, Realm* r) {
    if (c->xHi <= r->xLo || r->xHi <= c->xLo) return false;
    if (c->yHi <= r->yLo || r->yHi <= c->yLo) return false;
    if (c->zHi <= r->zLo || r->zHi <= c->zLo) return false;
    return true;
}

// Iterators
// first call with realmIterator 0
// successive calls return some next fno or vno, may not be ascending order!
// updates realmIterator to some private non-zero value
// returns -1 when no more found
// further calls will cause an error exit
// 
static void moveToNext(RealmIterator* realmIterator, Realm* realm) {
    
    RealmTreeNode const * n = (RealmTreeNode const *)realmIterator->p;
    unsigned long         i = realmIterator->i;

    // More in same leaf?
    //
    {
#ifdef REALM_UNIT_TEST
        if (!n->vnos) *(int*)-1 = 0;   // must be a leaf
#endif
        unsigned long c = i & leafIndexMask;
        c++;
        if (c < n->vnosSize) {
            realmIterator->i++;
            return;
        }
    }
        
    // This leaf is consumed, so search for the next non-empty leaf

GoUp:;
    // Go up from leaf or non-leaf until there is an unprocessed child
    //
    unsigned long c;
    do {
        i >>= ((n->vnos) ? leafIndexBits : childIndexBits);
        n = n->parent;
        c = (i & childIndexMask) + 1;
        if (i == 1) {               
            // exited the realm->deepestContainingNode
            n = NULL;
            goto Done;
        }
    } while (c == childrenSize);
    i += 1;     // Note: c adjusted above

GoDown:;
    // Find the first unprocessed child that might contribute more vno
    //
    RealmTreeNode const * child;
    for (;;) {
        child = n->children[c];
        if (!child->vnos || child->vnosSize)                // the child is a non-leaf or a leaf with children
            if (nodeIntersectsRealm(child, realm))          // and it might contribute to this realm
                break;                                      // it is what we are looking for
        c++;
        if (c >= childrenSize) {                            // no more siblings
            goto GoUp;
        }  
        i++;
    }

    // Go into the child
    //
    n = child;
    i <<= ((n->vnos) ? leafIndexBits : childIndexBits);
    c = 0;

    // If not a leaf, keep going down
    //
    if (!n->vnos) goto GoDown;

Done:;
    // Note for next time
    //
    realmIterator->i = i;
    realmIterator->p = (void*)n;
}

void initRealmIterator(RealmIterator* realmIterator, Realm* realm) {
    // The iterator implements a depth-first walk of the tree
    // It has to keep track of how far below the deepestContainingNode it is, and which of the children it is processing
    //
    // When i becomes 1, it is because we have popped the realm->deepestContainingNode
    //
    RealmTreeNode const * n = realm->deepestContainingNode;
    if (!n) {
        realmIterator->i = 1;
        realmIterator->p = NULL;
        return;
    }
    
    // use 1 to mark the underflow
    //
    unsigned long i = 1 << ((n->vnos) ? leafIndexBits : childIndexBits);
              
    // Down to the deepest leftmost descendent
    //
    while (!n->vnos) {
        n   = n->children[0];
        i <<= ((n->vnos) ? leafIndexBits : childIndexBits);
    }

    // Set up to access the first of these vno's, if any
    //    
    realmIterator->i = i;
    realmIterator->p = (void*)n;
    
    // If none there, pretend already returned
    //
    if (n->vnosSize == 0 || !nodeIntersectsRealm(n, realm)) moveToNext(realmIterator, realm);
}

int realmNextMightTouchFno(Realm* realm, RealmIterator* realmIterator) {
    fprintf(stderr, "%s:%d NYI\n", __FILE__, __LINE__);
    exit(1);
    return 0;
}

int realmNextMightTouchVno(Realm* realm, RealmIterator* realmIterator) {

    RealmTreeNode const * n = (RealmTreeNode const *)realmIterator->p;

    // If there isn't any more
    //
    if (!n) {
        return -1;
    }

    // Get this one
    //
    unsigned long i = realmIterator->i;
    unsigned long c = i & leafIndexMask;
    int const vno = n->vnos[c];

    // Step to the next one
    //
    moveToNext(realmIterator, realm);

    // Return this one
    //
    return vno;
}

static int numberOffnosHereAndDeeper(RealmTreeNode const* n, Realm* realm) {
    if (!nodeIntersectsRealm(n, realm)) return 0;
    int count = n->nFaces;
    if (!n->vnos) {
        int c;
        for (c = 0; c < childrenSize; c++) count += numberOffnosHereAndDeeper(n->children[c], realm);
    }
    return count; 
}

int realmNumberOfMightTouchFno(Realm* realm) {
    RealmTreeNode const* n = realm->deepestContainingNode;
    if (!n) return 0;
    int count = numberOffnosHereAndDeeper(n, realm);
    while ((n = n->parent)) count += n->nFaces;
    return count;  
}

static int fnosHere(RealmTree const* rt, RealmTreeNode const* n, int* fnos, int fnosCapacity, int fnosSize) {
    int fno = n->firstFnoPlus1 - 1;
    while (fno >= 0) {
#ifdef REALM_UNIT_TEST
        if (fnosCapacity <= fnosSize) {
            *(int*)-1 = 0;
        }
#endif
        fnos[fnosSize++] = fno;
        fno = rt->nextFnoPlus1[fno] - 1;
    }
    return fnosSize;
}

static int fnosHereAndDeeper(RealmTree const* rt, Realm* realm, RealmTreeNode const* n, int* fnos, int fnosCapacity, int fnosSize) {
    if (!nodeIntersectsRealm(n, realm)) return fnosSize;
    fnosSize = fnosHere(rt, n, fnos, fnosCapacity, fnosSize);
    if (!n->vnos) {
        int c;
        for (c = 0; c < childrenSize; c++) 
        fnosSize = fnosHereAndDeeper(rt, realm, n->children[c], fnos, fnosCapacity, fnosSize);
    }
    return fnosSize; 
}

int realmMightTouchFno(Realm* realm, int* fnos, int fnosCapacity) {
    RealmTreeNode const* n = realm->deepestContainingNode;
    if (!n) return 0;
    int written = fnosHereAndDeeper(realm->realmTree, realm, n, fnos, fnosCapacity, 0);
    while ((n = n->parent)) written = fnosHere(realm->realmTree, n, fnos, fnosCapacity, written);
    return written;
}

static void summarizeRealmTreeNodeIndent(RealmTreeNode const * n) {
    int i; 
    for (i = 0; i < n->depth; i++) fprintf(stdout,"   |");
}

static void summarizeRealmTreeNode(RealmTree const * realmTree, RealmTreeNode const * n) {
    if (!n) {
        fprintf(stdout,"not in the tree\n");
        return;
    }
    summarizeRealmTreeNodeIndent(n);
    fprintf(stdout,"x:%f..%f y:%f..%f z:%f..:%f nFaces:%d\n", n->xLo, n->xHi, n->yLo, n->yHi, n->zLo, n->zHi, n->nFaces);
    if (n->vnos) {
        summarizeRealmTreeNodeIndent(n);
        fprintf(stdout," vnosSize:%d vno:",n->vnosSize);
        int vi;
        for (vi = 0; vi < n->vnosSize; vi++) {
            fprintf(stdout," %d",n->vnos[vi]);
        }
        fprintf(stdout,"\n");
    }
    if (n->nFaces) {
        summarizeRealmTreeNodeIndent(n);
        fprintf(stdout," nFaces:%d fno:",n->nFaces);
        int entryFno = n->firstFnoPlus1 - 1;
        while (entryFno >= 0) {
            fprintf(stdout," %d",entryFno);
            entryFno = realmTree->nextFnoPlus1[entryFno] - 1;
        }
        fprintf(stdout,"\n");
    }
}

static int summarizeRealmTreeSubtree(RealmTree const * realmTree, RealmTreeNode const * n, int targetDepth) {
    int hasUnreachedChildren = 0;
    int atDepth = (n->depth == targetDepth);
    
    if (atDepth) summarizeRealmTreeNode(realmTree, n);
    if (n->vnos) {
        ;
    } else if (n->depth < targetDepth) {
        int c;
        for (c = 0; c < childrenSize; c++) 
            hasUnreachedChildren |= summarizeRealmTreeSubtree(realmTree, n->children[c], targetDepth);
    } else {
        hasUnreachedChildren=1;
    }
    return hasUnreachedChildren;
}

void summarizeRealmTree(RealmTree const * realmTree) {
    int depth;
    for (depth = 0; ; depth++)
        if (!summarizeRealmTreeSubtree(realmTree, &realmTree->root, depth)) break;
}

void summarizeRealmTreeVno(RealmTree const * realmTree, int vno) {
    fprintf(stdout,"vno:%d\n",vno);
    RealmTreeNode* n = realmTree->vnoToRealmTreeNode[vno];
    summarizeRealmTreeNode(realmTree, n);
}

void summarizeRealmTreeFno(RealmTree const * realmTree, int fno) {
    fprintf(stdout,"summarizeRealmTreeFno fno:%d\n",fno);
    RealmTreeNode* n = realmTree->fnoToRealmTreeNode[fno];
    summarizeRealmTreeNode(realmTree, n);
    FACE const * face = &realmTree->mris->faces[fno];
    fprintf(stdout,"vno");
    int vi;
    for (vi = 0; vi < VERTICES_PER_FACE; vi++) {
        fprintf(stdout," %d",face->v[vi]);
    }
    fprintf(stdout,"\n");
    for (vi = 0; vi < VERTICES_PER_FACE; vi++) {
        summarizeRealmTreeVno(realmTree, face->v[vi]); 
    }
}
static void sortIntSoFirstLo(int* a, int* b) {
    if (*a > *b) { int temp = *a; *a = *b; *b = temp; }
}

static void sortFloatSoFirstLo(float* a, float* b) {
    if (*a > *b) { float temp = *a; *a = *b; *b = temp; }
}




// Wrapper for realloc that 
//      adjusts the capacity
//      zero's the extension
//
static void growCapacity(int* p_capacity, int minCapacity) {
    int old_capacity = *p_capacity;
    int new_capacity = old_capacity + old_capacity/3;
    if (new_capacity < minCapacity) new_capacity = minCapacity;
    *p_capacity = new_capacity;
}

static void growInts(int** p_old, int old_capacity, int new_capacity) {
    int* old = *p_old;
    int* new = (int*)realloc(old, new_capacity*sizeof(int));
    bzero(&new[old_capacity], (new_capacity - old_capacity)*sizeof(int));
    *p_old = new;
}

static void growCapacityAndInts(int** p_old, int* p_capacity, int minCapacity) {
    int  old_capacity = *p_capacity;
    growCapacity(p_capacity, minCapacity);
    growInts(p_old, old_capacity, *p_capacity);
}

// The cells of the projection plane that contain the indexs for the great arcs that might intersect the cell
//
#define GRID_WIDTH  16
#define GRID_HEIGHT 16
#define CELLS_SIZE (GRID_WIDTH*GRID_HEIGHT + 1) // Need a grid, plus one cell for everything that has at least partial exceeds the grid
#define UNIVERSAL_CELL (CELLS_SIZE - 1)

typedef struct Cell {
    int  size, capacity;
    int* indexs;
} Cell;

static void finiCell(Cell* cell) {
    freeAndNULL(cell->indexs);
}

typedef struct Pair {                           // an entry in a chain off the GreatArcSet pairHeadsPlus1 hash table
    int nextPlus1;                              // index plus 1 (0 is end of list) of the next item
    int loVno;
    int hiVno;
} Pair;


#define PAIRS_HEADS_SIZE (1<<10)
#define PAIRS_HEADS_MASK (PAIRS_HEADS_SIZE-1)

struct GreatArcSet {
    MRIS*   mris;
    
    int     capacity, size,                     // the list of great arcs
            size_dispersed;                     // the size of the tail of the great arcs not yet in the cells
    int*    keys;                               // mrisurf gives a key to each great arc, to id the arc in the callback
    int*    loVnos;                             // beginning vno of the great arc
    int*    hiVnos;                             // ending vno of the great arc
    
    int*    passed;                             // used during one callback to note those passed to de-duplicate the multiple cells
    
    int     pairHeadsPlus1[PAIRS_HEADS_SIZE];   // a hash table, indexs into pairs
    int     pairsSize, pairsCapacity;           // chains off the hash table
    Pair*   pairs;                              //      the entries in the chains
    
    float   ax,ay,az,bx,by,bz,cx,cy,cz;         // used to project the vno into the rectangle h,w
    float   minH,minW,scaleH,scaleW;            // used to project the h,w into the cel grid

    Cell    cells[CELLS_SIZE];
};

static Cell* getCell(GreatArcSet* set, int w, int h) {
    // Does not include the UNIVERSAL_CELL
    return &set->cells[chkBnd(0,w,GRID_WIDTH)*GRID_HEIGHT + chkBnd(0,h,GRID_HEIGHT)];
} 


static void greatArcSet_project(GreatArcSet* set, float x, float y, float z, float* w, float* h, bool* universal_cell, bool trace) {

    // Project into the coordinate space invented for the defect
    //
    float
        px = set->cx*x + set->cy*y + set->cz*z,
        py = set->ax*x + set->ay*y + set->az*z,
        pz = set->bx*x + set->by*y + set->bz*z;

    if (trace) {
        fprintf(stdout, "%s:%d x:%g y:%g z:%g \n", __FILE__, __LINE__, x,y,z);
        fprintf(stdout, "  ax:%g ay:%g az:%g \n", set->ax,set->ay,set->az);
        fprintf(stdout, "  bx:%g by:%g bz:%g \n", set->bx,set->by,set->bz);
        fprintf(stdout, "  cx:%g cy:%g cz:%g \n", set->cx,set->cy,set->cz);
    }
    
    // Project into the perpendicular plane

    if (px <= 0.1) {                // This could be zero, but making it 0.1 reduce the range of the w h
        *w = -1.0f;                 // because our cx,... data is typically approx 100.0
        *h = -1.0f;
        *universal_cell = true;
        
        if (1) {
            static int once;
            if (!once) { once = 1;
                fprintf(stdout, "Large defect found - not a problem, but a curiousity\n");
            }
        }
    }
    
    *universal_cell = false;
    *w = pz/px;
    *h = py/px;

    if (trace) {
        fprintf(stdout, "  w:%g h:%g\n", *w, *h);
    }
}


static void greatArcSet_getCellCoords(GreatArcSet* set, float x, float y, float z, int* p_wI, int* p_hI, bool* universal_cell, bool trace) {

    float h,w;
    greatArcSet_project(set, x, y, z, &w, &h, universal_cell, trace);
    
    int wI = (int)((w - set->minW)*set->scaleW);
    int hI = (int)((h - set->minH)*set->scaleH);
   
    if (trace) {
        fprintf(stdout, "%s:%d wI:%d = (w:%g - minW:%g)*scaleW:%g\n",__FILE__,__LINE__,wI,w,set->minW,set->scaleW);
        fprintf(stdout, "%s:%d hI:%d = (w:%g - minH:%g)*scaleH:%g\n",__FILE__,__LINE__,hI,h,set->minH,set->scaleH);
    }
    
    if (*universal_cell
    || hI < 0 || GRID_HEIGHT <= hI
    || wI < 0 || GRID_WIDTH  <= wI) {
        if (0) {
            static int aFewTimes;
            if (aFewTimes < 100) { aFewTimes++;
                fprintf(stdout, "Large defect found - a curiousity, not a problem - *universal_cell:%d\n", *universal_cell);
            }
        }
        *universal_cell = true;
        *p_wI = -1;
        *p_hI = -1;
    } else {  
        *p_wI = wI;
        *p_hI = hI;
    }
}


static void greatArcSet_cellInsert(GreatArcSet* set, Cell* cell, int index) {
    if (cell->size == cell->capacity) growCapacityAndInts(&cell->indexs, &cell->capacity, 16);
    cell->indexs[chkBnd(0,cell->size++,cell->capacity)] = index;
}

void freeGreatArcSet(GreatArcSet** setPtr) {
    GreatArcSet* set = *setPtr; *setPtr = NULL;
    if (!set) return;
    { int i; for (i=0; i<CELLS_SIZE; i++) finiCell(&set->cells[i]); }
    freeAndNULL(set->pairs);
    freeAndNULL(set->passed);
    freeAndNULL(set->hiVnos);
    freeAndNULL(set->loVnos);
    freeAndNULL(set->keys);
}

GreatArcSet* makeGreatArcSet(MRIS* mris) {
    GreatArcSet* set = (GreatArcSet*)calloc(1, sizeof(GreatArcSet));
    set->mris = mris;
    return set;
}

static void growGreatArcBuffer(GreatArcSet* set) {
    size_t old_capacity = set->capacity;
    growCapacity(&set->capacity, 128);
    growInts(&set->keys,   old_capacity, set->capacity);
    growInts(&set->loVnos, old_capacity, set->capacity);
    growInts(&set->hiVnos, old_capacity, set->capacity);
    growInts(&set->passed, old_capacity, set->capacity);
}


static bool isInPairs(int* p_pairHeadsIndex, GreatArcSet* set, int vno0, int vno1) {
    int pairHeadsIndex = *p_pairHeadsIndex = chkBnd(0,((vno0*12797)^(vno1*97127)) & PAIRS_HEADS_MASK, PAIRS_HEADS_SIZE);
    int next = set->pairHeadsPlus1[pairHeadsIndex] - 1;
    while (next >= 0) {
        Pair* pair = &set->pairs[chkBnd(0,next,set->pairsSize)];
        if (pair->loVno == vno0 && pair->hiVno == vno1) 
            return true;
        next = pair->nextPlus1 - 1;
    }
    return false;
}


void insertGreatArc(GreatArcSet* set, int key, int vno0, int vno1) {

    // Order the vno's to help detect duplicate edges
    //
    sortIntSoFirstLo(&vno0, &vno1);

    // Check for duplicates, and ignore if present
    //
    int pairHeadsIndex;
    if (isInPairs(&pairHeadsIndex, set, vno0, vno1)) return;                            // duplicate, ignored
    
    // Insert
    //
    if (set->pairsSize == set->pairsCapacity) {
        growCapacity(&set->pairsCapacity,128);
        set->pairs = (Pair*)realloc(set->pairs, set->pairsCapacity * sizeof(Pair));      // note: not zeroed
    }

    Pair* pair = &set->pairs[set->pairsSize++];
    pair->nextPlus1 = set->pairHeadsPlus1[pairHeadsIndex];
    pair->loVno     = vno0;
    pair->hiVno     = vno1;
    
    set->pairHeadsPlus1[pairHeadsIndex] = set->pairsSize;                                  // +1 has been done by the set->pairsSize++

    if (set->size == set->capacity) growGreatArcBuffer(set);

    int i = set->size++;
    set->keys  [i] = key;
    set->loVnos[i] = vno0;
    set->hiVnos[i] = vno1;
}


static void decideProjection(GreatArcSet* set) {
    
    // Sum the coordinates of all (or a subset) of the vertexs
    //
    float sumX = 0, sumY = 0, sumZ = 0;
    {   int index;
        for (index = 0; index < set->size; index++) {
            VERTEX const * v0 = &set->mris->vertices[chkBnd(0,set->loVnos[index],set->mris->nvertices)];
            VERTEX const * v1 = &set->mris->vertices[chkBnd(0,set->hiVnos[index],set->mris->nvertices)];
            sumX += v0->cx + v1->cx;
            sumY += v0->cy + v1->cy;
            sumZ += v0->cz + v1->cz;
        }
    }

    // Compute the plane through the points that is going to hold the projection of the defect for search purposes
    // It is close enough to the sphere for small defects that it can guide the search for intersecting lines.
    //    
    // This is the center of the points, aka the vector from the origin to the center, 
    // so it is the normal to the plane that is tangent to the surface where the defect is centered.
    // Make it the unit vector to use as the axis.
    //
    float cx = sumX, cy = sumY, cz = sumZ;
    {
        float clen = sqrt(cx*cx+cy*cy+cz*cz);  if (clen == 0.0) clen = 1.0; clen = 1.0/clen; cx *= clen; cy *= clen; cz *= clen; 
    }

    // Construct some unit vector in the plane
    //    
    float ax,ay,az;                                         // construct a perpendicular to cx,cy,cz [make the dot product == 0.0]
    if (fabs(cx) >= fabs(cz) && fabs(cx) >= fabs(cz))       // cz is smallest
        ax = cy, ay = -cx, az = 0;
    else if (fabs(cx) >= fabs(cy) && fabs(cz) >= fabs(cy))  // cy is smallest
        ax = cz, ay = 0, az = -cx;
    else                                                    // cx is smallest
        ax = 0,  ay = cz, az = -cy;
    {
        float alen = sqrt(ax*ax+ay*ay+az*az);  if (alen == 0.0) alen = 1.0; alen = 1.0/alen; ax *= alen; ay *= alen; az *= alen; 
    }

    // Construct the third axis
    //
    float                                                   // cross product to get the third axis
        bx = ay*cz - az*cy,                                 // it will be a unit vector, but just in case tidy it up..
        by = az*cx - ax*cz,
        bz = ax*cy - ay*cx;
    {
        float blen = sqrt(bx*bx+by*by+bz*bz);  if (blen == 0.0) blen = 1.0; blen = 1.0/blen; bx *= blen; by *= blen; bz *= blen; 
    }

    // Save - used to project the lines later
    //
    set->ax = ax; set->ay = ay; set->az = az;
    set->bx = bx; set->by = by; set->bz = bz;
    set->cx = cx; set->cy = cy; set->cz = cz;

    // Decide the range of the projected data
    //
    float minH = +FLT_MAX, minW = +FLT_MAX;
    float maxH = -FLT_MAX, maxW = -FLT_MAX;
    {   int* vnos = set->loVnos;
        for (;;) {
            int index;
            for (index = 0; index < set->size; index++) {
                VERTEX const * v = &set->mris->vertices[chkBnd(0,vnos[index],set->mris->nvertices)];
                bool  universal_cell;
                float h,w;
                greatArcSet_project(set, v->cx, v->cy, v->cz, &h, &w, &universal_cell, false);
                if (!universal_cell) {
                    minH = MIN(minH,h); maxH = MAX(maxH,h);
                    minW = MIN(minW,w); maxW = MAX(maxW,w);
                }
            }
            if (vnos == set->hiVnos) break;
            vnos = set->hiVnos;
        }
    }

    set->minH = minH; if (maxH == minH) maxH = minH + 1; set->scaleH = GRID_HEIGHT/(maxH-minH);
    set->minW = minW; if (maxW == minW) maxW = minW + 1; set->scaleW = GRID_WIDTH /(maxW-minW);

    if (0) {
        bool universal_cell;
        int hI, wI;
        greatArcSet_getCellCoords(set, sumX/set->size, sumY/set->size, sumZ/set->size, &hI, &wI, &universal_cell, false);  
        fprintf(stderr, "%s:%d center -> (%d,%d) %s\n", __FILE__, __LINE__, hI,wI,universal_cell?"universal_cell":"");
    }
}


static void disperseGreatArcsIntoCells(GreatArcSet* set) 
{
    if (set->size_dispersed == set->size) return;

    if (set->size_dispersed == 0) decideProjection(set);

    // Insert the new arcs into the cells
    //
    int index;
    for (index = set->size_dispersed; index < set->size; index++) {

        int const vno0 = set->loVnos[index];
        int const vno1 = set->hiVnos[index];
        
        bool const trace = false; // vno0 == 496 && vno1 == 151637;
        
        VERTEX const * v0 = &set->mris->vertices[chkBnd(0,vno0,set->mris->nvertices)];
        VERTEX const * v1 = &set->mris->vertices[chkBnd(0,vno1,set->mris->nvertices)];

        bool universal_cell0, universal_cell1;
        int hI0, hI1, wI0, wI1;
        greatArcSet_getCellCoords(set, v0->cx, v0->cy, v0->cz, &wI0, &hI0, &universal_cell0, trace);
        greatArcSet_getCellCoords(set, v1->cx, v1->cy, v1->cz, &wI1, &hI1, &universal_cell1, trace);

        if (trace) {
            fprintf(stdout, "%s:%d v0:%d v1:%d dispersed to cell coords (%d..%d, %d..%d) %s\n", __FILE__, __LINE__,
                vno0,vno1, wI0,wI1,hI0,hI1,  (universal_cell0 || universal_cell1)?"universal_cell":"");
            fprintf(stdout, " v0 at (%g,%g,%g) v1 at (%g,%g,%g)\n", v0->cx,v0->cy,v0->cz, v1->cx,v1->cy,v1->cz);
        }

        if (universal_cell0 || universal_cell1) {
            greatArcSet_cellInsert(set, &set->cells[UNIVERSAL_CELL], index);
            continue;   
        }
        
        // There are many ways to do this, but most edges are within a cell or adjacent cells
        // so it is not worth being clever
        //
        sortIntSoFirstLo(&hI0,&hI1);
        sortIntSoFirstLo(&wI0,&wI1);
        
        int h,w;
        for (h = hI0; h <= hI1; h++)           
        for (w = wI0; w <= wI1; w++)
            greatArcSet_cellInsert(set, getCell(set,w,h), index);           
    }

    set->size_dispersed = set->size; 
}


static bool possiblyIntersectingCell(GreatArcSet* set,
    void* callbackCtx,
    bool (*callback)(
        void* callbackCtx, 
        int   key, 
        bool* isHit),               // returns true if should keep going, false if no more needed
    Cell* cell,
    bool  tracing,
    int*  pFirstPassedPlus1,
    int*  pCallBackCount,
    int*  pCallBackFirstFound) {

    int i;
    for (i = 0; i < cell->size; i++) {
        int index = chkBnd(0,cell->indexs[i],set->size);

        // eliminate duplicates
        //
        if (set->passed[index]) continue;           // in list, hence this is a duplicate
        set->passed[index] = *pFirstPassedPlus1;    // prepend to list
        *pFirstPassedPlus1 = index + 1;


        // NOTE: faster intersect code could go in here


        (*pCallBackCount)++;
        bool isHit = false;
        bool keepGoing = (*callback)(callbackCtx,set->keys[index],&isHit);
        if (isHit) {
            if (!*pCallBackFirstFound) *pCallBackFirstFound = *pCallBackCount;
        }
        if (!keepGoing) {
            return true;
        }
    }
    return false;
}

void possiblyIntersectingGreatArcs(GreatArcSet* set,
    void*        callbackCtx,
    bool (*callback)(
        void* callbackCtx, 
        int   key, 
        bool* isHit),               // returns true if should keep going, false if no more needed
    float x0, float y0, float z0,   // ends of the line, need not be a unit vector
    float x1, float y1, float z1,
    bool  tracing)
{
    tracing = false;

    if (tracing) {
        printf("possiblyIntersectingGreatArcs x:%g..%g y:%g..%g z:%g..%g\n", 
            x0,x1,y0,y1,z0,z1);
    }

    // Once all the arcs are available, they can be placed
    //
    disperseGreatArcsIntoCells(set);

    // Once placed, they can be looked for in the cells...
    //
    bool universal_cell0, universal_cell1;
    int wI0, wI1, hI0, hI1;
    greatArcSet_getCellCoords(set, x0, y0, z0, &wI0, &hI0, &universal_cell0, false);
    greatArcSet_getCellCoords(set, x1, y1, z1, &wI1, &hI1, &universal_cell1, false);

    if (universal_cell0 || universal_cell1) {
        wI0 = 0; wI1 = GRID_WIDTH -1;
        hI0 = 0; hI1 = GRID_HEIGHT-1;
    } else {
        sortIntSoFirstLo(&wI0,&wI1);
        sortIntSoFirstLo(&hI0,&hI1);          
    }
    
    // There are lots of ways to do this, but most edges are within a cell or adjacent cells
    // so it is not worth being clever, however if there is more than one cell, there is a risk of duplicates
    // 
    //
    int  callBackCount = 0, callBackFirstFound = -1;

    int firstPassedPlus1 = -1;  // 0 means not passed, -1 means end of list
    
    int w,h;
    for (w = wI0; w <= wI1; w++) 
    for (h = hI0; h <= hI1; h++) {
        Cell* cell = getCell(set,w,h);           
        if (possiblyIntersectingCell(set, callbackCtx, callback, cell, 
                tracing, &firstPassedPlus1, &callBackCount, &callBackFirstFound)) goto Found;
    }

    if (possiblyIntersectingCell(set, callbackCtx, callback, &set->cells[UNIVERSAL_CELL], 
                tracing, &firstPassedPlus1, &callBackCount, &callBackFirstFound)) goto Found;
        
Found:
    // Reset the passed list
    //
    {   int index;
        while ((index = firstPassedPlus1 - 1) >= 0) {
            firstPassedPlus1 = set->passed[index];
            set->passed[index] = 0;
        }
    }
    
    // Done
    //
    if (tracing) {
        fprintf(stderr, 
            " callBackCount:%d, callBackFirstFound:%d total size:%d\n", 
            callBackCount, callBackFirstFound, set->size);
    }
}


void possiblyIntersectingGreatArcs_Debug(                           // show how vno0..vno1 interacts with the arc
    GreatArcSet* set,
    float x0, float y0, float z0,   // ends of the arc, need not be a unit vector
    float x1, float y1, float z1,
    int vno0, int vno1)
{
    sortIntSoFirstLo(&vno0, &vno1);

    VERTEX const * v0 = &set->mris->vertices[chkBnd(0,vno0,set->mris->nvertices)];
    VERTEX const * v1 = &set->mris->vertices[chkBnd(0,vno1,set->mris->nvertices)];

    fprintf(stdout, "%s:%d possiblyIntersectingGreatArcs_Debug vno0:%d vno1:%d\n", __FILE__, __LINE__, vno0,vno1);
    fprintf(stdout, " v0 at (%g,%g,%g) v1 at (%g,%g,%g)\n", v0->cx,v0->cy,v0->cz, v1->cx,v1->cy,v1->cz);

    int pairHeadsIndex;
    if (isInPairs(&pairHeadsIndex, set, vno0, vno1)) {              // already in somewhere
        fprintf(stderr, "Is in GreatArcSet pairs\n");
    }
    
    int target_index;
    for (target_index = 0; target_index < set->size; target_index++) {
        if (set->loVnos[target_index] == vno0 && set->hiVnos[target_index] == vno1) break;
    }
    if (target_index == set->size) {
        fprintf(stderr, "Not in GreatArcSet list\n");
        return;
    }
    
    bool universal_cell0, universal_cell1;
    int wI0, wI1, hI0, hI1;
    greatArcSet_getCellCoords(set, x0, y0, z0, &wI0, &hI0, &universal_cell0, true);
    greatArcSet_getCellCoords(set, x1, y1, z1, &wI1, &hI1, &universal_cell1, true);

    if (universal_cell0 || universal_cell1) {
        wI0 = 0; wI1 = GRID_WIDTH-1;
        hI0 = 0; hI1 = GRID_HEIGHT-1;
    } else {
        sortIntSoFirstLo(&wI0,&wI1);
        sortIntSoFirstLo(&hI0,&hI1);          
    }
    
    fprintf(stderr, "Should be in (%d..%d,%d..%d) %s.  Found in ", wI0, wI1, hI0, hI1, (universal_cell0||universal_cell1) ? "universal cell":"");
    int w,h;
    for (w = 0; w < GRID_WIDTH;  w++)
    for (h = 0; h < GRID_HEIGHT; h++) {           
        Cell* cell = getCell(set,w,h);           
        int i;
        for (i = 0; i < cell->size; i++) {
            int index = chkBnd(0,cell->indexs[i],set->size);
            if (target_index == index) {
                fprintf(stderr, " (%d,%d)", w,h); 
            }
        }
    }
        Cell* cell = &set->cells[UNIVERSAL_CELL];
        int i;
        for (i = 0; i < cell->size; i++) {
            int index = chkBnd(0,cell->indexs[i],set->size);
            if (target_index == index) {
                fprintf(stderr, " cells[UNIVERSAL_CELL]"); 
            }
        }
    fprintf(stderr, "\n");
}
