#include "realm.h"

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

#include "fnv_hash.h"

#ifdef REALM_UNIT_TEST
    //
    // rm -f a.out ; gcc -o a.out -I ../include realm.c |& less ; ./a.out
    //
    
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
    
    static void bevins_break()
    {
        printf("bevins_break\n");
    }

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

    void test(int nvertices, int useDuplicates) {
        printf("Test nvertices:%d useDuplicates:%d\n", nvertices, useDuplicates);
        
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
                        if (fno >= mris.nfaces) bevins_break();
                        FACE* face = &mris.faces[fno];
                        face->v[0] = vnos[v0];
                        face->v[1] = vnos[v1]; 
                        face->v[2] = vnos[v2];
                        int vi;
                        for (vi = 0; vi < 3; vi++) {
                            int vno = face->v[vi];
                            VERTEX* v = &mris.vertices[vno];
                            if (v->num >= MAX_FACES_PER_VERTEX) bevins_break();
                            v->f[v->num++] = fno;
                        }
                        fno++;
                    }
                }
            }
        }
        mris.nfaces = fno;  // decrease to the correct number

        free(vnos ); vnos  = NULL;               
        free(ctx_x); ctx_x = NULL;
        free(ctx_y); ctx_y = NULL;
        free(ctx_z); ctx_z = NULL;
                
        // Make a tree
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
            printf("Moved vertices, now checking\n");
            checkRealmTree(realmTree, &mris, getSomeXYZ);
            printf("Checked realmTree now updating\n");
            updateRealmTree(realmTree, &mris, getSomeXYZ);
            printf("Updated realmTree, now checking\n");
            checkRealmTree(realmTree, &mris, getSomeXYZ);
            printf("Checked realmTree\n");
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
                printf("fCount:%d x:%f..%f y:%f.%f z:%f..%f\n", fCount, xLo, xHi, yLo, yHi, zLo, zHi);
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
                    printf("counter:%d ri.i:%ld ri.p:%p\n", counter, realmIterator.i, realmIterator.p); 
                    bevins_break();
                }
#endif
                vno = realmNextMightTouchVno(realm, &realmIterator);
#ifdef REALM_UNIT_TEST
                if (vno < -1 || vno >= mris.nvertices) {
                    printf("ERROR, vno:%d is illegal\n", vno); 
                    bevins_break();
                    exit(1);
                }
#endif
                if (0 > vno) break;
                if (counter == 0 || states[vno]) {
                    printf("ERROR, vno:%d reported again when counter:%d, was reported counter:%d\n", vno, counter, states[vno]); 
#ifdef REALM_UNIT_TEST
                    bevins_break();
#endif
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
                   &&  zLo <= v->someZ && v->someZ < zHi) printf("ERROR, vno:%d was not reported\n", vno);
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
                    printf("ERROR, fnos[fnosI]:%d fnos[fnosI+1]:%d\n", fnos[fnosI], fnos[fnosI + 1]);
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
                    printf("ERROR, fno:%d was not reported\n", fno);
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
                    printf("fnosSize:%d mris.nfaces:%d fNoBenefitCount:%d fHasBenefitCount:%d\n", 
                        fnosSize, mris.nfaces, fNoBenefitCount, fHasBenefitCount);
                }
            }
            
            // Done
            //
            free(fnos);
            free(states);
            freeRealm(&realm);
        }
            
        freeRealmTree(&realmTree);
        
        free(mris.vertices);        
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
            free(child);
        }
    } else {
        if (n->vnos != n->vnosBuff) free(n->vnos);
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
        bevins_break();
#endif

    float xMid = n->children[1]->xLo;
    float yMid = n->children[2]->yLo;
    float zMid = n->children[4]->zLo;
    
    int c = ((x < xMid) ? 0 : 1) + ((y < yMid) ? 0 : 2) + ((z < zMid) ? 0 : 4);
    
#ifdef REALM_UNIT_TEST
    if (!nodeContains(n->children[c], x,y,z)) 
        bevins_break();
#endif
    
    return c;
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
    Captured_VERTEX_xyz const * const captured_xyz = &realmTree->captured_VERTEX_xyz[vno];
    float const x = captured_xyz->x, y = captured_xyz->y, z = captured_xyz->z;

#ifdef REALM_UNIT_TEST
    MRIS const* mris = realmTree->mris;
    VERTEX const* v = &mris->vertices[vno];
    if (x != v->someX || y != v->someY || z != v->someZ) 
        bevins_break();
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
            if (n->vnos != n->vnosBuff) free(n->vnos);
            n->vnos = p;
        }
        
        // Can insert 
        //
        if (n->vnosSize < n->vnosCapacity) {
#ifdef REALM_UNIT_TEST    
            if (!nodeContains(n, x,y,z)) 
                bevins_break();
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
            bevins_break();
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
            if (child->depth >= maxDepth) bevins_break();
#endif
            if (c&1) child->xLo = xMid, child->xHi = n->xHi; else child->xLo = n->xLo, child->xHi = xMid; 
            if (c&2) child->yLo = yMid, child->yHi = n->yHi; else child->yLo = n->yLo, child->yHi = yMid;
            if (c&4) child->zLo = zMid, child->zHi = n->zHi; else child->zLo = n->zLo, child->zHi = zMid;
        }
        
        // Insert the saved vno into their child
        for (vi = 0; vi < vnosSize; vi++) {
            Captured_VERTEX_xyz* cxyz = &realmTree->captured_VERTEX_xyz[vnos[vi]];
            insertIntoChild(realmTree, n, vnos[vi]);
        }
    }

    // Insert this vno into the right child
    //
    return insertIntoChild(realmTree, n, vno);
}

static RealmTreeNode* insertVnoNear(RealmTree* realmTree, RealmTreeNode* n, int vno) {
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
    
    free(rt->nextFnoToUpdatePlus1);
    free(rt->nextVnoToUpdatePlus1);
    free(rt->nextFnoPlus1);
    free(rt->fnoToRealmTreeNode);
    free(rt->vnoToRealmTreeNode);
    destroyRealmTreeNode(&rt->root);
    
    free(rt);
}

static float widenHi(float lo, float hi) {
    float step = FLT_MIN;
    while (hi + step == hi) {
        step *= 2.0;
    }
    return hi + step;
}

static RealmTreeNode * chooseRealmTreeNodeForFno(MRIS const * const mris, RealmTree const * const rt, int const fno) {
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


static void insertFnoIntoRealmTreeNode(RealmTree* realmTree, RealmTreeNode* n, int fno) {
    realmTree->fnoToRealmTreeNode[fno] = n;
    realmTree->nextFnoPlus1[fno]       = n->firstFnoPlus1;
    n->firstFnoPlus1                   = fno + 1;
    // adjust the count
    n->nFaces++;
}

static void removeFnoFromRealmTree(RealmTree* realmTree, int fno) {
    RealmTreeNode* n = realmTree->fnoToRealmTreeNode[fno];
    realmTree->fnoToRealmTreeNode[fno] = NULL;

    // find in the list
    int prevFno = -1;                                           // this list is usually very small
    int entryFno = n->firstFnoPlus1 - 1;                        // should this search be a performance problem
    while (entryFno != fno) {                                   //      change to a per-node btree
        prevFno = entryFno;
        entryFno = realmTree->nextFnoPlus1[entryFno] - 1;
    }

    // remove from the list
    if (prevFno < 0) {
        n->firstFnoPlus1 = realmTree->nextFnoPlus1[fno];
    } else {
        realmTree->nextFnoPlus1[prevFno] = realmTree->nextFnoPlus1[fno];
    }

    // be tidy...
    realmTree->nextFnoPlus1[fno] = 0;

    // adjust the count
    n->nFaces--;
}

RealmTree* makeRealmTree(MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    // Fills in the tree using the existing position of 
    // the vertices and faces
    RealmTree* rt = (RealmTree*)calloc(1, sizeof(RealmTree));
    constructRealmTreeNode(&rt->root, NULL);
    rt->mris                = mris;
    rt->captured_VERTEX_xyz = (Captured_VERTEX_xyz*)calloc(mris->nvertices, sizeof(Captured_VERTEX_xyz));
    rt->vnoToRealmTreeNode  = (RealmTreeNode**     )calloc(mris->nvertices, sizeof(RealmTreeNode*));
    rt->fnoToRealmTreeNode  = (RealmTreeNode**     )calloc(mris->nfaces,    sizeof(RealmTreeNode*));

    rt->nextFnoPlus1        = (int*                )calloc(mris->nfaces,    sizeof(int           ));

    rt->firstVnoToUpdatePlus1 = -1; // end of list marker, since none pending
    rt->nextVnoToUpdatePlus1  = (int*              )calloc(mris->nvertices, sizeof(int           ));

    rt->nextFnoToUpdatePlus1  = (int*              )calloc(mris->nfaces,    sizeof(int           ));

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
    xHi = widenHi(xLo,xHi);
    yHi = widenHi(yLo,yHi);
    zHi = widenHi(zLo,zHi);

    RealmTreeNode* recentNode  = &rt->root;
    recentNode->xLo = xLo; recentNode->yLo = yLo; recentNode->zLo = zLo;
    recentNode->xHi = xHi; recentNode->yHi = yHi; recentNode->zHi = zHi;

    // Place all the vertices into nodes.  recentNode tries to speed up by assuming some locality.
    // 
    for (vno = 0; vno < mris->nvertices; vno++) {
        captured_xyz = &rt->captured_VERTEX_xyz[vno];
        recentNode = insertVnoNear(rt, recentNode, vno);
    }

    // Place all the faces into nodes
    //
    int fno;
    for (fno = 0; fno < mris->nfaces; fno++) {
        insertFnoIntoRealmTreeNode(rt, chooseRealmTreeNodeForFno(mris, rt, fno), fno);
    }
        
    if (0) {
        printf("%s:%d summarizeRealmTree after made\n", __FILE__, __LINE__);
        summarizeRealmTree(rt);
    }
    
    return rt;
}

void noteIfXYZChangedRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ, int vno) {
    VERTEX const * vertex = &mris->vertices[vno];

    float x,y,z;
    getXYZ(vertex, &x, &y, &z);
        // Get the new values
        
    Captured_VERTEX_xyz* c = &realmTree->captured_VERTEX_xyz[vno];
        // Get the old values

    if (x == c->x && y == c->y && z == c->z) {
        // ignore if has not moved
        // printf("noteIfXYZChangedRealmTree vno:%d has not moved\n", vno);       // this happens a lot
        return;
    }
    
    if (x < realmTree->root.xLo || realmTree->root.xHi <= x || 
        y < realmTree->root.yLo || realmTree->root.yHi <= y ||
        z < realmTree->root.zLo || realmTree->root.zHi <= z 
    ) {
        printf("noteIfXYZChangedRealmTree vno:%d has ", vno);
        printf("moved outside root\n");                                 // this almost never happens
    } else {
        // printf("stayed inside root\n");                              // this happens a few tens of times
    }
    
    // rather than updating now, batch them 
    //  so faces only get moved once when several of their vertexs move
    //
    if (realmTree->nextVnoToUpdatePlus1[vno] == 0) {                                // make sure not already in list
        realmTree->nextVnoToUpdatePlus1[vno] = realmTree->firstVnoToUpdatePlus1;    // append existing list to this node
        realmTree->firstVnoToUpdatePlus1 = vno + 1;                                 // make this vno the start of the list
    }
}

void updateRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    
    // process all the pending vno and build the pending face list
    //
    int firstFnoToUpdatePlus1 = -1;
    int vno = realmTree->firstVnoToUpdatePlus1 - 1;
    while (vno >= 0) {
        VERTEX const * const vertex = &mris->vertices[vno];

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
            { 
                int fi; 
                for (fi = 0; fi < vertex->num; fi++) {
                    int fno = vertex->f[fi];
                    if (realmTree->nextFnoToUpdatePlus1[fno] != 0) continue;        // already in list
                    realmTree->nextFnoToUpdatePlus1[fno] = firstFnoToUpdatePlus1;   // add to list
                    firstFnoToUpdatePlus1 = fno + 1;
                }
            }

            // remove it from the existing realmTreeNode
            //
            realmTree->vnoToRealmTreeNode[vno] = NULL;
            {
                int vi;
                for (vi = n->vnosSize - 1; n->vnos[vi] != vno; vi--);                    // find it, backwards since the active ones are at end
                do { n->vnos[vi] = n->vnos[vi+1]; } while (++vi < n->vnosSize - 1);     // remove it
            }

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

    // process all the pending fno, resetting the links to zero
    //    
    int fno = firstFnoToUpdatePlus1 - 1;
    while (fno >= 0) {

        RealmTreeNode * chosenForFno = chooseRealmTreeNodeForFno(mris, realmTree, fno);
    
        if (chosenForFno != realmTree->fnoToRealmTreeNode[fno]) {
            removeFnoFromRealmTree(realmTree, fno);
            insertFnoIntoRealmTreeNode(realmTree, chosenForFno, fno);
        }
                
        int next_fno = realmTree->nextFnoToUpdatePlus1[fno] - 1;        // get the next fno, if any
        realmTree->nextFnoToUpdatePlus1[fno] = 0;                       // clear the link for next time
        fno = next_fno;
    }
}

void checkRealmTree(RealmTree const * realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ) {
    int count = countXYZChanges(realmTree, mris, getXYZ);
    if (count > 0) {
        fprintf(stderr, "%s:%d mris %d vertex xyz have changed\n", __FILE__, __LINE__, count);
        // DON'T EXIT FOR NOW  exit(1);
    }
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
    free(*realmPtr);
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
        if (!n->vnos) bevins_break();   // must be a leaf
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

static int numberOffnosHereAndDeeper(RealmTreeNode const* n) {
    int count = n->nFaces;
    if (!n->vnos) {
        int c;
        for (c = 0; c < childrenSize; c++) count += numberOffnosHereAndDeeper(n->children[c]);
    }
    return count; 
}

int realmNumberOfMightTouchFno(Realm* realm) {
    RealmTreeNode const* n = realm->deepestContainingNode;
    if (!n) return 0;
    int count = numberOffnosHereAndDeeper(n);
    while ((n = n->parent)) count += n->nFaces;
    return count;  
}

static int fnosHere(RealmTree const* rt, RealmTreeNode const* n, int* fnos, int fnosCapacity, int fnosSize) {
    int fno = n->firstFnoPlus1 - 1;
    while (fno >= 0) {
#ifdef REALM_UNIT_TEST
        if (fnosCapacity <= fnosSize) {
            bevins_break();
        }
#endif
        fnos[fnosSize++] = fno;
        fno = rt->nextFnoPlus1[fno] - 1;
    }
    return fnosSize;
}

static int fnosHereAndDeeper(RealmTree const* rt, RealmTreeNode const* n, int* fnos, int fnosCapacity, int fnosSize) {
    fnosSize = fnosHere(rt, n, fnos, fnosCapacity, fnosSize);
    if (!n->vnos) {
        int c;
        for (c = 0; c < childrenSize; c++) fnosSize = fnosHereAndDeeper(rt, n->children[c], fnos, fnosCapacity, fnosSize);
    }
    return fnosSize; 
}

int realmMightTouchFno(Realm* realm, int* fnos, int fnosCapacity) {
    RealmTreeNode const* n = realm->deepestContainingNode;
    if (!n) return 0;
    int written = fnosHereAndDeeper(realm->realmTree, n, fnos, fnosCapacity, 0);
    while ((n = n->parent)) written = fnosHere(realm->realmTree, n, fnos, fnosCapacity, written);
    return written;
}

static void summarizeRealmTreeNode(RealmTreeNode const * n) {
    int i; for (i = 0; i < n->depth; i++) printf("   |");
    printf("x:%f..%f y:%f..%f z:%f..:%f nFaces:%d", n->xLo, n->xHi, n->yLo, n->yHi, n->zLo, n->zHi, n->nFaces);
    if (n->vnos) {
        printf(" nosSize:%d\n",n->vnosSize);
    } else {
        printf("\n");
        int c;
        for (c = 0; c < childrenSize; c++) summarizeRealmTreeNode(n->children[c]);
    }
}

void summarizeRealmTree(RealmTree const * realmTree) {
    summarizeRealmTreeNode(&realmTree->root);
}
