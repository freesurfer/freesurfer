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

#include "fnv_hash.h"


#ifdef REALM_UNIT_TEST
    //
    // rm -f a.out ; gcc -o a.out -I ../include realm.c ; ./a.out
    //
    typedef struct VERTEX {
        float x,y,z;
    } VERTEX;

    struct MRIS {
      int     nvertices;
      VERTEX* vertices;
    };
    
    static float MIN(float lhs, float rhs) { return (lhs < rhs) ? lhs : rhs; }
    static float MAX(float lhs, float rhs) { return (lhs > rhs) ? lhs : rhs; }
    
    static void bevins_break()
    {
        printf("bevins_break\n");
    }

    void test(int nvertices, int useDuplicates) {
        printf("Test nvertices:%d useDuplicates:%d\n", nvertices, useDuplicates);
        
        MRIS mris;
        mris.nvertices = nvertices;
        mris.vertices = (VERTEX*)calloc(mris.nvertices, sizeof(VERTEX));
        
        
        int vno;
        for (vno = 0; vno < mris.nvertices; vno++) {
            int key = vno > useDuplicates ? vno : 936; 
            VERTEX* v = &mris.vertices[vno];
            v->x = (key*321)%51; 
            v->y = (key*7321)%71; 
            v->z = (key*17321)%91;
        }
        vno = 0;
        float xMin = mris.vertices[vno].x, xMax = xMin,
              yMin = mris.vertices[vno].y, yMax = yMin, 
              zMin = mris.vertices[vno].z, zMax = zMin;
        for (vno = 1; vno < mris.nvertices; vno++) {
            VERTEX* v = &mris.vertices[vno];
            xMin = MIN(xMin, v->x);
            yMin = MIN(yMin, v->y);
            zMin = MIN(zMin, v->z);
            xMax = MAX(xMax, v->x);
            yMax = MAX(yMax, v->y);
            zMax = MAX(zMax, v->z);
        }
        
        RealmTree* realmTree = makeRealmTree(&mris);

        int fLimit = 1;
        int fCount = 0;
        float xfLo, xfHi;
        float yfLo, yfHi;
        float zfLo, zfHi;
        for (xfLo = 0; xfLo <= 1; xfLo += 0.1)
        for (xfHi = 0; xfHi <= 1; xfHi += 0.1)
        for (yfLo = 0; yfLo <= 1; yfLo += 0.1)
        for (yfHi = 0; yfHi <= 1; yfHi += 0.1)
        for (zfLo = 0; zfLo <= 1; zfLo += 0.1)
        for (zfHi = 0; zfHi <= 1; zfHi += 0.1)
        {
            float xLo = xMin + xfLo*(xMax-xMin);
            float xHi = xMax - xfHi*(xMax-xLo);
            float yLo = yMin + yfLo*(yMax-yMin);
            float yHi = yMax - yfHi*(yMax-yLo);
            float zLo = zMin + zfLo*(zMax-zMin);
            float zHi = zMax - zfHi*(zMax-zLo);

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
                if (false && (counter == 1 || counter == 122)) {
                    printf("counter:%d ri.i:%ld ri.p:%p\n", counter, realmIterator.i, realmIterator.p); 
                    bevins_break();
                }
                vno = realmNextMightTouchVno(realm, &realmIterator);
                if (vno < -1 || vno >= mris.nvertices) {
                    printf("ERROR, vno:%d is illegal\n", vno); 
                    bevins_break();
                    exit(1);
                }
                if (0 > vno) break;
                if (counter == 0 || states[vno]) {
                    printf("ERROR, vno:%d reported again when counter:%d, was reported counter:%d\n", vno, counter, states[vno]); 
                    bevins_break();
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
                   if (xLo <= v->x && v->x < xHi 
                   &&  yLo <= v->y && v->y < yHi
                   &&  zLo <= v->z && v->z < zHi) printf("ERROR, vno:%d was not reported\n", vno);
                }
            }
            
            free(states);
            freeRealm(&realm);
        }
            
        freeRealmTree(&realmTree);
        
        free(mris.vertices);        
    }
    
    int main() {
        int useDuplicates;
        for (useDuplicates = 0; useDuplicates <= 1000; useDuplicates += 200) {
            test(1, useDuplicates);
            test(2, useDuplicates);
            test(3, useDuplicates);
            test(4, useDuplicates);
            test(100, useDuplicates);
            test(1000, useDuplicates);
            test(10000, useDuplicates);
            test(100000, useDuplicates);
            test(0, useDuplicates);
        }
        return 0;
    }
#endif



// RealmTree construction and destruction
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
            int         vnosSize;                    // above 2 assumes these are the same as the vnosBuff elements
            int         vnosCapacity;
            int         vnosBuff[vnosBuffSize];
        };
    };
};
static const unsigned long childIndexBits =      childrenSizeLog2;
static const unsigned long childIndexMask = ((1<<childrenSizeLog2) - 1);
static const unsigned long leafIndexBits  =      maxVnosSizeLog2;
static const unsigned long leafIndexMask  = ((1<<maxVnosSizeLog2 ) - 1);

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

static int chooseChild(
    RealmTreeNode const * n,
    float x, float y,float z) {

    if (!nodeContains(n, x,y,z)) 
        bevins_break();

    float xMid = n->children[1]->xLo;
    float yMid = n->children[2]->yLo;
    float zMid = n->children[4]->zLo;
    
    int c = ((x < xMid) ? 0 : 1) + ((y < yMid) ? 0 : 2) + ((z < zMid) ? 0 : 4);
    
    if (!nodeContains(n->children[c], x,y,z)) 
        bevins_break();
    
    return c;
}

static RealmTreeNode const * deepestContainingNode(RealmTreeNode const * n, 
    float const x, float const y, float const z) {
    n = upUntilContainsNode(n, x, y, z);
    while (n && !n->vnos) {
        int c = chooseChild(n, x, y, z);
        n = n->children[c];
    }
    return n;
}

static RealmTreeNode* insertIntoNode(
    RealmTree*      const realmTree,
    RealmTreeNode*  const n, 
    int const vno, float const x, float const y, float const z);

static RealmTreeNode* insertIntoChild(
    RealmTree*     realmTree,
    RealmTreeNode* n,
    int vno, float x, float y,float z) {
    int c = chooseChild(n, x, y, z);
    return insertIntoNode(realmTree, n->children[c], vno,x,y,z);
}

struct RealmTree {
    MRIS const  *   mris;
    unsigned long   fnv_hash;
    RealmTreeNode   root;
    RealmTreeNode** vnoToRealmTreeNode;
};

static RealmTreeNode* insertIntoNode(
    RealmTree*      const realmTree,
    RealmTreeNode*  const n, 
    int const vno, float const x, float const y, float const z)
{
    MRIS const* mris = realmTree->mris;

    VERTEX const* v = &mris->vertices[vno];
    if (x != v->x || y != v->y || z != v->z) 
            bevins_break();
    
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
        
        // Insert the saved vno into their right child
        for (vi = 0; vi < vnosSize; vi++) {
            VERTEX const * vertex = &mris->vertices[vnos[vi]];
            insertIntoChild(realmTree, n, vnos[vi], vertex->x, vertex->y, vertex->z);
        }
    }

    // Insert this vno into the right child
    //
    return insertIntoChild(realmTree, n, vno, x, y, z);
}

unsigned long computeRealmTreeHash(MRIS const * mris) {
    unsigned long hash = fnv_init();
    int vno;
    for (vno = 1; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        float f[3]; f[0] = vertex->x, f[1] = vertex->y, f[2] = vertex->z;
        hash = fnv_add(hash, (const unsigned char*)f, sizeof(f));
    }
    return hash;
}


void freeRealmTree(RealmTree** realmTreePtr) {
    RealmTree* rt = *realmTreePtr; *realmTreePtr = NULL;
    if (!rt) return;
    destroyRealmTreeNode(&rt->root);
    free(rt->vnoToRealmTreeNode);
    free(rt);
}

static float widenHi(float lo, float hi) {
    float step = FLT_MIN;
    while (hi + step == hi) {
        step *= 2.0;
    }
    return hi + step;
}

RealmTree* makeRealmTree(MRIS const * mris) {
    // Fills in the tree using the existing position of 
    // the vertices and faces
    RealmTree* rt = (RealmTree*)calloc(1, sizeof(RealmTree));
    rt->mris     = mris;
    rt->fnv_hash = computeRealmTreeHash(mris);
    rt->vnoToRealmTreeNode = (RealmTreeNode**)calloc(mris->nvertices, sizeof(RealmTreeNode*));
    constructRealmTreeNode(&rt->root, NULL);

    if (mris->nvertices == 0) return rt;
    
    // Calculate the outer box
    //
    int vno = 0;
    VERTEX const * vertex0 = &mris->vertices[vno];
    float xLo = vertex0->x, xHi = xLo,
          yLo = vertex0->y, yHi = yLo,
          zLo = vertex0->z, zHi = zLo;
    for (vno = 1; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        float const x = vertex->x, y = vertex->y, z = vertex->z;
        xLo = MIN(xLo, x); yLo = MIN(yLo, y); zLo = MIN(zLo, z); 
        xHi = MAX(xHi, x); yHi = MAX(yHi, y); zHi = MAX(zHi, z); 
    }
    
    // Since contains is xLo <= x < xHi etc. the bounds need to be slightly wider than Hi
    // so that the Hi is in a Node
    //
    xHi = widenHi(xLo,xHi);
    yHi = widenHi(yLo,yHi);
    zHi = widenHi(zLo,zHi);
    
    RealmTreeNode* recentNode  = &rt->root;
    recentNode->xLo = xLo; recentNode->yLo = yLo; recentNode->zLo = zLo;
    recentNode->xHi = xHi; recentNode->yHi = yHi; recentNode->zHi = zHi;

    // Place all the vertices into boxes
    // 
    for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        int const x = vertex->x, y = vertex->y, z = vertex->z;
        // Find the right subtree
        while (!nodeContains(recentNode, x,y,z)) {
            recentNode = recentNode->parent;
        }
        // Insert here, or deeper
        recentNode = insertIntoNode(rt, recentNode, vno, x,y,z);
    }

    return rt;
}

void checkRealmTree(RealmTree* realTree, MRIS const * mris) {
    unsigned long hash_now = computeRealmTreeHash(mris);
    if (realTree->fnv_hash != hash_now) {
        fprintf(stderr, "%s:%d mris some vertex xyz has changed\n", __FILE__, __LINE__);
        exit(1);
    }
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
    // that MIGHT intersect this brick

    Realm* r = (Realm*)calloc(1, sizeof(Realm));
    r->realmTree = realmTree;
    r->xLo = xLo; r->yLo = yLo; r->zLo = zLo;
    r->xHi = xHi; r->yHi = yHi; r->zHi = zHi;
    
    r->deepestContainingNode = 
        upUntilContainsNode(
            deepestContainingNode(&realmTree->root, xLo, yLo, zLo), 
            xHi, yHi, zHi);

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

bool realmMightTouchFno(Realm const * realm, int fno) {
    fprintf(stderr, "%s:%d NYI\n", __FILE__, __LINE__);
    exit(1);
    return false;
}

bool realmMightTouchVno(Realm const * realm, int vno) {
    RealmTreeNode const * realmNode = realm->deepestContainingNode;
    RealmTreeNode const * vnoNode   = realm->realmTree->vnoToRealmTreeNode[vno];
    while (vnoNode) {
        if (vnoNode == realmNode) return true;
        vnoNode = vnoNode->parent;
    }
    return false;
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
