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
#include <stdbool.h>
#include "fnv_hash.h"

// RealmTree construction and destruction
//
typedef struct RealmTreeNode RealmTreeNode;
struct RealmTreeNode {
    float xLo, xHi, yLo, yHi, zLo, zHi;
    RealmTreeNode* parent;
    
#define VNOS_CAPACITY 8
    int vnoSize;
    union {
        RealmTreeNode* children[8]; // 2x2x2, when vnoSize  > VNOS_CAPACITY
        int vnos[VNOS_CAPACITY];    //        when vnoSize <= VNOS_CAPACITY
    };
};

static bool nodeContains(
    RealmTreeNode*  const n, 
    float const x, float const y, float const z) {
    return  n->xLo <= x && x < n->xHi &&
            n->yLo <= y && y < n->yHi &&
            n->zLo <= z && z < n->zHi;
}

static RealmTreeNode* upUntilContainsNode(RealmTreeNode* const n, 
    float const x, float const y, float const z) {
    while (!nodeContains(n,x,y,z)) n = n->parent;
    return n;
}

static int chooseChild(
    RealmTreeNode* n,
    float x, float y,float z) {
    int c = ((x < n->children[1]->xLo) ? 0 : 1) +
            ((x < n->children[2]->yLo) ? 0 : 1) +
            ((x < n->children[4]->zLo) ? 0 : 1);
    return c;
}

static RealmTreeNode* deepestContainingNode(RealmTreeNode* const n, 
    float const x, float const y, float const z) {
    n = upUntilContainsNode(n, x, y, z);
    while (n->vnoSize > VNOS_CAPACITY) {
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

static void insertIntoNode(
    RealmTree*      const realmTree,
    RealmTreeNode*  const n, 
    int const vno, float const x, float const y, float const z)
{
    MRIS const* mris = realmTree->mris;
    
    // Can fit in this node
    if (n->vnoSize < VNOS_CAPACITY) {
        n->vnos[n->vnoSize++] = vno;
        realmTree->vnoToRealmTreeNode[vno] = n;
        return n;
    }

    // Has already or need to subdivide
    if (n->vnoSize == VNOS_CAPACITY) {

        // Chose the splitting values
        float xMid = (n->xLo + n->xHi)/2;
        float yMid = (n->yLo + n->yHi)/2;
        float zMid = (n->zLo + n->zHi)/2;

        // Save the vnos, since the next step overwrites them
        //
        int vnos[VNOS_CAPACITY];
        int vi;
        for (vi = 0; vi < VNOS_CAPACITY; vi++) {
            vnos[vi] = n->children[vi];
        } 
        
        // Make the children
        int c;
        for (c = 0; c < 8; c++) { 
            RealmTreeNode* child = n->children[c] = 
                (RealmTreeNode*)calloc(1, sizeof(RealmTreeNode));
            if (c&1) c->xLo = xMid, c->xHi = n->xHi; else c->xLo = xLo, c->xHi = xMid; 
            if (c&2) c->yLo = yMid, c->yHi = n->yHi; else c->yLo = yLo, c->yHi = yMid;
            if (c&4) c->zLo = zMid, c->zHi = n->zHi; else c->zLo = zLo, c->zHi = zMid;
        }
        n->vnoSize = VNOS_CAPACITY + 1;
        
        // Insert the saved vno into their right child
        int vi;
        for (vi = 0; vi < 8; vi++) {
            RealmTreeNode* tmpN     = n; 
            int            tmpDepth = depth; 
            VERTEX const * vertex = &mris->vertices[vnos[vi]];
            insertIntoChild(realmTree, n, vertex->x, vertex->y, vertex->z);
        }
    }

    // Insert this vno into the right child
    return insertIntoChild(realmTree, n, x, y, z);
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
    free(rt->vnoToRealmTreeNode);
    free(rt);
}

RealmTree* makeRealmTree(MRIS const * mris) {
    // Fills in the tree using the existing position of 
    // the vertices and faces
    RealmTree* rt = (RealmTree*)calloc(1, sizeof(RealmTree));
    rt->mris     = mris;
    rt->fnv_hash = computeRealmTreeHash(mris);
    rt->vnoToRealmTreeNode = (RealmTreeNode**)calloc(mris->nvertices, sizeof(RealmTreeNode*));

    if (mris->nvertices == 0) return;
    
    // Calculate the outer box
    //
    int vno = 0;
    VERTEX const * vertex0 = &mris->vertices[vno];
    int xLo = vertex->x, xHi = xLo,
        yLo = vertex->y, yHi = yLo,
        zLo = vertex->z, zHi = zHi;
    for (vno = 1; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        int const x = vertex->x, y = vertex->y, z = vertex->z;
        xLo = MIN(xLo, x); yLo = MIN(yLo, y); zLo = MIN(zLo, z); 
        xHi = MAX(xHi, x); yHi = MAX(yHi, y); zHi = MAX(zHi, z); 
    }

    RealmTreeNode* recentNode  = &rt->root;
    recentNode->xLo = xLo; recentNode->yLo = yLo; recentNode->zLo = zLo;
    recentNode->xHi = xHi; recentNode->yHi = yHi; recentNode->zHi = zHi;

    // Place all the vertices into boxes
    // There are up to 8 boxes at each level
    // Assuming 100,000 vno's we want to subdivide into about 5 levels
    //      about   8x8     64
    //             x8x8   4096
    //             x8    32768  boxes, so about 3 per box
    // 
    for (vno = 0; vno < mris->nvertices; vno++) {
        VERTEX const * vertex = &mris->vertices[vno];
        int const x = vertex->x, y = vertex->y, z = vertex->z;
        // Find the right subtree
        while (!nodeContains(recentNode, x,y,z)) {
            recentNode = recentNode->parent;
        }
        // Insert here, or deeper
        recentNode = insertIntoNode(recentNode, vno, x,y,z);
    }

    return rt;
}

void checkRealmTree(RealmTree* realTree, MRIS const * mris) {
    unsigned long hash_now = computeRealmTreeHash(mris);
    if (rt->fnv_hash != hash_now) {
        printf(stderr, "%s:%d mris some vertex xyz has changed\n");
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
}


// Quick tests
//
bool realmMightTouchFno(Realm const * realm, int fno) {
    fprintf(stderr, "%s:%d NYI\n", __FILE__, __LINE__);
    exit(1);
    return false;
}

bool realmMightTouchVno(Realm const * realm, int vno) {
    RealmTreeNode* realmNode = realm->deepestContainingNode;
    RealmTreeNode* vnoNode   = realm->realmTree->vnoToRealmTreeNode[vno];
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
int realmNextMightTouchFno(Realm* realm, int & realmIterator) {
    TBD;
}

int realmNextMightTouchVno(Realm* realm, int & realmIterator) {
    TBD;
}
