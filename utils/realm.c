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
    RealmTreeNode* children[8]; // 2x2x2
}

struct RealmTree {
    MRIS const  *  mris;
    bool           inited;
    unsigned long  fnv_hash;
    RealmTreeNode  root;
}

void freeRealmTree(RealmTree** realmTreePtr) {
    free(*realmTreePtr);
    *realmTreePtr = NULL;
}

RealmTree* makeRealmTree(MRIS const * mris) {
    // Fills in the tree using the existing position of 
    // the vertices and faces
    RealmTree* rt = (RealmTree*)calloc(1, sizeof(RealmTree));
    rt->mris   = mris;
    rt->inited = false;
    return rt;
}

void checkRealmTree(RealmTree* realTree, MRIS const * mris) {
    TBD;
}

// Realm construction and destruction
//
struct Realm {
    RealmTree const * realmTree;
    float xLo, xHi, yLo, yHi, zLo, zHi;
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
}


// Quick tests
//
bool realmMightTouchFno(Realm const * realm, int fno) {
    TBD;
}

bool realmMightTouchVno(Realm const * realm, int vno) {
    TBD;
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
