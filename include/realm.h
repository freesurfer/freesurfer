#ifndef REALM_H
#define REALM_H

/**
 * @file  realm.h
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

#define REALM_UNIT_TEST

#include <stdbool.h>

#ifndef REALM_UNIT_TEST
#include "mrisurf.h"
#else
    typedef struct VERTEX VERTEX;
    typedef struct MRIS MRIS;
#endif

typedef struct RealmTree RealmTree;
typedef void (*GetXYZ_FunctionType)(VERTEX const * vertex, float* x, float* y, float* z);
void freeRealmTree(RealmTree** realmTreePtr);
RealmTree* makeRealmTree(MRIS const * mris, 
    GetXYZ_FunctionType getXYZ  // This lets realms be on x,y,z, origx,origy,origz, or anything else...
    );
void checkRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ);
    //
    // Fills in the tree using the existing position of 
    // the vertices and faces.  The check version verifies
    // that the faces and vertices have not moved since they were 
    // used to make the tree.

void noteIfXYZChangedRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ, int vno);
void updateRealmTree(RealmTree* realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ);

void getRealmTreeBnds(
    RealmTree* realmTree, float* xLo, float*xHi, float* yLo, float* yHi, float* zLo, float* zHi);

typedef struct Realm Realm;
void freeRealm(Realm** realmPtr);
Realm* makeRealm(RealmTree const * realmTree, 
    float xLo, float xHi, 
    float yLo, float yHi,
    float zLo, float zHi);
    //
    // Creates a realm that can be used to quickly find vertices and faces
    // that MIGHT intersect this brick
    
typedef struct RealmIterator {  // can be assigned safely
    unsigned long i;
    void*         p;
} RealmIterator;
void initRealmIterator(RealmIterator* realmIterator, Realm* realm);
    // no fini needed

int realmNextMightTouchVno(Realm* realm, RealmIterator* realmIterator);
    // first call with realmIterator init'ed by initRealmIterator
    // successive calls return some next fno or vno, may not be ascending order!
    // updates realmIterator to some private state
    // returns -1 when no more found, and on further calls after this

int realmNumberOfMightTouchFno(Realm* realm);
int realmMightTouchFno(Realm* realm, int* fnos, int fnosCapacity);

void summarizeRealmTree(RealmTree const * rt);
    // Writes a summary to stdout


#endif
