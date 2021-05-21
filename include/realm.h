#ifndef REALM_H
#define REALM_H

/**
 * @brief support quickly scanning all the vertices or faces on an MRI for
 *        for those that might intersect a brick
 *
 * 
 */
/*
 * Original Author: Bevin R Brett
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


//  I looked at many synonyms for a region, zone, section, sector, subdivision, partition, ... 
//  of a volume into subvolumes but searching the sources forthe other choices found them 
//  being used already.  This particular choice had only one hit, in one comment.
//  
//  Summarized, a RealmTree data structure is populated just before fixing the defects.  
//  It is adjusted as vertexs are moved during defect repair, and discarded when no longer needed.
//  
//  A RealmTree the xyz volume of the surface is recursively broken into 2x2x2 subvolumes, 
//  each vertex is listed in the deepest subvolume it fits in.  Each face is listed in the 
//  deepest subvolume that all its vertexs fit in.  
//  
//  Prior to fixing the defects a RealmTree is built.  The cost of the build is therefore 
//  amortized over many uses.   
//  
//  A Realm is a (xLo..xHi, ...) subset of the whole coordinate space.  It does not have to be
//  a specific node in the RealmTree.
//  
//  To find all the faces that might intersect a realm, a walk of the RealmTree visits only those
//  RealmTree nodes that intersect the Realm.   Typically this produces a lisst that is about 
//  1/20th of the all faces, thus dramatically speeding up the loop that chooses which faces to 
//  consider.


//  #define REALM_UNIT_TEST

#include "base.h"

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
int checkRealmTree(RealmTree const * realmTree, MRIS const * mris, GetXYZ_FunctionType getXYZ);
    //
    // Fills in the tree using the existing position of 
    // the vertices and faces.  The check version verifies
    // that the faces and vertices have not moved since they were 
    // used to make the tree.
    //
    // Returns the number of errors

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

void summarizeRealmTree   (RealmTree const * realmTree);
void summarizeRealmTreeFno(RealmTree const * realmTree, int fno);
void summarizeRealmTreeVno(RealmTree const * realmTree, int vno);
    // Writes a summary to stderr


// GreatArcs supports projecting a set of edges onto the unit circle, 
// and then quickly asking whether one projected edge intersects any others
//
// This is needed when fixing defects in mris_fix_topology
//
typedef struct GreatArcSet GreatArcSet;

void freeGreatArcSet(GreatArcSet** setPtr);
GreatArcSet* makeGreatArcSet(MRIS* mris);

void insertGreatArc(GreatArcSet* set, 
    int key,                        // key must not be in the set already
    int vno0,
    int vno1);
    
void possiblyIntersectingGreatArcs(
    GreatArcSet* set,
    void*        callbackCtx,
    bool (*callback)(void* callbackCtx, int key, bool* isHit),      // returns true if should keep going, false if no more needed
    int vno0, int vno1,
    float x0, float y0, float z0,   // ends of the arc, need not be a unit vector
    float x1, float y1, float z1,
    bool tracing);

void possiblyIntersectingGreatArcs_Debug(                           // show how vno0..vno1 interacts with the arc
    GreatArcSet* set,
    float x0, float y0, float z0,   // ends of the arc, need not be a unit vector
    float x1, float y1, float z1,
    int vno0, int vno1);
    
#endif
