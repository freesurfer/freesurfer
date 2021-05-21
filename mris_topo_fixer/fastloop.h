/*
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


#ifndef TOPOLOGY_FASTLOOP_H
#define TOPOLOGY_FASTLOOP_H

//#pragma  warning( disable : 4702 )
#include <queue>
#include <functional>
#include <climits>
//#pragma  warning( default : 4702 )

#include "globals.h"
#include "surface.h"
#include "segment.h"
#include "loop.h"

class FastLoop
{
private:
  typedef enum {eAlive=0, eTrial=1, eFar=2, eForbidden=3, eTemporary=4 } eState;

  //////////////////////////////////////////////////////////////////////////////
  //information about the surface
  class FaceData
  {
  public:
    int v[3]; // the 3 neighboring vertices
    int f[3]; // the 3 neighboring faces
    double x,y,z;
    bool border; // bordering face or not
    //Fast Marching on the Faces
    int fmState; // state of the face in the fast marching (eAlive, eTrial, eFar or eForbidden
    int nfather; // number of face fathers in the fast marching
    int fmFather; //the first father in the fast marching
    double val;
    double val2;
    // Fast Segmentation Into Connected Components
    int fccState; // state of the face in the fast segmentation into connected components (eAlive, eTrial, eFar or eForbidden
    int fccFather; // the first father in the fast segmentation into connected components
    int fccLabel;
    //constructor
    FaceData()
    {
      fmState = eForbidden;
      fccState = eForbidden;
      border = 1;
    };
  };
  class VertexData
  {
  public:
    int fmFather;
    int fmState;
    VertexData()
    {
      fmState = eForbidden;
      fmFather = -2;
    }
  };

  ////////////////////////////////////////////////////////////////////////////
  // Fast Marching class
class HeapCompare : std::binary_function<int,int,bool>
  {
  protected:
    FaceData *f;
  public:
    HeapCompare(FaceData *_f) : f(_f)
    {}
    bool operator() (const int a, int b) const
    {
      return (f[a].val>f[b].val);
    }
  };

  typedef std::priority_queue<int,std::vector<int>,HeapCompare> FaceHeap;

  //////////////////////////////////////////////////////////////////////////////
  //information about the surface
  Surface *surface;
  //Vertices *vertices;
  //additional information about the surface
  FaceData *facedata;
  VertexData *vertexdata;

  //list of defect faces
  int *defect_faces;
  int ndefect_faces;

  //Fast Marching
  FaceHeap *FM_trial_heap;

  //Fast Segmentation
  FaceHeap *FCC_trial_heap;
  int final_face[2];
  int nsegments[2];
  Segment segments[2];

  //Functions
  void _InitFaceData(void);
  void _UpdateFace(int fdst, int fsrc);
  void _UpdateSegmentFace(int fdst, int fsrc);
  void _InitSegment(int which_segment, int fn);
  int _InitSegmentation(int fn);
  int _FindSeedFaces(int conflicting_face, int& init_fn1, int& init_fn2);
  void _UpdateSegmentFaces(int which_segment);
  double _Distance(int fdst, int fsrc);
  int _CheckAdjacency(int fn);
  int _FastSegmentation();
  int _AddAliveFace(int fno);
  int _AddTrialFace(int fno);
  int _Run(int &stopping_face);
  void _InitDefect();
  int _ExtractFirstLoop(Loop& loop, int init_fn1, int init_fn2);
  int _FindFacePath(Loop &loop,int init_fn1,int init_fn2);
  int _FindCommonVertex(int init_fn1,int init_fn2);
  int _FindNextFace(int next_fn,int vno);
  int _ExtractSecondLoop(Loop& loop, int init_fn);
  int _SimplifyLoop(Loop &loop);
  int _OrderLoop(Loop &loop);
  double _GetLoopLength(Loop &loop);

public:
  //////////////////////////////////////////////////////////////////////////////
  //constructor / destructor
  FastLoop(Surface &s);
  ~FastLoop(void);
  void Init(void);
  void SetSeed(int seed);
  void SetDefectList(int nfaces, int *list_of_faces);
  Loop* FindLoop(int seed);
  void FindMinimalLoop(Loop &minimial_loop , int max_init_face=-1 , int nattempts = 10);
};

#endif
