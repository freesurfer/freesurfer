#ifndef ShortestPathFinder_h
#define ShortestPathFinder_h


#include <list>
#include "Point2.h"
#include "Array2.h"
#include "DebugReporter.h"




class listElement {
 public:
  listElement *Prev; 
  listElement *Next; 
  Point2<int> Coord;
  listElement() {this->Prev=NULL; this->Next=NULL;}; 
};


class linkedList : public Array2<listElement>{
 public:
  linkedList ( int izX, int izY );
};
 
class circularQueue {
 public:  
  circularQueue ( int izX, int izY, int icBuckets );  
  ~circularQueue();

  void Insert ( Point2<int>& iLocation, int iCost );
  void Remove ( Point2<int>& iLocation );
  void Remove ( listElement *el );
  listElement *GetListElement ( int iCost );

 private:
  int GetBucket ( int iCost );
  int FindMinBucket ( int iCost );

  linkedList *A;
  listElement *Circle;
  int C;
};



class ShortestPathFinder : public DebugReporter {

 public:

  ShortestPathFinder();
  virtual ~ShortestPathFinder();

  void SetDimensions ( int izX, int izY, int iLongestEdge );

  void SetStraightBias ( float iBias ) { mStraightBias = iBias; }
  void SetEdgeBias ( float iBias ) { mEdgeBias = iBias; }

  void FindPath ( Point2<int>& iStartPoint, Point2<int>& iEndPoint,
		  std::list<Point2<int> >& ioPoints );
  
  virtual float GetEdgeCost ( Point2<int>& iPoint ) { return 1.0; }

  void SetDebug ( bool iDebug ) { mDebug = iDebug; }
  

 protected:
  float mStraightBias;
  float mEdgeBias;
  int mzX, mzY;
  float mLongestEdge;
  circularQueue *mQueue;
  Array2<float> *maCost;
  Array2<int> *maDir;
  Array2<bool> *maDone;
  bool mDebug;
};


#endif
