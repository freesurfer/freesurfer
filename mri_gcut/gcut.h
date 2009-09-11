/*
 * Original Author: Vitali Zagorodnov, ZHU Jiaqi (September, 2009)
 * CVS Revision Info:
 *    $Author:
 *    $Date:
 *    $Revision:
 *
 * Copyright (C) 2009
 * Nanyang Technological University, Singapore
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

typedef struct
{
  int x;
  int y;
  int z;
}
gc_POS;

class CCubeNode
{
public:
  CCubeNode();
  CCubeNode(int subX, int subY, int subZ, int index, double mean);
  CCubeNode(int subX, int subY, int subZ, int index, double mean, double variance, double variance8V);
  ~CCubeNode();

  void SetMean(double mean);
  void SetVariance(double variance);
  void SetVariance8V(double variance8V);
  void SetNextNode(CCubeNode* ptrNextNode);

  int GetSubX(void);
  int GetSubY(void);
  int GetSubZ(void);
  int GetIndex(void);
  double GetMean(void);
  double GetVariance(void);
  double GetVariance8V(void);
  CCubeNode* GetNextNode(void);

private:
  int thisSubX;
  int thisSubY;
  int thisSubZ;
  int thisIndex;

  double thisMean;
  double thisVariance;
  double thisVariance8V;

  CCubeNode *thisPtrNextNode;
};

class CQueue3D
{
public:
  CQueue3D();
  ~CQueue3D();

  void Enqueue(CCubeNode* ptrNewNode);
  CCubeNode* Dequeue(void);
  CCubeNode* GetHead(void);
  CCubeNode* GetTail(void);

private:
  CCubeNode* thisHead;
  CCubeNode* thisTail;
};

