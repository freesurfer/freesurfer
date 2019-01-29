/**
 * @file  gcut.h
 * @brief 
 *
 */
/*
 * Original Author: Vitali Zagorodnov, ZHU Jiaqi (Sept. 2009, revised Feb. 2010)
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:16 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2009-2010 Nanyang Technological University, Singapore
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

typedef struct
{
  int x;
  int y;
  int z;
}
gc_POS;

typedef struct {
	int x;
	int y;
	int z;
	int current_child;
}lc_Component;

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

