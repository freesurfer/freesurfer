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

#include "gcut.h"

CCubeNode::CCubeNode()
{
  thisSubX = 0;
  thisSubY = 0;
  thisSubZ = 0;
  thisIndex = 0;

  thisMean = 0;
  thisVariance = 0;
  thisVariance8V = 0;

  thisPtrNextNode = NULL;
}

CCubeNode::CCubeNode(int subX, int subY, int subZ, int index, double mean)
{
  thisSubX = subX;
  thisSubY = subY;
  thisSubZ = subZ;
  thisIndex = index;

  thisMean = mean;
  thisVariance = 0;
  thisVariance8V = 0;

  thisPtrNextNode = NULL;
}

CCubeNode::CCubeNode(int subX, int subY, int subZ, 
                     int index, double mean, 
                     double variance, double variance8V)
{
  thisSubX = subX;
  thisSubY = subY;
  thisSubZ = subZ;
  thisIndex = index;

  thisMean = mean;
  thisVariance = variance;
  thisVariance8V = variance8V;

  thisPtrNextNode = NULL;
}

CCubeNode::~CCubeNode()
{ }

void CCubeNode:: SetMean(double mean)
{
  if (mean>0)
    thisMean = mean;
  else
    printf("Invalid mean value. It must greater than 0\n");
}

void CCubeNode::SetVariance(double variance)
{
  if (variance>0)
    thisVariance = variance;
  else
    printf("Invalid variance value. It must greater than 0\n");
}

void CCubeNode::SetVariance8V(double variance8V)
{
  if (variance8V>0)
    thisVariance8V = variance8V;
  else
    printf("Invalid variance value. It must greater than 0\n");
}

void CCubeNode::SetNextNode(CCubeNode* ptrNextNode)
{
  thisPtrNextNode = ptrNextNode;
}

int CCubeNode::GetSubX(void)
{
  return thisSubX;
}

int CCubeNode::GetSubY(void)
{
  return thisSubY;
}

int CCubeNode::GetSubZ(void)
{
  return thisSubZ;
}

int CCubeNode::GetIndex(void)
{
  return thisIndex;
}

double CCubeNode::GetMean(void)
{
  return thisMean;
}

double CCubeNode::GetVariance(void)
{
  return thisVariance;
}

double CCubeNode::GetVariance8V(void)
{
  return thisVariance8V;
}

CCubeNode* CCubeNode::GetNextNode(void)
{
  return thisPtrNextNode;
}

CQueue3D::CQueue3D()
{
  thisHead = NULL;
  thisTail = NULL;
}

CQueue3D::~CQueue3D()
{}

void CQueue3D::Enqueue(CCubeNode* ptrNewNode)
{
  if (thisTail != NULL)
  {
    (*thisTail).SetNextNode(ptrNewNode);
    thisTail = ptrNewNode;
  }
  else
  {
    thisHead = ptrNewNode;
    thisTail = ptrNewNode;
  }
}

CCubeNode* CQueue3D::Dequeue(void)
{
  if (thisHead != NULL)
  {
    CCubeNode* tempNode = thisHead;

    if (thisHead == thisTail)
    {
      thisHead = NULL;
      thisTail = NULL;
    }
    else
    {
      thisHead = (*thisHead).GetNextNode();
    }

    return tempNode;
  }
  else
  {
    return NULL;
  }
}

CCubeNode* CQueue3D::GetHead(void)
{
  return thisHead;
}


// -- Function prototypes
double Mean5x5x5(int startX, int startY, int startZ, 
                 int iSizeX, int iSizeY, 
                 unsigned char ***Mat);
double Variance5x5x5(int startX, int startY, int startZ, 
                     int iSizeX, int iSizeY, 
                     double mean, unsigned char ***Mat);
double Variance8V(int startX, int startY, int startZ, 
                  int iSizeX, int iSizeY, 
                  unsigned char ***Mat, int size);
double Mean3x3x3(int startX, int startY, int startZ, 
                 int iSizeX, int iSizeY, unsigned char ***Mat);
double Variance3x3x3(int startX, int startY, int startZ, 
                     int iSizeX, int iSizeY, 
                     double mean, unsigned char ***Mat);
void   SetRegion3x3x3(int startX, int startY, int startZ, 
                      int iSizeX, int iSizeY, 
                      unsigned char ***Mat, short value);
int    Sub2Ind3D(int sX, int sY, int sZ, 
                 int iSizeX, int iSizeY);

// -- Main function
int AutoSelectSeed(gc_POS *iSeed, unsigned char ***Mat, 
                   int xVol, int yVol, int zVol, 
                   double mT, double vT)
{
  // -- Get the pointer to the input argument(s)
  double MeanThreshold = mT;
  double VarianceThreshold = vT;

  // -- Get the size of the image
  int iSizeX = xVol;
  int iSizeY = yVol;
  int iSizeZ = zVol;

  double mean = 0;
  double variance = 0;
  double variance8V = 0;
  double weightedVariance = 100000;

  int tempX = 0;
  int tempY = 0;
  int tempZ = 0;
  double tempMean = 0;
  double tempVariance = 0;

  bool qualifiedFlag = true;

  const int size = 5;

  double shortListMean[10] =
    {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
  double shortListVariance[10] =
    {
      100000, 100000, 100000, 100000, 100000, 
      100000, 100000, 100000, 100000, 100000
    };
  int shortListSeedX[10] =
    {
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
    };
  int shortListSeedY[10] =
    {
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
    };
  int shortListSeedZ[10] =
    {
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
    };

  int neighbor6[6][3] =
    {
      {
        size,0,0
      }
      , {-size,0,0}, {0,size,0}, {0,-size,0}, {0,0,size}, {0,0,-size}
    };

  // -- Divide the whole image to 5x5x5 cubes and select one as the seed
  for (int sX=size; sX<(iSizeX-size); sX=sX+size)
  {
    for (int sY=size; sY<(iSizeY-size); sY=sY+size)
    {
      for (int sZ=size; sZ<(iSizeZ-size); sZ=sZ+size)
      {
        mean = Mean5x5x5(sX, sY, sZ, iSizeX, iSizeY, Mat);
        variance = Variance5x5x5(sX, sY, sZ, iSizeX, iSizeY, mean, Mat);
        variance8V = Variance8V(sX, sY, sZ, iSizeX, iSizeY, Mat, 5);

        if ( (mean>MeanThreshold) && (variance<VarianceThreshold) )
        {
          // -- check whether its neighbors are qualified for seed or not
          for (int n1=0; n1<6; n1++)
          {
            tempX = sX + neighbor6[n1][0];
            tempY = sY + neighbor6[n1][1];
            tempZ = sZ + neighbor6[n1][2];
            //boundary check again!
            if ( tempX >= iSizeX-size || 
                 tempY >= iSizeY-size || 
                 tempZ >= iSizeZ-size )
            {
              qualifiedFlag = false;
              break;
            }
            tempMean = Mean5x5x5(tempX, tempY, tempZ, iSizeX, iSizeY, Mat);
            tempVariance = Variance5x5x5(tempX, tempY, tempZ, 
                                         iSizeX, iSizeY, tempMean, Mat);

            if (!( (tempMean>MeanThreshold) && 
                   (tempVariance<VarianceThreshold) ))
            {
              qualifiedFlag = false;
              break;
            }
            else
            {
              qualifiedFlag = true;
            }
          }

          // -- IF the currently examined seed qualifies to be considered
          if (qualifiedFlag == true)
          {
            weightedVariance = (0.5*variance) + (0.5*variance8V);
            for (int i=0; i<10; i++)
            {
              if (shortListVariance[i]>=weightedVariance)
              {
                for (int j=9; j>i;j--)
                {
                  shortListMean[j] = shortListMean[j-1];
                  shortListVariance[j] = shortListVariance[j-1];
                  shortListSeedX[j] = shortListSeedX[j-1];
                  shortListSeedY[j] = shortListSeedY[j-1];
                  shortListSeedZ[j] = shortListSeedZ[j-1];
                }
                shortListMean[i] = mean;
                shortListVariance[i] = weightedVariance;
                shortListSeedX[i] = sX;
                shortListSeedY[i] = sY;
                shortListSeedZ[i] = sZ;
                break;
              }
            }
          }
        }
      }
    }
  }

  // -- choose the one with lowest intensity
  double lowestMean = 10000000;
  int target = 0;
  for (int k=0; k<10; k++)
  {
    if (shortListMean[k]<lowestMean)
    {
      lowestMean = shortListMean[k];
      target = k;
    }
  }

  int seedX = shortListSeedX[target];
  int seedY = shortListSeedY[target];
  int seedZ = shortListSeedZ[target];

  // -- Create the array for the output argument(s)
  iSeed ->x = seedX;
  iSeed ->y = seedY;
  iSeed ->z = seedZ;

  return 0;
}

// -- Main function
int RegionGrowing(gc_POS iSeed, 
                  unsigned char ***MImage, 
                  unsigned char ***MLabel, 
                  int xVol, int yVol, int zVol, 
                  double lmdT, 
                  double umdT, 
                  double nmdT, 
                  double vT)
{

  double ptrLMeanDiffThreshold = lmdT;
  double ptrUMeanDiffThreshold = umdT;
  double ptrNMeanDiffThreshold = nmdT;
  double ptrVarianceThreshold  = vT;

  // -- Get the size of the image
  int iSizeX = xVol;
  int iSizeY = yVol;
  int iSizeZ = zVol;

  // -- Get the seed position of the image
  int iSeedPosX = iSeed.x;
  int iSeedPosY = iSeed.y;
  int iSeedPosZ = iSeed.z;

  // -- create the label array, which is to be returned as the result of region growing
  //double *ptrILabel = (double*) malloc(iSizeX*iSizeY*iSizeZ*sizeof(double));

  unsigned char ***ptrVisited;
  ptrVisited = new unsigned char**[iSizeZ];
  for (int i = 0; i < iSizeZ; i++)
  {
    ptrVisited[i] = new unsigned char*[iSizeY];
    for (int j = 0; j < iSizeY; j++)
    {
      ptrVisited[i][j] = new unsigned char[iSizeX];
      for (int k = 0; k < iSizeX; k++)
      {
        ptrVisited[i][j][k] = 0;
      }
    }
  }

  // -- Retrieve the seed's parameters and put it in a queue
  double meanSeed = Mean3x3x3(iSeedPosX, iSeedPosY, iSeedPosZ, 
                              iSizeX, iSizeY, MImage);
  double varianceSeed = Variance3x3x3(iSeedPosX, iSeedPosY, iSeedPosZ, 
                                      iSizeX, iSizeY, meanSeed, MImage);
  double variance8VSeed = Variance8V(iSeedPosX, iSeedPosY, iSeedPosZ, 
                                     iSizeX, iSizeY, MImage, 3);
  int iSeedIndex = Sub2Ind3D(iSeedPosX, iSeedPosY, iSeedPosZ, 
                             iSizeX, iSizeY);

  CCubeNode *seedNode = new CCubeNode(iSeedPosX, iSeedPosY, iSeedPosZ, 
                                      iSeedIndex, meanSeed, 
                                      varianceSeed, variance8VSeed);
  CQueue3D  *ptrQueue = new CQueue3D();
  ptrQueue->Enqueue(seedNode);

  // -- Set as visited
  SetRegion3x3x3(iSeedPosX, iSeedPosY, iSeedPosZ, 
                 iSizeX, iSizeY, MLabel, 1);
  SetRegion3x3x3(iSeedPosX, iSeedPosY, iSeedPosZ, 
                 iSizeX, iSizeY, ptrVisited, 1);

  // -- compute the max and min intensity of the pixels
  double minMean = meanSeed - ptrLMeanDiffThreshold;
  double maxMean = meanSeed + ptrUMeanDiffThreshold;

  // -- Traverse the tree structure using the queue
  CCubeNode *ptrCurrentNode = ptrQueue->Dequeue();

  double meanCurrentNode = 0;
  //double varianceCurrentNode = 0;
  //double variance8VCurrentNode = 0;
  int subXCurrentNode = 0;
  int subYCurrentNode = 0;
  int subZCurrentNode = 0;
  //int indexCurrentNode = 0;

  double meanNeighborNode = 0;
  double varianceNeighborNode = 0;
  double variance8VNeighborNode = 0;
  int subXNeighborNode = 0;
  int subYNeighborNode = 0;
  int subZNeighborNode = 0;
  int indexNeighborNode = 0;

  //double weightedVariance = 0;
  //int loopControlCount = 0;

  int sixNeighbourOffset[6][3] = {
    {-1,0,0},
    {1,0,0},
    {0,-1,0},
    {0,1,0},
    {0,0,-1},
    {0,0,1}};

  while (ptrCurrentNode != NULL)
  {

    // -- process the current node
    meanCurrentNode = (*ptrCurrentNode).GetMean();
    subXCurrentNode = (*ptrCurrentNode).GetSubX();
    subYCurrentNode = (*ptrCurrentNode).GetSubY();
    subZCurrentNode = (*ptrCurrentNode).GetSubZ();

    // -- expand to the 6-neighbor of the current node
    for (int n1=0; n1<6; n1++)
    {
      subXNeighborNode = subXCurrentNode + sixNeighbourOffset[n1][0];
      subYNeighborNode = subYCurrentNode + sixNeighbourOffset[n1][1];
      subZNeighborNode = subZCurrentNode + sixNeighbourOffset[n1][2];

      if ( subXNeighborNode <= 3 || 
           subXNeighborNode >= iSizeX-3 || 
           subYNeighborNode <= 3 ||
           subYNeighborNode >= iSizeY-3 || 
           subZNeighborNode <= 3 || 
           subZNeighborNode >= iSizeZ-3 ) // skip the pixel out of the image
        continue;

      indexNeighborNode = Sub2Ind3D(subXNeighborNode, 
                                    subYNeighborNode, 
                                    subZNeighborNode, 
                                    iSizeX, iSizeY);

      if ( ptrVisited[subZNeighborNode][subYNeighborNode][subXNeighborNode] 
           != 1)
      {
        meanNeighborNode = Mean3x3x3(subXNeighborNode, 
                                     subYNeighborNode, 
                                     subZNeighborNode, 
                                     iSizeX, iSizeY, MImage);
        varianceNeighborNode = Variance3x3x3(subXNeighborNode, 
                                             subYNeighborNode, 
                                             subZNeighborNode, 
                                             iSizeX, iSizeY, 
                                             meanNeighborNode, MImage);
        variance8VNeighborNode = Variance8V(subXNeighborNode, 
                                            subYNeighborNode, 
                                            subZNeighborNode, 
                                            iSizeX, iSizeY, 
                                            MImage, 3);

        //weightedVariance = (0.5*varianceNeighborNode) + (0.5*variance8VNeighborNode);
        SetRegion3x3x3(subXNeighborNode, 
                       subYNeighborNode, 
                       subZNeighborNode, 
                       iSizeX, iSizeY, 
                       ptrVisited, 1);

        if ( fabs(meanCurrentNode-meanNeighborNode) < ptrNMeanDiffThreshold
             && varianceNeighborNode < ptrVarianceThreshold
             //&& variance4VNeighborNode < *ptrVarianceThreshold
             //&& weightedVariance < *ptrVarianceThreshold
             && meanNeighborNode > minMean
             && meanNeighborNode < maxMean )
        {
          meanNeighborNode = (meanNeighborNode+meanCurrentNode) / 2;
          ptrCurrentNode = new CCubeNode(subXNeighborNode, 
                                         subYNeighborNode, 
                                         subZNeighborNode, 
                                         indexNeighborNode, 
                                         meanNeighborNode, 
                                         varianceNeighborNode, 
                                         variance8VNeighborNode);
          ptrQueue->Enqueue(ptrCurrentNode);
          SetRegion3x3x3(subXNeighborNode, 
                         subYNeighborNode, 
                         subZNeighborNode, 
                         iSizeX, iSizeY, 
                         MLabel, 1);
        }
      }
    }

    ptrCurrentNode = ptrQueue->Dequeue();
  }

  // -- free memory
  for (int i = 0; i < iSizeZ; i++)
  {
    for (int j = 0; j < iSizeY; j++)
    {
      delete[] ptrVisited[i][j];
    }
    delete[] ptrVisited[i];
  }
  delete[] ptrVisited;
  if (ptrCurrentNode) delete ptrCurrentNode;
  ptrCurrentNode = ptrQueue->Dequeue();
  while (ptrCurrentNode)
  {
    delete ptrCurrentNode;
    ptrCurrentNode = ptrQueue->Dequeue();
  }
  delete ptrQueue;
  delete seedNode;

  return 0;
}

// -- the main pre_processing function
double pre_porocessing(unsigned char ***image, 
                       unsigned char ***label, 
                       int xVol, int yVol, int zVol)
{
  gc_POS iSeed;
  double meanThreshold = 0.45 * 255;
  double varianceThreshold = 0.004 * 255 * 255;

  AutoSelectSeed(&iSeed, image, xVol, yVol, zVol, 
                 meanThreshold, varianceThreshold);
  int repetition = 1;
  while ( iSeed.x == -1 || iSeed.y == -1 || iSeed.z == -1 )
  {
    meanThreshold = meanThreshold - 0.05 * 255;
    AutoSelectSeed(&iSeed, image, xVol, yVol, zVol, 
                   meanThreshold, varianceThreshold);
    if ( repetition++ >= 8 )
      break;
  }
  if ( iSeed.x == -1 || iSeed.y == -1 || iSeed.z == -1 )
  {
    printf("AutoSelectSeed failed!");
    return -1;
  }
  iSeed.x += 2;
  iSeed.y += 2;
  iSeed.z += 2;
  // -- grow
  double meanDiffThreshold = 0.03 * 255;
  //double varianceThreshold = 0.004 * 255 * 255;
  double lmdT = 0.05 * 255;
  double umdT = 0.25 * 255;
  RegionGrowing(iSeed, image, label, xVol, yVol, zVol, 
                lmdT, umdT, meanDiffThreshold, varianceThreshold);
  //mean of white matte
  double whitemean = 0;
  int numSeed = 0;
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        if ( label[z][y][x] == 1 )
        {
          whitemean += image[z][y][x];
          numSeed++;
        }
      }
    }
  }
  if ( numSeed == 0 || whitemean <= 1 )
    printf("no white mean!!\n");
  else
  {
    whitemean = whitemean / numSeed;
    printf("white mean: %f; seed number %d\n", whitemean, numSeed);
  }
  return whitemean;
}

// -- Subroutine to compute the mean of a 5x5x5 cube
double Mean5x5x5(int startX, int startY, int startZ, 
                 int iSizeX, int iSizeY, unsigned char ***Mat)
{
  const int size = 5;
  double sum = 0;
  //int ind = 0;

  for (int x=startX; x<(startX+size); x++)
  {
    for (int y=startY; y<(startY+size); y++)
    {
      for (int z=startZ; z<(startZ+size); z++)
      {
        sum = sum + Mat[z][y][x];
      }
    }
  }

  double mean = sum/(size*size*size);
  return mean;
}

// -- Subroutine to compute the variance of a 5x5x5 cube
double Variance5x5x5(int startX, int startY, int startZ, 
                     int iSizeX, int iSizeY, 
                     double mean, unsigned char ***Mat)
{
  const int size = 5;
  double sum = 0;
  //int ind = 0;
  double pixel = 0;

  for (int x=startX; x<(startX+size); x++)
  {
    for (int y=startY; y<(startY+size); y++)
    {
      for (int z=startZ; z<(startZ+size); z++)
      {
        pixel = Mat[z][y][x];
        sum = sum + ((pixel-mean)*(pixel-mean));
      }
    }
  }

  double variance = sum/(size*size*size);
  return variance;
}

// -- subrountine to compute the variance of the 8 vertices of a 5x5x5 cube
double Variance8V(int startX, int startY, int startZ, 
                  int iSizeX, int iSizeY, 
                  unsigned char ***Mat, int size)
{
  //const int size = 5;
  double sum = 0;
  int i;

  int ind[8];
  double pixel[8];

  ind[0] = Mat[startZ][startY][startX];
  if (startX+size)
    ind[1] = Mat[startZ][startY][startX+size];
  ind[2] = Mat[startZ][startY+size][startX];
  ind[3] = Mat[startZ][startY+size][startX+size];
  ind[4] = Mat[startZ+size][startY][startX];
  ind[5] = Mat[startZ+size][startY][startX+size];
  ind[6] = Mat[startZ+size][startY+size][startX];
  ind[7] = Mat[startZ+size][startY+size][startX+size];

  for (i=0; i<8; i++)
  {
    pixel[i] = ind[i];
    sum = sum + pixel[i];
  }

  double mean = sum/8;

  double variance = 0;
  for (i=0; i<8; i++)
  {
    variance = variance + ((pixel[i]-mean)*(pixel[i]-mean));
  }
  variance = variance/8;
  return variance;
}


// -- Subroutine to compute the mean of a 3x3x3 cube
double Mean3x3x3(int startX, int startY, int startZ, 
                 int iSizeX, int iSizeY, unsigned char ***Mat)
{
  const int size = 3;
  double sum = 0;
  //int ind = 0;

  for (int x=startX; x<(startX+size); x++)
  {
    for (int y=startY; y<(startY+size); y++)
    {
      for (int z=startZ; z<(startZ+size); z++)
      {
        sum = sum + Mat[z][y][x];
      }
    }
  }

  double mean = sum/(size*size*size);
  return mean;
}

// -- Subroutine to compute the variance of a 5x5x5 cube
double Variance3x3x3(int startX, int startY, int startZ, 
                     int iSizeX, int iSizeY, 
                     double mean, unsigned char ***Mat)
{
  const int size = 3;
  double sum = 0;
  //int ind = 0;
  double pixel = 0;

  for (int x=startX; x<(startX+size); x++)
  {
    for (int y=startY; y<(startY+size); y++)
    {
      for (int z=startZ; z<(startZ+size); z++)
      {
        pixel = Mat[z][y][x];
        sum = sum + ((pixel-mean)*(pixel-mean));
      }
    }
  }

  double variance = sum/(size*size*size);
  return variance;
}

// -- Subroutine to set the region as ROI (true);
void SetRegion3x3x3(int startX, int startY, int startZ, 
                    int iSizeX, int iSizeY, 
                    unsigned char ***Mat, short value)
{
  //const int size = 3;
  Mat[startZ][startY][startX] = value;
  return;
}

// -- Subroutine to convert the subscript to linear index
int Sub2Ind3D(int sX, int sY, int sZ, int iSizeX, int iSizeY)
{
  int ind = sX + (sY*iSizeX) + (sZ*iSizeX*iSizeY);
  return ind;
}

