/**
 * @brief 
 *
 */
/*
 * Original Author: Vitali Zagorodnov, ZHU Jiaqi
 *
 * Copyright © 2009 Nanyang Technological University, Singapore
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

CCubeNode::CCubeNode(int subX, int subY, int subZ, int index,
                     double mean, double variance, double variance8V)
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
                 int iSizeX, int iSizeY,
                 unsigned char ***Mat);
double Variance3x3x3(int startX, int startY, int startZ,
                     int iSizeX, int iSizeY,
                     double mean, unsigned char ***Mat);
void   SetRegion3x3x3(int startX, int startY, int startZ,
                      int iSizeX, int iSizeY,
                      unsigned char ***Mat, short value);
int    Sub2Ind3D(int sX, int sY, int sZ, int iSizeX, int iSizeY);

// -- Main function
int AutoSelectSeed(gc_POS *iSeed,
                   unsigned char ***Mat,
                   int xVol, int yVol, int zVol,
                   double mT,
                   double vT)
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
                  double lmdT, double umdT, double nmdT, double vT)
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

  // -- create the label array, which is to be returned as
  // the result of region growing
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

  int sixNeighbourOffset[6][3] = {{-1,0,0},
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
                                            iSizeX, iSizeY, MImage, 3);

        //weightedVariance = (0.5*varianceNeighborNode) + (0.5*variance8VNeighborNode);
        SetRegion3x3x3(subXNeighborNode, subYNeighborNode, subZNeighborNode,
                       iSizeX, iSizeY, ptrVisited, 1);

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
                         iSizeX, iSizeY, MLabel, 1);
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
double pre_processing(unsigned char ***image,
                      unsigned char ***label,
                      int xVol, int yVol, int zVol)
{
  gc_POS iSeed;
  double meanThreshold = 0.45 * 160;
  double varianceThreshold = 0.004 * 160 * 160;

  AutoSelectSeed(&iSeed, image, xVol, yVol, zVol,
                 meanThreshold, varianceThreshold);
  int repetition = 1;
  while ( iSeed.x == -1 || iSeed.y == -1 || iSeed.z == -1 )
  {
    meanThreshold = meanThreshold - 0.05 * 160;
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
  double meanDiffThreshold = 0.03 * 160;
  double varianceDiffThreshold = 0.004 * 160 * 160;
  double lmdT = 0.05 * 160;
  double umdT = 0.25 * 160;
  RegionGrowing(iSeed, image, label, xVol, yVol, zVol,
                lmdT, umdT, meanDiffThreshold, varianceDiffThreshold);
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

// direction map function: 3-dimensional 3*3*3
void map(lc_Component p, int & x, int & y, int & z)
{
  int sign = p.current_child;
  z = sign / 9;
  y = (sign - z * 9) / 3;
  x = (sign - z * 9 - y * 3);
  z--;
  y--;
  x--;
  z = z + p.z;
  y = y + p.y;
  x = x + p.x;
}

void push(lc_Component* p, int & current, int x, int y, int z)
{
  current++;
  p[current].x = x;
  p[current].y = y;
  p[current].z = z;
  p[current].current_child = 0;//0~26, skip 13
}

void pop(lc_Component* p, int & current)
{
  current--;
  if (current >= 0)
  {
    p[current].current_child++;
  }
}

//checking function
bool IsLeaf(lc_Component* p, int current, int & x, int & y, int & z,
            unsigned char ***label, int ***marked,
            int xl, int xu, int yl, int yu, int zl, int zu)
{
  if (p[current].current_child >= 27)//all children are done
    return true;

  for (int i = p[current].current_child; i <= 26; i++)
  {
    if (i == 13)
    {
      p[current].current_child = i + 1;
      continue;
    }
    map(p[current], x, y, z);
    if (z >= zl && z <= zu && y >= yl && y <= yu && x >= xl && x <= xu)
    {
      if (label[z][y][x] == 1 && marked[z][y][x] == 0)
        return false;
    }

    p[current].current_child = i + 1;
  }
  return true;
}

//locate a seed from within the LCC
int locateSeedfromLCC(gc_POS & iSeed,
                      unsigned char ***label,
                      int xVol, int yVol, int zVol)
{
  //int xStart, yStart, zStart, xEnd, yEnd, zEnd;
  int bExit = 0;
  int offset = 10;
  //start
  for (int z = offset; z < zVol - offset && bExit == 0; z++)
  {
    for (int y = offset; y < yVol - offset && bExit == 0; y++)
    {
      for (int x = offset; x < xVol - offset && bExit == 0; x++)
      {
        if (label[z][y][x] == 1)
        {
          iSeed.x = x;
          iSeed.y = y;
          iSeed.z = z;
          bExit = 1;
        }
      }
    }
  }

  //reset lable
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        label[z][y][x] = 0;
      }
    }
  }

  return bExit;
}

// -- finding the largest connected component function: when -110 is used
int LCC_function(unsigned char ***image,
                 unsigned char ***label,
                 int xVol, int yVol, int zVol,
                 double & whitemean)
{
  //110 voxels
  int voxelCount = 0;
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        if ( image[z][y][x] == 110 )
        {
          label[z][y][x] = 1;
          voxelCount++;
        }
      }
    }
  }
  //largest connected component
  int ***marked;
  marked = new int**[zVol];
  for (int i = 0; i < zVol; i++)
  {
    marked[i] = new int*[yVol];
    for (int j = 0; j < yVol; j++)
    {
      marked[i][j] = new int[xVol];
      for (int k = 0; k < xVol; k++)
      {
        marked[i][j][k] = 0;
      }
    }
  }
  //extract components
  lc_Component* component_list = new lc_Component[voxelCount];
  int current_comNum = -1;
  int group_id = 0;
  int xPos = -1, yPos = -1, zPos = -1;
  int _limit = 10;
  for (int z = 1; z < zVol - 1; z++)
  {
    for (int y = 1; y < yVol - 1; y++)
    {
      for (int x = 1; x < xVol - 1; x++)
      {
        if ( label[z][y][x] == 1 && marked[z][y][x] == 0 )
        {
          group_id++;
          //LC_func(label, marked, group_id, x, y, z, 1, xVol-10, 1, yVol-1, 1, zVol-1);
          marked[z][y][x] = group_id;
          push(component_list, current_comNum, x, y, z);
          while (current_comNum >= 0)
          {
            while ( !IsLeaf(component_list,
                            current_comNum,
                            xPos, yPos, zPos,
                            label,
                            marked,
                            _limit,
                            xVol - _limit,
                            _limit,
                            yVol - _limit,
                            _limit,
                            zVol - _limit) )
            {
              marked[zPos][yPos][xPos] = group_id;
              push(component_list, current_comNum, xPos, yPos, zPos);
            }
            pop(component_list, current_comNum);
          }
        }
      }
    }
  }
  //find max
  int *max_array = new int[group_id + 1];
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        if (marked[z][y][x] > 0)
        {
          max_array[marked[z][y][x]]++;
        }
      }
    }
  }
  int temp = -1;
  int max_group = 0;
  //find max
  for (int i = 1; i <= group_id; i++)
  {
    if (max_array[i] > temp)
    {
      temp = max_array[i];
      max_group = i;
    }
  }
  //check percentage
  int bNeedGrow = 0;
  if ( max_array[max_group] < 150000 )
    bNeedGrow = 1;
  //LCC
  for (int z = 0; z < zVol; z++)
  {
    for (int y = 0; y < yVol; y++)
    {
      for (int x = 0; x < xVol; x++)
      {
        if ( marked[z][y][x] != max_group )
        {
          label[z][y][x] = 0;
        }
      }
    }
  }
  //free
  delete[] max_array;
  delete[] component_list;
  for (int i = 0; i < zVol; i++)
  {
    for (int j = 0; j < yVol; j++)
    {
      delete[] marked[i][j];
    }
    delete[] marked[i];
  }
  delete[] marked;

  if ( bNeedGrow == 1 )//still need region grow
  {
    gc_POS iSeed;
    iSeed.x = iSeed.y = iSeed.z = 0;
    if ( locateSeedfromLCC(iSeed, label, xVol, yVol, zVol) == 1 )
    {
      // -- grow
      double meanDiffThreshold = 0.03 * 160;
      double varianceDiffThreshold = 0.004 * 160 * 160;
      double lmdT = 0.05 * 160;
      double umdT = 0.25 * 160;
      RegionGrowing(iSeed,
                    image,
                    label,
                    xVol, yVol, zVol,
                    lmdT,
                    umdT,
                    meanDiffThreshold,
                    varianceDiffThreshold);
      //mean of white matte
      //double whitemean = 0;
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
        printf("white mean: %f\n", whitemean);
      }
    }
  }
  return bNeedGrow;

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
  else
    ind[1] = 0; // added to quiet compiler
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
                 int iSizeX, int iSizeY,
                 unsigned char ***Mat)
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

