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


#ifndef sparse2DGausDistMatrix_h
#define sparse2DGausDistMatrix_h


#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

#include "mri_ca_gausDist.h"
//using namespace std;


/*

  Main Design Goals for Sparse Matrix:
   1 Must use zero bytes for zero-valued entries
   2 Want to minimize garbage collection thus no reallocating of arrays

   3 Want access to be as fast as possible
   4 Next want insertions to be as fast as possible

   5 Want to be able to insert new values at any point in the list
   NOTE: Key values will not change in value

  Answer:
   Use the STL map container because it meets goals 1,2,3 (logarithmic time based on number of nonzero entries in the Z direction)
    That is have a map for each Z value (x*y maps are needed)
   Goal 4 (log time insertion), 5 is satisfied

  Next Design goals:
   write out the sparse matrix in binary form which takes little space on disk and can be read in quickly.
   write out in an ascii form which can be loaded by matlab: see SparseLib++ from NIST for some code which does this
    Cant use their library because they dont allow new non-zero entries to be inserted into the sparse matrix




*/

typedef map<unsigned short, CGaussianDistribution> TypeGausDistMap;


typedef unsigned long ulong;
typedef unsigned short ushort;

class CSparse2DGausDistMatrix
{
public:
  // copy ctor
  CSparse2DGausDistMatrix(const CSparse2DGausDistMatrix& srcSparseMatrix)
  {
    nGausDistMeanVectorDim=srcSparseMatrix.nGausDistMeanVectorDim;
    nXMAX=srcSparseMatrix.nXMAX;
    nYMAX=srcSparseMatrix.nYMAX;

    //if (arrMap!=NULL)
    //   delete[] arrMap;

    arrMap=new TypeGausDistMap[nXMAX];

    // copy the src's map into this's map
    for (int nWhichMap=0; nWhichMap<nXMAX; nWhichMap++)
    {
      arrMap[nWhichMap]=srcSparseMatrix.arrMap[nWhichMap];
    }

  }

  // declare the null ctor to satisfy map<> class. Do not use directly
  CSparse2DGausDistMatrix()
  {
    nGausDistMeanVectorDim=0;
    nXMAX=0;
    nYMAX=0;
    arrMap=NULL;
  }

  CSparse2DGausDistMatrix(ushort nXDim, ushort nYDim, ushort nNewGausDistMeanVectorDim)
  {
    nGausDistMeanVectorDim=nNewGausDistMeanVectorDim;
    nXMAX=nXDim;
    nYMAX=nYDim;  // you could alter the yDim Size but not the XDim size
    arrMap=new TypeGausDistMap[nXMAX];
    for (int i=0; i<nXMAX; i++)
    {
      arrMap[i].clear();
    }

  }

  ~CSparse2DGausDistMatrix()
  {
    if (arrMap!=NULL)
      delete[] arrMap;
  }


  // increment values at specific location
  bool insertMeasure(const TypeVectorFloat& measure, ushort x, ushort y, float fFractionOfAMeasure=1.0)
  {
    bool bInserted=false;
    TypeGausDistMap::iterator it;
    if (x<=nXMAX-1)  // && (x>=0)
    {
      if ( (it=arrMap[x].find(y)) == arrMap[x].end())
      { // Init a new gaus dist
        CGaussianDistribution newGausDist(nGausDistMeanVectorDim);
        newGausDist.insertMeasure(measure,fFractionOfAMeasure);
        arrMap[x].insert(TypeGausDistMap::value_type(y,newGausDist));
      }
      else
      { // update the existing gaus dist
        arrMap[x][y].insertMeasure(measure,fFractionOfAMeasure);
      }

      if (y+1>nYMAX)
      {
        nYMAX=y+1;
      }

      bInserted=true;
    }

    return(bInserted);
  }

  // Fwd to the GaussianDistribution class.
  // retrieve values at specific location
  // If there is no gaussian dist for a given location, an empty gauss dist
  //  is not created and queried, instead the density zero is returned directly
  double density(const TypeVectorFloat& measure,ushort x, ushort y)
  {
    double dValue=0;
    TypeGausDistMap::iterator it;

    if (x<=nXMAX-1)  // && (x>=0)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        dValue=arrMap[x][y].density(measure);
      }
    }

    return dValue;
  }

  double extractMeanComponent(ushort x, ushort y, int nWhichComponent)
  {
    double dValue=0;
    TypeGausDistMap::iterator it;

    if (x<=nXMAX-1)  // && (x>=0)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        dValue=it->second.meanVector()[nWhichComponent]; // no component checking
      }
    }

    return dValue;
  }

  // returns the (nWhichComponent,nWhichComponent) diagonal entry of the Cov Matrix
  double extractVarianceComponent(ushort x, ushort y, int nWhichComponent)
  {
    double dValue=0;
    TypeGausDistMap::iterator it;

    if (x<=nXMAX-1)  // && (x>=0)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        dValue=it->second.covMatrix()[nWhichComponent][nWhichComponent]; // no component checking
      }
    }
    return dValue;
  }
#if 0
  double extractMeanComponentWithNegForNonExistant(ushort x, ushort y, int nWhichComponent)
  {
    double dValue=-1;
    TypeGausDistMap::iterator it;

    if (x<=nXMAX-1)  // && (x>=0)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        dValue=it->second.meanVector()[nWhichComponent]; // no component checking
      }
    }

    return dValue;
  }

  // returns the (nWhichComponent,nWhichComponent) diagonal entry of the Cov Matrix
  double extractVarianceComponentWithNegForNonExistant(ushort x, ushort y, int nWhichComponent)
  {
    double dValue=-1;
    TypeGausDistMap::iterator it;

    if (x<=nXMAX-1)  // && (x>=0)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        dValue=it->second.covMatrix()[nWhichComponent][nWhichComponent]; // no component checking
      }
    }
    return dValue;
  }
#endif

  bool setGausParams(ushort x, ushort y,
                     TypeVectorFloat& vectFloatMeasureMean,
                     TypeMatrixFloat& matrixFloatMeasureVariance,
                     float floatNewNumMeasuresEntered)
  {
    bool bSet=false;
    TypeGausDistMap::iterator it;
    if (x<=nXMAX-1)
    {
      if ( (it=arrMap[x].find(y)) == arrMap[x].end())
      { // Init a new gaus dist
        CGaussianDistribution newGausDist;
        newGausDist.init(vectFloatMeasureMean,matrixFloatMeasureVariance,floatNewNumMeasuresEntered);
        arrMap[x].insert(TypeGausDistMap::value_type(y,newGausDist));
      }
      else
      { // update the existing gaus dist
        arrMap[x][y].init(vectFloatMeasureMean,matrixFloatMeasureVariance,floatNewNumMeasuresEntered);
      }

      if (y+1>nYMAX)
      {
        nYMAX=y+1;
      }

      bSet=true;
    }


    return(bSet);
  }


  /* [[]] to be implemented another day. Requires removing the null ctor so that an appropriately sized
   // zeroDist can be created
  // retrieve the gaussian distribution at a specific location
  // if that location does not have a distribution then a const reference
  // to the internal zero dist is returned
  const CSparse2DGausDistMatrix& operator()(ushort x, ushort y)
    {
      TypeGausDistMap::iterator it;

      if (x<=nXMAX-1)  // && (x>=0)
  {
    if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
  return(arrMap[x][y]);
      }
  }
      else return(gausDistZeroDist);
    }
  */

private:
  ushort nXMAX, nYMAX;
  TypeGausDistMap* arrMap;
  // CGaussianDistribution gausDistZeroDist;
  ushort nGausDistMeanVectorDim;

  friend istream& operator>>(istream&,CSparse2DGausDistMatrix&);
  friend ostream& operator<<(ostream&,CSparse2DGausDistMatrix&);
  friend ostream& print(ostream&, CSparse2DGausDistMatrix&);

};

istream& operator>>(istream& is, CSparse2DGausDistMatrix& mat)
{
  // Read in from the stream as binary data
  is.read(&mat.nXMAX, sizeof(mat.nXMAX));
  is.read(&mat.nYMAX, sizeof(mat.nYMAX));
  is.read(&mat.nGausDistMeanVectorDim, sizeof(mat.nGausDistMeanVectorDim));

  // Dynamic 2D array allocation
  // allocate row arrays
  mat.arrMap=new TypeGausDistMap[mat.nXMAX];

  ushort nMapIndex;
  ushort x;
  ushort nKey;
  CGaussianDistribution newGausDist(mat.nGausDistMeanVectorDim);
  for (x=0; x<mat.nXMAX; x++)
  {
    ulong nMapSize;
    is.read(&nMapSize, sizeof(nMapSize));
    // iterate over the values in the map and read in the (key,value) pairs
    for (nMapIndex=0; nMapIndex<nMapSize; nMapIndex++)
    {
      is.read(&nKey, sizeof(nKey));

      // since all the contents are not saved, read in only what was saved from the Gaus Dist.
      // is.read(&newGausDist, sizeof(newGausDist));
      is >> newGausDist;

      mat.arrMap[x][nKey]=newGausDist;
    }
  }

  return is;
}

ostream& operator<<(ostream& os, CSparse2DGausDistMatrix& mat)
{
  // write out to the stream as binary data
  os.write((unsigned char *)&(mat.nXMAX), sizeof(mat.nXMAX));
  os.write((unsigned char *)&(mat.nYMAX), sizeof(mat.nYMAX));
  os.write((unsigned char *)&(mat.nGausDistMeanVectorDim), sizeof(mat.nGausDistMeanVectorDim));


  TypeGausDistMap::iterator it;
  ushort x;
  for (x=0; x<mat.nXMAX; x++)
  {
    ulong nMapSize=mat.arrMap[x].size();
    os.write((unsigned char *)&nMapSize, sizeof(nMapSize));
    // iterate over the values in the map and write out the (key,value) pairs
    for (it=mat.arrMap[x].begin(); it!=mat.arrMap[x].end(); it++)
    {
      os.write((unsigned char *)&(it->first), sizeof(it->first));

      //[[]] debug
#if 0
      if ( (x==50) && (it->first==67))
      {
        char cstrMsg[500];
        sprintf(cstrMsg, "neg var ");
        cout << cstrMsg << "\n";

      }
#endif

      // instead of just writing out the whole gaus dist, write out only what is needed
      // os.write((unsigned char *)&(it->second), sizeof(it->second));
      os << it->second;
    }

  }
  return os;
}


ostream& print(ostream& os, CSparse2DGausDistMatrix& mat)
{
  // write out to the stream as binary data
  os << "Dimensions: (x,y) = (" << mat.nXMAX << ", " << mat.nYMAX << ")  \n";
  os << "Dim of mean gaus vector = " << mat.nGausDistMeanVectorDim << "\n";

  TypeGausDistMap::iterator it;
  ushort x;
  //  for (x=0; x<mat.nXMAX; x++)
  for (x=0; x<3; x++) // truncate the huge dumps
  { ulong nMapSize=mat.arrMap[x].size();
    os << "  (" << x << ") has " << nMapSize << "  {y, <GausDist>} elements:";

    // iterate over the values in the map and write out the (key,value) pairs
    for (it=mat.arrMap[x].begin(); it!=mat.arrMap[x].end(); it++)
    {
      os << " {" << it->first << ", ";
      print(os,it->second);
      os << " } ";
    }
    os << "\n";

  }
  return os;
}



#endif



