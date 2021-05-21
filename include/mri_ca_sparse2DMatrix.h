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


#ifndef sparse2DMatrix_h
#define sparse2DMatrix_h


#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

using namespace std;


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

typedef map<unsigned short, float> TypeFloatMap;
typedef TypeFloatMap* TypeFloatMapPtr;

typedef unsigned long ulong;
typedef unsigned short ushort;
/* STL map methods
map<key type, value type>
mymapType::iterator it;
insertion:
  myMap[a key]=a value;
  insert(mymapType::value_type(a key, a value)

accessing values:
 if myMap.find(a key)!=myMap.end()  then
    (myMap.find(a key))->fist=the key
    (myMap.find(a key))->second=the value
Info:
 size()
 empty()
 clear()
*/
class CSparse2DMatrix
{
public:
  // copy ctor
  CSparse2DMatrix(const CSparse2DMatrix& srcSparseMatrix)
  {
    nXMAX=srcSparseMatrix.nXMAX;
    nYMAX=srcSparseMatrix.nYMAX;
    //if (arrMap!=NULL)
    //   delete[] arrMap;

    arrMap=new TypeFloatMap[nXMAX];

    // copy the src's map into this's map
    for (int nWhichMap=0; nWhichMap<nXMAX; nWhichMap++)
    {
      arrMap[nWhichMap]=srcSparseMatrix.arrMap[nWhichMap];
    }

  }

  // declare the null ctor to satisfy map<> class. Do not use directly
  CSparse2DMatrix()
  {
    nXMAX=0;
    nYMAX=0;
    arrMap=NULL;
  }

  CSparse2DMatrix(ushort nXDim, ushort nYDim)
  {
    nXMAX=nXDim;
    nYMAX=nYDim;  // you could alter the yDim Size but not the XDim size
    arrMap=new TypeFloatMap[nXMAX];
    for (int i=0; i<nXMAX; i++)
    {
      arrMap[i].clear();
    }

  }

  ~CSparse2DMatrix()
  {
    if (arrMap!=NULL)
      delete[] arrMap;
  }


  // increment values at specific location
  float increment(ushort x, ushort y, float fIncrementAmount=1.0)
  {
    TypeFloatMap::iterator it;
    float fValue=0;
    float fOldValue;
    if (x<=nXMAX-1)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        fOldValue=arrMap[x][y];
        fValue=arrMap[x][y]=fOldValue+fIncrementAmount;
      }
      else
      {
        fValue=arrMap[x][y]=fIncrementAmount;
      }


      if (y+1>nYMAX)
      {
        nYMAX=y+1;
      }
    }
    return fValue;

  }

  void setCount(ushort x, ushort y, float floatNewCount)
  {
    if (x<=nXMAX-1)
    {
      arrMap[x][y]=floatNewCount;

      if (y+1>nYMAX)
      {
        nYMAX=y+1;
      }
    }
  }

  // retrieve values at specific location
  float operator()(ushort x, ushort y)
  {
    float fValue=0;
    TypeFloatMap::iterator it;

    if (x<=nXMAX-1)
    {
      if ( (it=arrMap[x].find(y)) != arrMap[x].end())
      {
        fValue=arrMap[x][y];
      }
    }
    return fValue;
  }




private:
  ushort nXMAX, nYMAX;
  TypeFloatMap* arrMap;

  friend istream& operator>>(istream&,CSparse2DMatrix&);
  friend ostream& operator<<(ostream&,CSparse2DMatrix&);
  friend ostream& print(ostream&, CSparse2DMatrix&);

};


/*
on Linux / intel 686:
Size of pointer=4

Size of short=2
Size of u-short=2
Size of int=4
Size of unsigned int=4
Size of long int=4
Size of float=4
Size of double =8
 */






istream& operator>>(istream& is, CSparse2DMatrix& mat)
{
  // Read in from the stream as binary data
  is.read(&mat.nXMAX, sizeof(mat.nXMAX));
  is.read(&mat.nYMAX, sizeof(mat.nYMAX));


  // Dynamic 2D array allocation
  // allocate row arrays
  mat.arrMap=new TypeFloatMap[mat.nXMAX]; // returns TypeFloatMapPtr[nXMAX]

  ushort nMapIndex;
  ushort x;
  ushort nKey;
  float  fValue;
  for (x=0; x<mat.nXMAX; x++)
  {
    ulong nMapSize;
    is.read(&nMapSize, sizeof(nMapSize));
    // iterate over the values in the map and read in the (key,value) pairs
    for (nMapIndex=0; nMapIndex<nMapSize; nMapIndex++)
    {
      is.read(&nKey, sizeof(nKey));
      is.read(&fValue, sizeof(fValue));
      mat.arrMap[x][nKey]=fValue;
    }
  }

  return is;
}

ostream& operator<<(ostream& os, CSparse2DMatrix& mat)
{
  // write out to the stream as binary data
  os.write((unsigned char *) &(mat.nXMAX), sizeof(mat.nXMAX));
  os.write((unsigned char *)&(mat.nYMAX), sizeof(mat.nYMAX));

  TypeFloatMap::iterator it;
  ushort x;
  for (x=0; x<mat.nXMAX; x++)
  {
    ulong nMapSize=mat.arrMap[x].size();
    os.write((unsigned char *)&nMapSize, sizeof(nMapSize));
    // iterate over the values in the map and write out the (key,value) pairs
    for (it=mat.arrMap[x].begin(); it!=mat.arrMap[x].end(); it++)
    {
      os.write((unsigned char *)&(it->first), sizeof(it->first));
      os.write((unsigned char *)&(it->second), sizeof(it->second));
    }

  }
  return os;
}

ostream& print(ostream& os, CSparse2DMatrix& mat)
{
  // write out to the stream as binary data
  os << "Dimensions: (x,y) = (" << mat.nXMAX << ", " << mat.nYMAX << ")  \n";


  TypeFloatMap::iterator it;
  ushort x;
  for (x=0; x<mat.nXMAX; x++)
  {
    ulong nMapSize=mat.arrMap[x].size();
    os << "  (" << x << ") has " << nMapSize << "  {y, val} elements:";

    // iterate over the values in the map and write out the (key,value) pairs
    for (it=mat.arrMap[x].begin(); it!=mat.arrMap[x].end(); it++)
    {
      os << " {" << it->first << ", " << it->second << "} ";
    }
    os << "\n";

  }
  return os;
}



#endif


