/**
 * @file  mri_ca_sparse3DMatrix.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.3 $
 *
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


#ifndef sparse3DMatrix_h
#define sparse3DMatrix_h

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

typedef map<unsigned short, unsigned long> mapType;
typedef mapType* mapTypePtr;

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
class CSparse3DMatrix
{
public:

  /*
  CSparse3DMatrix()
    {
      nXMAX=0;
      nYMAX=0;
      nZMAX=0;
    }
  */

  CSparse3DMatrix(ushort nXDim, ushort nYDim, ushort nZDim)
  {
    nXMAX=nXDim;
    nYMAX=nYDim;
    nZMAX=nZDim;  // only the first two dimensions must remain fixed
    //arrMap=new mapType[nXMAX][nYMAX]; // wont work b/c it returns mapType (*)[nYMax]

    // Dynamic 2D array allocation
    // allocate row arrays
    arrMap=new mapTypePtr[nXMAX]; // returns mapTypePtr[nXMAX]
    for (ushort x=0; x<nXMAX; x++)
    {
      arrMap[x]=new mapType[nYMAX];
    }

  }

  ~CSparse3DMatrix()
  {
    for (ushort x=0; x<nXMAX; x++)
    {
      delete[] arrMap[x];
    }
    delete[] arrMap;
  }


  // increment values at specific location
  ulong increment(ushort x, ushort y, ushort z)
  {
    mapType::iterator it;
    ulong nValue=0;
    if ((x<=nXMAX-1) && (y<=nYMAX-1))  // && (y>=0) && (x>=0)
    {
      if ( (it=arrMap[x][y].find(z)) != arrMap[x][y].end())
      {
        nValue=++arrMap[x][y][z];
      }
      else
      {
        nValue=arrMap[x][y][z]=1;
      }


      if (z+1>nZMAX)
      {
        nZMAX=z+1;
      }
    }
    return nValue;

  }


  // retrieve values at specific location
  ulong operator()(ushort x, ushort y, ushort z)
  {
    ulong nValue=0;
    mapType::iterator it;

    if ((x<=nXMAX-1) && (y<=nYMAX-1))  // && (y>=0) && (x>=0)
    {
      if ( (it=arrMap[x][y].find(z)) != arrMap[x][y].end())
      {
        nValue=arrMap[x][y][z];
      }
    }
    return nValue;
  }




private:
  ushort nXMAX, nYMAX, nZMAX;
  mapType** arrMap;

  friend istream& operator>>(istream&,CSparse3DMatrix&);
  friend ostream& operator<<(ostream&,CSparse3DMatrix&);
  friend ostream& print(ostream&, CSparse3DMatrix&);

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






istream& operator>>(istream& is, CSparse3DMatrix& mat)
{
  // Read in from the stream as binary data
  is.read(&mat.nXMAX, sizeof(mat.nXMAX));
  is.read(&mat.nYMAX, sizeof(mat.nYMAX));
  is.read(&mat.nZMAX, sizeof(mat.nZMAX));

  // Dynamic 2D array allocation
  // allocate row arrays
  mat.arrMap=new mapTypePtr[mat.nXMAX]; // returns mapTypePtr[nXMAX]
  for (ushort x=0; x<mat.nXMAX; x++)
  {
    mat.arrMap[x]=new mapType[mat.nYMAX];
  }


  ushort nMapIndex;
  ushort x,y;
  ushort nKey;
  ulong  nValue;
  for (x=0; x<mat.nXMAX; x++)
    for (y=0; y<mat.nYMAX; y++)
    {
      ulong nMapSize;
      is.read(&nMapSize, sizeof(nMapSize));
      // iterate over the values in the map and read in the (key,value) pairs
      for (nMapIndex=0; nMapIndex<nMapSize; nMapIndex++)
      {
        is.read(&nKey, sizeof(nKey));
        is.read(&nValue, sizeof(nValue));
        mat.arrMap[x][y][nKey]=nValue;
      }
    }

  return is;
}

ostream& operator<<(ostream& os, CSparse3DMatrix& mat)
{
  // write out to the stream as binary data
  os.write((unsigned char *) &(mat.nXMAX), sizeof(mat.nXMAX));
  os.write((unsigned char *)&(mat.nYMAX), sizeof(mat.nYMAX));
  os.write((unsigned char *)&(mat.nZMAX), sizeof(mat.nZMAX));

  mapType::iterator it;
  ushort x,y;
  for (x=0; x<mat.nXMAX; x++)
    for (y=0; y<mat.nYMAX; y++)
    {
      ulong nMapSize=mat.arrMap[x][y].size();
      os.write((unsigned char *)&nMapSize, sizeof(nMapSize));
      // iterate over the values in the map and write out the (key,value) pairs
      for (it=mat.arrMap[x][y].begin(); it!=mat.arrMap[x][y].end(); it++)
      {
        os.write((unsigned char *)&(it->first), sizeof(it->first));
        os.write((unsigned char *)&(it->second), sizeof(it->second));
      }

    }
  return os;
}

ostream& print(ostream& os, CSparse3DMatrix& mat)
{
  // write out to the stream as binary data
  os << "Dimensions: (x,y,z) = (" << mat.nXMAX << ", " << mat.nYMAX << ", " << mat.nZMAX << ")  \n";


  mapType::iterator it;
  ushort x,y;
  for (x=0; x<mat.nXMAX; x++)
    for (y=0; y<mat.nYMAX; y++)
    {
      ulong nMapSize=mat.arrMap[x][y].size();
      os << "  (" << x <<", " << y  << ") has " << nMapSize << "  {z, val} elements:";

      // iterate over the values in the map and write out the (key,value) pairs
      for (it=mat.arrMap[x][y].begin(); it!=mat.arrMap[x][y].end(); it++)
      {
        os << " {" << it->first << ", " << it->second << "} ";
      }
      os << "\n";

    }
  return os;
}



/*
// a PriorDist holds the probability of a label at a given location
class CPriorDist
{

  probability(x,y,z)



 public:



}


// now define a collection (vector) of prior distributions
class CPriorDists
{

}

*/


#endif


