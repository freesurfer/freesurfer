//
// CostFunctions is a class that makes available many different
// cost functions for images
//
// written by Martin Reuter
// March 20th ,2009
//

#ifndef CostFunctions_H
#define CostFunctions_H

#ifdef __cplusplus
extern "C" {
#endif
#include "mri.h"
#ifdef __cplusplus
}
#endif

//#include <utility>
//#include <string>
//#include <vector>
#include <iostream>

class CostFunctions
{
    public:
    double mean(MRI *i);
    double var(MRI *i);
    double sdev(MRI *i) {return sqrt(var(i));};
    double leastSquares(MRI * i1, MRI * i2 = NULL);
    double normalizedCorrelation(MRI * i1, MRI * i2);
    double mutualInformation(MRI * i1, MRI * i2 = NULL);
    double normalizedMutualInformation(MRI * i1, MRI * i2 = NULL);
    double woods(MRI * i1, MRI * i2 = NULL);
    double correlationRatio(MRI * i1, MRI * i2 = NULL);

};

class MRIiterator
{
   public:
	MRIiterator(MRI * i);
	
	void begin();
	bool isEnd();
	MRIiterator& operator++(int);
	float operator*();
	

   protected:
      
        float fromUCHAR(void);
	float fromSHORT(void);
	float fromINT(void);
	float fromLONG(void);
	float fromFLOAT(void);
	
	MRIiterator& opincchunk(int);
	MRIiterator& opincnochunk(int);
   
   	MRI * img;
	unsigned char * pos;
	unsigned char * end;
	float (MRIiterator::*getVal)(void);
	MRIiterator& (MRIiterator::*opinc)(int);
	int x,y,z;
	int bytes_per_voxel;
};

inline MRIiterator::MRIiterator(MRI * i):img(i)
{
   // set current pos to begin
   // and initialize end pointer
   begin();
   
    switch (img->type)
    {
    case MRI_UCHAR:
      getVal = &MRIiterator::fromUCHAR;
      bytes_per_voxel = sizeof(unsigned char);
      break;
    case MRI_SHORT:
      getVal = &MRIiterator::fromSHORT;
      bytes_per_voxel = sizeof(short);
      break;
    case MRI_INT:
      getVal = &MRIiterator::fromINT;
      bytes_per_voxel = sizeof(int);
      break;
    case MRI_LONG:
      getVal = &MRIiterator::fromLONG;
      bytes_per_voxel = sizeof(long);
      break;
    case MRI_FLOAT:
      getVal = &MRIiterator::fromFLOAT;
      bytes_per_voxel = sizeof(float);
      break;
    default:
       std::cerr << "MRIiterator: Type not supported: " << img->type << std::endl;
       exit(1);
    }
   
   
}

inline void MRIiterator::begin()
// set pos to first element
{
   if (img->ischunked)
   {
      pos = (unsigned char*) img->chunk;
      end = (unsigned char*) img->chunk + img->bytes_total;
      opinc = &MRIiterator::opincchunk;
   }
   else
   {
      x = 0; y=0; z=0;
      pos = (unsigned char*) img->slices[0][0];
      end = NULL;
      opinc = &MRIiterator::opincnochunk;
   }
}

inline bool MRIiterator::isEnd()
{
    //    if(pos > end && end != 0)
//	{
//	   std::cerr << "MRIiterator::isEnd outside data???" << std::endl;
//	   exit(1);
//	}
	return (pos == end);
}

inline MRIiterator& MRIiterator::operator++(int i)
{
   return (this->*opinc)(i);
}

inline MRIiterator& MRIiterator::opincchunk(int)
{
//   if (pos < end)
      pos += img->bytes_per_vox;
      return *this;
}

inline MRIiterator& MRIiterator::opincnochunk(int)
{
   x++;
   if (x == img->width)
   {
      x = 0; y++;
      if (y == img->height)
      {
         y=0; z++;
	 if (z== img->depth)
	 {
	    z=0;
	    pos = NULL;
	    return *this;
	 }
      }
      pos = (unsigned char*) img->slices[z][y];
   }
   else pos += bytes_per_voxel;
   return *this;
}

inline float MRIiterator::fromUCHAR()
{return ((float) *(unsigned char *)pos);}

inline float MRIiterator::fromSHORT()
{return ((float) *(short *)pos);}

inline float MRIiterator::fromINT()
{return ((float) *(int *)pos);}

inline float MRIiterator::fromLONG()
{return ((float) *(long *)pos);}

inline float MRIiterator::fromFLOAT()
{return ((float) *(float *)pos);}

inline float MRIiterator::operator*()
{
//   if (pos < end && pos >= img->chunk)
		return (this->*getVal)();
}

// example:
// MRIiterator it(mri);
// for (it.begin(); !it.isEnd(); it++) 
// {
//    std::cout << *it << std::endl;
////    *it = 0;
// }


#endif
