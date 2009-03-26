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
   
   	MRI * img;
	unsigned char * pos;
	unsigned char * end;
	float (MRIiterator::*getVal)(void);
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
      break;
    case MRI_SHORT:
      getVal = &MRIiterator::fromSHORT;
      break;
    case MRI_INT:
      getVal = &MRIiterator::fromINT;
      break;
    case MRI_LONG:
      getVal = &MRIiterator::fromLONG;
      break;
    case MRI_FLOAT:
      getVal = &MRIiterator::fromFLOAT;
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
      end = (unsigned char*)img->chunk + img->bytes_total;
   }
   else
   {
      std::cerr << "MRIiterator:begin(): Need chunck MRI, set environment: FS_USE_MRI_CHUNK" << std::endl;
      exit(1);
   }
}

inline bool MRIiterator::isEnd()
{
        if(pos > end)
	{
	   std::cerr << "MRIiterator::isEnd outside data???" << std::endl;
	   exit(1);
	}
	return (pos == end);
}

inline MRIiterator &MRIiterator::operator++(int)
{
//   if (pos < end)
//   {
      pos += img->bytes_per_vox;
//   }
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
//   {
		return (this->*getVal)();
//   }
}

// example:
// MRIiterator it<float>(mri);
// for (it.begin(); !it.isEnd(); it++) 
// {
//    std::cout << *it.val() << std::endl;
//    it.val() = 0;
// }


#endif
