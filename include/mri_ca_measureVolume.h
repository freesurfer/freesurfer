#ifndef MeasureVolume_h
#define MeasureVolume_h

#include <string>
extern "C" 
{
#include "mri.h" 
}

using namespace std;




class CMeasureVolume
{
 public:
    enum MeasureType { T1, T2, PD, EEG, MEG, fMRI, Unknown};


    MeasureType measureType;
    MRI* pVolume;
    string strMeasureFileDir;    
};


#endif
