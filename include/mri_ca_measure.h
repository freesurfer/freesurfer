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


#ifndef Measure_h
#define Measure_h

#include <stdlib.h>
#include "mri_ca_configFile.h"
#include "mri_ca_util.h"
using namespace std ;

class CMeasure;
typedef vector<class CMeasure> TypeVectorMeasures;

class CMeasure
{

public:
  // Default ctor, dtor, copy and assign OK

  CMeasure()
  {}

  void get(CConfigFile& configFile, string strSectionName)
  {
    if (strstr(strSectionName.c_str(), "T1")!=NULL)
      measureType=T1;
    else if (strstr(strSectionName.c_str(), "T2")!=NULL)
      measureType=T2;
    else if (strstr(strSectionName.c_str(), "PD")!=NULL)
      measureType=PD;
    else if (strstr(strSectionName.c_str(), "EEG")!=NULL)
      measureType=EEG;
    else if (strstr(strSectionName.c_str(), "MEG")!=NULL)
      measureType=MEG;
    else if (strstr(strSectionName.c_str(), "fMRI")!=NULL)
      measureType=fMRI;
    configFile.get(vectDoubleIntensityScales,  strSectionName, "IntensityScales");
    configFile.get(vectDoubleGradientScales,  strSectionName, "GradientScales");
    configFile.get(vectDoubleCurvatureScales,  strSectionName, "CurvatureScales");

  }


  void write(CConfigFile& configFile, bool bEchoToStdOut=false)
  {
    string strSectionName="MEASURE: ";
    switch (measureType)
    {
    case T1:
      strSectionName.append("T1");
      break;
    case T2:
      strSectionName.append("T2");
      break;
    case PD:
      strSectionName.append("PD");
      break;
    case EEG:
      strSectionName.append("EEG");
      break;
    case MEG:
      strSectionName.append("MEG");
      break;
    case fMRI:
      strSectionName.append("fMRI");
      break;
    }

    configFile.writeSection(strSectionName,bEchoToStdOut);
    configFile.write(vectDoubleIntensityScales,   "IntensityScales",bEchoToStdOut);
    configFile.write(vectDoubleGradientScales,   "GradientScales",bEchoToStdOut);
    configFile.write(vectDoubleCurvatureScales,  "CurvatureScales",bEchoToStdOut);

#if 0
    configFile.write(vectDoubleIntensityScales,  strSectionName, "IntensityScales",bEchoToStdOut);
    configFile.write(vectDoubleGradientScales,  strSectionName, "GradientScales",bEchoToStdOut);
    configFile.write(vectDoubleCurvatureScales,  strSectionName, "CurvatureScales",bEchoToStdOut);
#endif
  }



public:
  // use a class scoped enumeration type
  enum MeasureType { T1, T2, PD, EEG, MEG, fMRI};


  MeasureType      measureType;
  TypeVectorDouble vectDoubleIntensityScales;
  TypeVectorDouble vectDoubleGradientScales;
  TypeVectorDouble vectDoubleCurvatureScales;
};


class CMeasures
{
public:

  void get(CConfigFile& configFile)
  {
    TypeVectorString vectStrMeasures;
    configFile.find(vectStrMeasures, "MEASURE:");

    int nMaxMeasureNameIndex=vectStrMeasures.size()-1;
    for (int nWhichMeasureNameIndex=0; nWhichMeasureNameIndex<=nMaxMeasureNameIndex; nWhichMeasureNameIndex++)
    {
      CMeasure measure;
      measure.get(configFile, vectStrMeasures[nWhichMeasureNameIndex]);
      measures.push_back(measure); // add the measure read in to the list of measures
    }


  }

  void write(CConfigFile& configFile, bool bEchoToStdOut=false)
  {
    int nMaxMeasureIndex=measures.size()-1;
    for (int nWhichMeasureIndex=0; nWhichMeasureIndex<=nMaxMeasureIndex; nWhichMeasureIndex++)
    {
      measures[nWhichMeasureIndex].write(configFile,bEchoToStdOut);
    }
  }


  void clear()
  {
    measures.clear();

  }

public:
  TypeVectorMeasures measures;

};


#endif
