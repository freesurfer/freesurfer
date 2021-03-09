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


#ifndef TrainingParameters_h
#define TrainingParameters_h

using namespace std ;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "mri_ca_configFile.h"
#include "mri_ca_measure.h"
#include "mri_ca_util.h"



enum LabeledVolumeSourceType { enumCMA=0, enumNMRAutoFixer, enumKilliany, enumUnknown};

LabeledVolumeSourceType stringToLabeledVolumeSourceType(string strLabeledVolumeSourceType)
{
  LabeledVolumeSourceType newEnumLabeledVolumeSourceType;
  if (strLabeledVolumeSourceType=="CMA")
  {
    newEnumLabeledVolumeSourceType=enumCMA;
  }
  else if (strLabeledVolumeSourceType=="NMRAutoFixer")
  {
    newEnumLabeledVolumeSourceType=enumNMRAutoFixer;
  }
  else if (strLabeledVolumeSourceType=="Killiany")
  {
    newEnumLabeledVolumeSourceType=enumKilliany;
  }
  else
  {
    newEnumLabeledVolumeSourceType=enumUnknown;
  };
  return(newEnumLabeledVolumeSourceType);
}

class CTrainingParameters
{

public:



  void read()
  { // init the data members to be empty
    strStatisticsVolumeHeaderPath.erase();
    arr2DstrSubjectsToTrain.clear();
    measures.clear();

    // get the values from the config file
    trainingConfigFile.get(strLabelLUT, "GeneralTrainingParameters", "LabelLUT");
    trainingConfigFile.get(strLabeledVolumeSourceType, "GeneralTrainingParameters", "LabeledVolumeSourceType");
    enumLabeledVolumeSourceType=stringToLabeledVolumeSourceType(strLabeledVolumeSourceType);
    trainingConfigFile.get(strStatisticsVolumeHeaderPath, "GeneralTrainingParameters", "StatisticsVolumeHeaderPath");
    trainingConfigFile.get(arr2DstrSubjectsToTrain, "SubjectsToTrain", "");
    measures.get(trainingConfigFile);


    // test(); // call this to test the CConfigFile class
  }

  void write(bool bEchoToStdOut=false)
  {
    trainingConfigFile.writeSection("GeneralTrainingParameters",bEchoToStdOut);
    trainingConfigFile.write(strLabelLUT, "LabelLUT",bEchoToStdOut);
    trainingConfigFile.write(strLabeledVolumeSourceType,  "LabeledVolumeSourceType",bEchoToStdOut);
    trainingConfigFile.write(strStatisticsVolumeHeaderPath, "StatisticsVolumeHeaderPath",bEchoToStdOut);
    trainingConfigFile.writeSection("SubjectsToTrain",bEchoToStdOut);
    trainingConfigFile.write(arr2DstrSubjectsToTrain, "",bEchoToStdOut);
    measures.write(trainingConfigFile,bEchoToStdOut);

#if 0
    trainingConfigFile.write(strLabelLUT, "GeneralTrainingParameters", "LabelLUT",bEchoToStdOut);
    trainingConfigFile.write(strLabeledVolumeSourceType, "GeneralTrainingParameters", "LabeledVolumeSourceType",bEchoToStdOut);
    trainingConfigFile.write(strStatisticsVolumeHeaderPath,"GeneralTrainingParameters", "StatisticsVolumeHeaderPath",bEchoToStdOut);
    trainingConfigFile.write(arr2DstrSubjectsToTrain, "SubjectsToTrain", "",bEchoToStdOut);
    measures.write(trainingConfigFile,bEchoToStdOut);
#endif


  }

  // default dtor, copy , assign OK
  CTrainingParameters(string strFilepath, int nProt=0)
  {
    trainingConfigFile.init(strFilepath, nProt);
    strStatisticsVolumeHeaderPath.erase();
    arr2DstrSubjectsToTrain.clear();
    measures.clear();
  }


public:  // changed to public access so that these could be accessed in train.cxx

  LabeledVolumeSourceType enumLabeledVolumeSourceType;
  string strLabeledVolumeSourceType;
  string strStatisticsVolumeHeaderPath;
  TypeMatrixString arr2DstrSubjectsToTrain;
  CMeasures measures;
  string strLabelLUT;
  //--
  CConfigFile trainingConfigFile;

private:
  CTrainingParameters()
  {};  // disable default ctor
};

#endif





