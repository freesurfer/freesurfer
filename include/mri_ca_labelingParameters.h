#ifndef LabelingParameters_h
#define LabelingParameters_h

using namespace std ;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "mri_ca_configFile.h"
#include "mri_ca_measure.h"
#include "mri_ca_util.h"
#include "mri_ca_trainingParameters.h"


class CLabelingParameters 
{

 public: 
  void read()
     { // init the data members to be empty
       strStatisticsVolumeHeaderPath.erase();
       arr2DstrSubjectsToLabel.clear();

       // get the values from the config file
       labelingConfigFile.get(strLabelLUT, "GeneralLabelingParameters", "LabelLUT");
       labelingConfigFile.get(strLabeledVolumeSourceType, "GeneralLabelingParameters", "LabeledVolumeSourceType");
       if (strLabeledVolumeSourceType=="CMA")
	 {
          enumLabeledVolumeSourceType=enumCMA;
	 }
       else if (strLabeledVolumeSourceType=="NMRAutoFixer")
	 {enumLabeledVolumeSourceType=enumNMRAutoFixer;
	 }
       else if (strLabeledVolumeSourceType=="Killiany")
	 {enumLabeledVolumeSourceType=enumKilliany;
	 }
       else
	 {enumLabeledVolumeSourceType=enumUnknown;
	 };   
       labelingConfigFile.get(strStatisticsVolumeHeaderPath, "GeneralLabelingParameters", "StatisticsVolumeHeaderPath");
       labelingConfigFile.get(nNumberOfBlocks, "GeneralLabelingParameters", "nNumberOfBlocks");
       if (nNumberOfBlocks<=0) 
	 {
           cout << "Error: nNumberOfBlocks must be a positive integer\n";
           exit(1);
	 }

       labelingConfigFile.get(arr2DstrSubjectsToLabel, "SubjectsToLabel", "");
     }

  void write(bool bEchoToStdOut=false)
    { 
      labelingConfigFile.writeSection("GeneralLabelingParameters", bEchoToStdOut);
        labelingConfigFile.write(strLabelLUT,  "LabelLUT",bEchoToStdOut);
        labelingConfigFile.write(strLabeledVolumeSourceType,  "LabeledVolumeSourceType",bEchoToStdOut);
        labelingConfigFile.write(strStatisticsVolumeHeaderPath,  "StatisticsVolumeHeaderPath",bEchoToStdOut);
        labelingConfigFile.write(nNumberOfBlocks,  "nNumberOfBlocks",bEchoToStdOut);
      labelingConfigFile.writeSection("SubjectsToLabel", bEchoToStdOut);
        labelingConfigFile.write(arr2DstrSubjectsToLabel, "",bEchoToStdOut);
    }

   // default dtor, copy , assign OK
  CLabelingParameters(string strFilepath, int nProt=0)
    {   
      labelingConfigFile.init(strFilepath, nProt);
       strStatisticsVolumeHeaderPath.erase();
       arr2DstrSubjectsToLabel.clear();
    }


  public:  

   LabeledVolumeSourceType enumLabeledVolumeSourceType;
   string strLabeledVolumeSourceType;   
   string strStatisticsVolumeHeaderPath;
   TypeMatrixString arr2DstrSubjectsToLabel;
   string strLabelLUT;
   int nNumberOfBlocks;
   //--
   CConfigFile labelingConfigFile;

 private:
   CLabelingParameters() {};  // disable default ctor
};

#endif
