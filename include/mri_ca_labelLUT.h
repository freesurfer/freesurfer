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


#ifndef labelLUT_h
#define labelLUT_h

#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "mri_ca_configFile.h"
#include "mri_ca_util.h"

using namespace std;

typedef map<int, string> TypeMapLUT;




class CLabelLUT
{
public:
  // default dtor, copy , assign OK
  CLabelLUT(string strFilepath, int nProt=0)
  {
    labelLUTConfigFile.init(strFilepath, nProt);
    mapLut.clear();
    read();
  }

  void read()
  { // init the data members to be empty
    mapLut.clear();

    TypeMatrixString arr2DLUTStrings;
    // get the values from the config file
    labelLUTConfigFile.get(arr2DLUTStrings, "LabelLUT", "");
    for (int i=0; i<(int)arr2DLUTStrings.size(); i++)
    {
      if (arr2DLUTStrings[i].size()>=2)
      {
        int nKey=atoi(arr2DLUTStrings[i][0].c_str());
        string strLabel=arr2DLUTStrings[i][1];
        mapLut[nKey]=strLabel;
      }
    }
  }

  void write(bool bEchoToStdOut=false)
  {
    TypeMatrixString arr2DLUTStrings;

    TypeMapLUT::iterator it;
    for (it=mapLut.begin(); it!=mapLut.end(); it++)
    {
      int nKey=it->first;
      string strLabel=it->second;
      char cstrKey[500];
      sprintf(cstrKey,"%d",nKey);
      string strKey=cstrKey;
      TypeVectorString vectorString;
      vectorString.push_back(strKey);
      vectorString.push_back(strLabel);
      arr2DLUTStrings.push_back(vectorString);
    }

    // getwrite  the values to the config file
    labelLUTConfigFile.writeSection("LabelLUT",bEchoToStdOut);
    labelLUTConfigFile.write(arr2DLUTStrings, "",bEchoToStdOut);
  }


  int size()
  {
    return(mapLut.size());
  }

  string label(int nKey)
  {
    return(mapLut[nKey]);
  }

  TypeMapLUT& getMapLUT()
  {
    return(mapLut);
  }


private:

  TypeMapLUT mapLut;

  //friend istream& operator>>(istream&,CLabelLUT&);
  //friend ostream& operator<<(ostream&,CLabelLUT&);

  CConfigFile labelLUTConfigFile;

};

/*
 istream& operator>>(istream& is, CLabelLUT& labelLUT)
 {
   labelLUT.mapLut.clear();
   ulong nMapSize;
   is >> nMapSize;
   for (ulong nMapIndex=0; nMapIndex<nMapSize; nMapIndex++)
     {  int nKey;
 string strLabel;
 is >> nKey;
 is >> strLabel;
 labelLUT.mapLut[nKey]=strLabel;
      }

   return is;
 }


 ostream& operator<<(ostream& os, CLabelLUT& labelLUT)
 {

   ulong uMapSize= labelLUT.mapLut.size();
   os << uMapSize;

   TypeMapLUT::iterator it;
   for (it=labelLUT.mapLut.begin(); it!=labelLUT.mapLut.end(); it++)
     {  int nKey = it->first;
 string strLabel= it->second;
 os << nKey << " " << strLabel << "\n";
      }

   return os;
 }
*/

#endif






















