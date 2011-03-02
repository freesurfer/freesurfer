/**
 * @file  mri_ca_labelMapper.h
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


#ifndef LabelMapper_h
#define LabelMapper_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>

#include "mri_ca_configFile.h"
#include "mri_ca_util.h"
#include "mri_ca_sparse2DMatrix.h"
#include "mri_ca_labelLUT.h"
#include "mri_ca_statisticsVolumeHeader.h"

extern "C"
{
#include "mri.h"
}


using namespace std;

typedef map<int, CSparse2DMatrix> TypeMapSparse2DMatrix;
typedef vector<CSparse2DMatrix> TypeVectorSparse2DMatrix;


enum enumVolVizMethod { eThresholdedLabelProbabilities=0, eTranslucentLabelProbabilities};

class CLabel2PriorMapper
{



public:
  CLabel2PriorMapper(CLabelLUT& theLabelLUT, CStatisticsVolumeHeader& newStatisticsVolumeHeader) : labelLUT(theLabelLUT), statisticsVolumeHeader(newStatisticsVolumeHeader)
  {};




  TypeMapSparse2DMatrix& loadSlicePriors(int nSliceNumber, int nWidth, int nHeight)
  {
    priorSlices.clear();
    nCurrSliceNum=nSliceNumber;
    // iterate over all labels in labelLUT
    TypeMapLUT::iterator it;

    for (it=labelLUT.getMapLUT().begin(); it!=labelLUT.getMapLUT().end(); it++)
    {
      // form the filename
      char strFilename[500];
      int nWhichLabel=it->first;
      sprintf(strFilename,"%s/Prior-Slc%03d-Lab%03d.spmat",statisticsVolumeHeader.strStatisticsVolumeDirPath.c_str(),nCurrSliceNum,nWhichLabel);

      // load the sparse matrix
      //  create an initially zero prior and then use >> to update  it if the file exists
      CSparse2DMatrix zeroPriorProb(nWidth,nHeight);
      //priorSlices[nWhichLabel]=zeroPriorProb;
      priorSlices.insert(TypeMapSparse2DMatrix::value_type(nWhichLabel, zeroPriorProb));
      fstream ifs;
      ifs.open(strFilename,ios::in | ios::nocreate);
      if (ifs.good())
      {
        ifs >> priorSlices[nWhichLabel];
        ifs.close();
      }
    }
    return(priorSlices);
  }

  // -----------------------------For Baysian Classifier------------------------------------------------------------------------------
  bool loadLabelProbabilityVolume(int nWhichLabel, bool bDebug=false)
  {
    bool bLoadedVolumeOK=true;
    vectorPriorSlices.clear();

    //  create an initially zero prior and then use >> to update  it if the file exists
    CSparse2DMatrix zeroPriorProb(statisticsVolumeHeader.nXDIM,statisticsVolumeHeader.nYDIM);

    for (int nWhichSlice=0; nWhichSlice<statisticsVolumeHeader.nNumSlices; nWhichSlice++)
    {
      // form the filename for this slice
      char strFilename[500];
      sprintf(strFilename,"%s/Prior-Slc%03d-Lab%03d.spmat",statisticsVolumeHeader.strStatisticsVolumeDirPath.c_str(),nWhichSlice,nWhichLabel);
      if (bDebug==false)
      {
        cout << "    Loading " << strFilename << "\n" << flush;
      }
      else
      {
        cout << "." << flush;
        if (nWhichSlice==statisticsVolumeHeader.nNumSlices-1) cout << "\n";
      }
      // load the sparse matrix
      vectorPriorSlices.push_back(zeroPriorProb);
      fstream ifs;
      ifs.open(strFilename,ios::in | ios::nocreate);
      if (ifs.good())
      {
        ifs >> vectorPriorSlices[nWhichSlice];
        ifs.close();
      }
      else
      {
        bLoadedVolumeOK=false;
      }
    }
    return(bLoadedVolumeOK);
  }

  float probability(int nTiX,int nTiY,int nTiZ)
  {
    return((float)vectorPriorSlices[nTiZ](nTiX,nTiY)/(float)statisticsVolumeHeader.nNumSubjects);
    // [[]] could speed by not normalizing return(vectorPriorSlices[nTiZ](nTiX,nTiY));
  }

  // interpolate the probability at the given talairach point for the label loaded with loadMeasureProbabilityVolume()
  // NOTE: vectorCanonicalPointsForInterpolation contains the 8 strongly connected locations in talairach space
  float probability(float fTiX, float fTiY, float fTiZ, const TypeMatrixFloat& vectorCanonicalPointsForInterpolation)
  {
    // compute the value at each of the 8 strongly connected locations
    TypeVectorFloat vectorValueAtTriplets;
    for (int i=0; i<(int)vectorCanonicalPointsForInterpolation.size(); i++)
    {
      vectorValueAtTriplets.push_back(probability(vectorCanonicalPointsForInterpolation[i][0],
                                      vectorCanonicalPointsForInterpolation[i][1],
                                      vectorCanonicalPointsForInterpolation[i][2]));
    }

    // interpolate these 8 strongly connected locations to get the value for the current floating point lcoation in talairach space
    return(statisticsVolumeHeader.interpolateLinear(fTiX,fTiY,fTiZ,vectorCanonicalPointsForInterpolation,vectorValueAtTriplets));
  }


  bool freeLabelProbabilityVolume()
  {
    vectorPriorSlices.clear();
    return(true);
  }

  // ==============================End For Baysian Classifie==========================================================================


  bool writeSlicePriors()
  {
    // loop over all priorSlices
    TypeMapSparse2DMatrix::iterator it;
    for (it=priorSlices.begin(); it!=priorSlices.end(); it++)
    {
      // form the filename
      char strFilename[500];
      int nWhichLabel=it->first;
      sprintf(strFilename,"%s/Prior-Slc%03d-Lab%03d.spmat",statisticsVolumeHeader.strStatisticsVolumeDirPath.c_str(),nCurrSliceNum,nWhichLabel);

      fstream ofs;
      ofs.open(strFilename,ios::out);
      // write the sparse matrix out to disk
      ofs << priorSlices[nWhichLabel];
      ofs.close();
    }

    return(true);
  }


  bool updateSlicePriors(MRI* pVolume)
  {
    // loop over the pixels in slice and increment the corresponding prior location
    int nTiZ=nCurrSliceNum;
    TypeMapFractionalAreaAndPoints mapFractionalAreaAndPoints;
    TypeMapFractionalAreaAndPoints::iterator it;


    //for (int nTiX=76; nTiX<77; nTiX++)
    //  for (int nTiY=65; nTiY<66; nTiY++)
    for (int nTiX=0; nTiX<statisticsVolumeHeader.nXDIM; nTiX++)
      for (int nTiY=0; nTiY<statisticsVolumeHeader.nYDIM; nTiY++)
      {
        if (statisticsVolumeHeader.AreasAndLabelsForVoxel(pVolume, nTiX,nTiY, nTiZ, mapFractionalAreaAndPoints))
        {
          for (it=mapFractionalAreaAndPoints.begin(); it!=mapFractionalAreaAndPoints.end(); it++)
          {
            int nLabel=it->first;
            float fFractionalVolume= it->second.fFractionalArea;
            priorSlices[nLabel].increment(nTiX,nTiY, fFractionalVolume);
          }
        }
      };

    return(true);
  }

  // makeTestLabelPriorSlices generates a Test label prior vol
  bool makeTestLabelPriorSlices()
  {
    unsigned char ucharFirstLabel=1;
    unsigned char ucharSecondLabel=2;
    float floatCountFirstLabel;
    float floatCountSecondLabel;

    for (int x=0; x<statisticsVolumeHeader.nXDIM; x++)
    {
      floatCountFirstLabel=1;
      floatCountSecondLabel=1;
      // make the opposing ramps for the two labels
      if (x<=64)
      {
        floatCountFirstLabel=fabs(64-x);
      }
      else
      {
        floatCountSecondLabel=fabs(64-x);
      }

      for (int y=0; y<statisticsVolumeHeader.nYDIM; y++)
      {
        priorSlices[ucharFirstLabel].setCount(x,y,floatCountFirstLabel);
        priorSlices[ucharSecondLabel].setCount(x,y,floatCountSecondLabel);
      };
    }
    return(true);
  }

  unsigned char  labelIndexMostAboveThreshold(ushort nCol, ushort nRow, float fThreshold, float& fProbability)
  {
    // iterate over all labels in labelLUT, find the label with the max probability
    TypeMapLUT::iterator it;
    TypeMapLUT::iterator itEnd=labelLUT.getMapLUT().end();
    it=labelLUT.getMapLUT().begin();
    it++; // don't include the unknown label in the comparisons
    unsigned char ucharLabelWithMaxPrior=it->first;
    float floatMaxCount=0;

    float floatCurrLabelCount;
    int nWhichLabel;


    for ( ; it!=itEnd; it++)
    {
      nWhichLabel=it->first;
      floatCurrLabelCount=priorSlices[nWhichLabel](nCol, nRow);
      if (floatCurrLabelCount>floatMaxCount)
      {
        ucharLabelWithMaxPrior=nWhichLabel;
        floatMaxCount=floatCurrLabelCount;
      }
    }

    fProbability= (float)floatMaxCount/(float)statisticsVolumeHeader.nNumSubjects;
    if (fProbability>fThreshold)
      return(ucharLabelWithMaxPrior);
    else
    {
      fProbability=1.0-fProbability;
      return(0);  //  index 0 is transparent index
    }
  }



  bool freeSlicePriors()
  {
    priorSlices.clear();
    return(true);
  }




private:
  CLabelLUT& labelLUT;
  CStatisticsVolumeHeader& statisticsVolumeHeader;
  TypeMapSparse2DMatrix priorSlices;  // holds one slice of each labels prior label probability
  int nCurrSliceNum;

  // --- used for labeling after priors are computed, holds one label's prior label probability volume
  TypeVectorSparse2DMatrix vectorPriorSlices;


};



#endif









