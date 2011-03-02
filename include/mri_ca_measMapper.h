/**
 * @file  mri_ca_measMapper.h
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


#ifndef MeasMapper_h
#define MeasMapper_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>


#include "mri_ca_measure.h"
#include "mri_ca_util.h"
#include "mri_ca_sparse2DGausDistMatrix.h"
#include "mri_ca_labelLUT.h"
#include "mri_ca_measureVolume.h"
#include "mri_ca_configFile.h"
#include "mri_ca_trainingParameters.h"
#include "mri_ca_statisticsVolumeHeader.h"
#include <map>
#include <vector>
extern "C"
{
#include "mri.h"
}

using namespace std;


typedef vector<CMeasureVolume> TypeVectorMeasureVolume;
typedef map<int, CSparse2DGausDistMatrix> TypeMapSparse2DGausDistMatrix;
typedef vector<CSparse2DGausDistMatrix> TypeVectorSparse2DGausDistMatrix;

class CMeas2DensMapper
{



public:


  CMeas2DensMapper(CLabelLUT& labelNewLUT, CStatisticsVolumeHeader& newStatisticsVolumeHeader) : labelLUT(labelNewLUT), statisticsVolumeHeader(newStatisticsVolumeHeader)
  {};

  bool setNumberOfMeasures(int nNumOfMeasures)
  {
    CMeasureVolume newMeasureVolume;
    newMeasureVolume.pVolume=NULL;
    newMeasureVolume.strMeasureFileDir="";
    newMeasureVolume.measureType=CMeasureVolume::Unknown;

    measureVolumes.clear();
    for (int i=0; i<nNumOfMeasures; i++)
    {
      measureVolumes.push_back(newMeasureVolume);
    }

    return(true);
  }


  // load a measure and return the MRI pointer if desired (ppMeasuredVolume!=NULL)
  bool loadMeasure(int nMeasureVolumeIndex, string strMeasureDir, CMeasureVolume::MeasureType theMeasureType, MRI** ppMeasuredVolume=NULL)
  {
    bool bMeasureLoaded=false;
    if (nMeasureVolumeIndex<(int)measureVolumes.size())
    {
      measureVolumes[nMeasureVolumeIndex].strMeasureFileDir=strMeasureDir;
      measureVolumes[nMeasureVolumeIndex].measureType=theMeasureType;

      char strMeasDir[500];
      sprintf(strMeasDir,"%s",strMeasureDir.c_str());
      measureVolumes[nMeasureVolumeIndex].pVolume=MRIread(strMeasDir);
      if (measureVolumes[nMeasureVolumeIndex].pVolume!=NULL)
      {
        if (ppMeasuredVolume!=NULL)
        {
          *ppMeasuredVolume=measureVolumes[nMeasureVolumeIndex].pVolume;
        }
        bMeasureLoaded=true;
      }
      else
      {
        bMeasureLoaded=false;

      }
    }

    return(bMeasureLoaded);
  }

  // -----------------------------For Baysian Classifier------------------------------------------------------------------------------
  bool loadMeasureDensityVolume(int nWhichLabel, bool bDebug=false)
  {
    bool bLoadedVolumeOK=true;
    vectorDensitySlices.clear();

    //  create an initially zero prior and then use >> to update  it if the file exists
    CSparse2DGausDistMatrix zeroDensityMatrix(statisticsVolumeHeader.nXDIM,statisticsVolumeHeader.nYDIM,statisticsVolumeHeader.nNumMeasures);

    for (int nWhichSlice=0; nWhichSlice<statisticsVolumeHeader.nNumSlices; nWhichSlice++)
    {
      // form the filename for this slice
      char strFilename[500];
      sprintf(strFilename,"%s/Density-Slc%03d-Lab%03d.spmat",statisticsVolumeHeader.strStatisticsVolumeDirPath.c_str(),nWhichSlice,nWhichLabel);

      if (bDebug==false)
      {
        cout << "    Loading " << strFilename << "\n"<< flush;
      }
      else
      {
        cout << "." << flush;
        if (nWhichSlice==statisticsVolumeHeader.nNumSlices-1) cout << "\n";
      }
      // load the sparse matrix
      vectorDensitySlices.push_back(zeroDensityMatrix);
      fstream ifs;
      ifs.open(strFilename,ios::in | ios::nocreate);
      if (ifs.good())
      {
        ifs >> vectorDensitySlices[nWhichSlice];
        ifs.close();
      }
      else
      {
        bLoadedVolumeOK=false;
      }
    }
    return(bLoadedVolumeOK);
  }


  void formMeasureVector(int nNiX,int nNiY,int nNiZ, TypeVectorFloat& vectorMeasurement)
  {
    vectorMeasurement.clear();
    int nMaxMeasureNum=measureVolumes.size();
    for (int nWhichMeasure=0; nWhichMeasure<nMaxMeasureNum; nWhichMeasure++) // create a zero vector of the approp size
    {
      vectorMeasurement.push_back(MRIvox(measureVolumes[nWhichMeasure].pVolume, nNiX,nNiY, nNiZ));
    }
  }

  float density(const TypeVectorFloat& measure, int nTiX,int nTiY,int nTiZ, bool bDebugMeanVariance=false)
  {
    float fDensity=vectorDensitySlices[nTiZ].density(measure,nTiX,nTiY);
    if (bDebugMeanVariance)
    {
      char cstrMsg[500];
      float fMean=vectorDensitySlices[nTiZ].extractMeanComponent(nTiX,nTiY,0);
      float fVariance=vectorDensitySlices[nTiZ].extractVarianceComponent(nTiX,nTiY,0);
      sprintf(cstrMsg,"  Canonical voxel (%d,%d,%d) has mean %g variance %g, --> intensity %g has density %g\n", nTiX, nTiY, nTiZ,fMean, fVariance, measure[0], fDensity);
      cout << cstrMsg;

    }
    return(fDensity);
  }

  // interpolate the density at the given talairach point for the measure given by the given native index point using the
  // measures loaded and for the label loaded with loadMeasureDensityVolume()
  float density(int nNiX,int nNiY,int nNiZ, float fTiX, float fTiY, float fTiZ, const TypeMatrixFloat& vectorCanonicalPointsForInterpolation, bool bDebugMeanVariance=false)
  {
    bool bDebug=false;
    // get the measurement vector
    TypeVectorFloat vectorMeasurement;
    formMeasureVector(nNiX,nNiY,nNiZ,vectorMeasurement);

    if (bDebugMeanVariance)
    {
      char cstrMsg[500];
      sprintf(cstrMsg,"Point (%d,%d,%d) has its density interpolated from these voxels:\n", nNiX, nNiY, nNiZ);
      cout << cstrMsg;
    }

    // compute the value at each of the 8 strongly connected locations
    TypeVectorFloat vectorValueAtTriplets;
    for (int i=0; i<(int)vectorCanonicalPointsForInterpolation.size(); i++)
    {
      vectorValueAtTriplets.push_back(density(vectorMeasurement,
                                              vectorCanonicalPointsForInterpolation[i][0],
                                              vectorCanonicalPointsForInterpolation[i][1],
                                              vectorCanonicalPointsForInterpolation[i][2], bDebugMeanVariance));

    }
#if 0
    if (((nNiX==153) && (nNiY==131)) ||
        ((nNiX==89) && (nNiY==116))  ||
        ((nNiX==69) && (nNiY==32))   ||
        ((nNiX==124) && (nNiY==133))      )
    {
      cout << "vectorValueAtTriplets=\n";
      printType(vectorValueAtTriplets);
      bDebug=true;
    }

#endif




    // interpolate these 8 strongly connected locations to get the value for the current floating point lcoation in talairach space
    float fDensity=statisticsVolumeHeader.interpolateLinear(fTiX,fTiY,fTiZ,vectorCanonicalPointsForInterpolation,vectorValueAtTriplets, bDebug);
    if (bDebugMeanVariance)
    {
      char cstrMsg[500];
      sprintf(cstrMsg,"Overall interpolated conditional density is %g\n", fDensity);
      cout << cstrMsg;
    }

    return(fDensity);

  }

  bool freeMeasureDensityVolume()
  {
    vectorDensitySlices.clear();
    return(true);
  }

  // ==============================End For Baysian Classifie==========================================================================

  // NOTE: uses the number of measures defined in statisticsVolumeHeader to avoid forcing the user of this library to call
  //setNumberOfMeasures first, plus the value is stored in the statistics volume header file anyway
  bool loadSliceMeasureDensities(int nSliceNumber, int nWidth, int nHeight)
  {
    mapDensitySlices.clear();
    nCurrSliceNum=nSliceNumber;
    // iterate over all slice numbers in labelLUT
    TypeMapLUT::iterator it;

    for (it=labelLUT.getMapLUT().begin(); it!=labelLUT.getMapLUT().end(); it++)
    {
      // form the filename
      char strFilename[500];
      int nWhichLabel=it->first;
      sprintf(strFilename,"%s/Density-Slc%03d-Lab%03d.spmat",statisticsVolumeHeader.strStatisticsVolumeDirPath.c_str(),nCurrSliceNum,nWhichLabel);

      // load the sparse matrix
      //  create an initially zero prior and then use >> to update  it if the file exists
      CSparse2DGausDistMatrix zeroDensityMatrix(nWidth,nHeight,statisticsVolumeHeader.nNumMeasures);

      mapDensitySlices.insert(TypeMapSparse2DGausDistMatrix::value_type(nWhichLabel, zeroDensityMatrix));
      fstream ifs;
      ifs.open(strFilename,ios::in | ios::nocreate);
      if (ifs.good())
      {
        ifs >> mapDensitySlices[nWhichLabel];
        ifs.close();
      }
    }
    return(true);
  }



  bool writeSliceMeasureDensities()
  {

    // loop over all slicePriors
    TypeMapSparse2DGausDistMatrix::iterator it;
    for (it=mapDensitySlices.begin(); it!=mapDensitySlices.end(); it++)
    {
      // form the filename
      char strFilename[500];
      int nWhichLabel=it->first;
      sprintf(strFilename,"%s/Density-Slc%03d-Lab%03d.spmat",statisticsVolumeHeader.strStatisticsVolumeDirPath.c_str(),nCurrSliceNum,nWhichLabel);

      fstream ofs;
      ofs.open(strFilename,ios::out);
      // write the sparse matrix out to disk
      ofs << mapDensitySlices[nWhichLabel];
      ofs.close();
    }

    return(true);
  }





  bool freeSliceMeasureDensities()
  {
    mapDensitySlices.clear();
    return(true);
  }




  TypeVectorFloat averageMeasureVectorAtPoints(TypeVectorNativeIndexPoints& vectorNativeIndexPoints)
  {
    TypeVectorFloat vectorAverageMeasure;
    vectorAverageMeasure.clear();
    int nMaxMeasureNum=measureVolumes.size();
    for (int nWhichMeasure=0; nWhichMeasure<nMaxMeasureNum; nWhichMeasure++) // create a zero vector of the approp size
    {
      vectorAverageMeasure.push_back(0.0);
    }

    int nVectorNativeIndexPointsSize=vectorNativeIndexPoints.size();
    // add the measures up and then compute the average measure.
    for (int i=0; i<nVectorNativeIndexPointsSize; i++) // loop over all points
    {
      for (int nWhichMeasure=0; nWhichMeasure<nMaxMeasureNum; nWhichMeasure++)
      {
        double dMeasureValue=(double)MRIvox(measureVolumes[nWhichMeasure].pVolume,
                                            vectorNativeIndexPoints[i].nNiX,
                                            vectorNativeIndexPoints[i].nNiY,
                                            vectorNativeIndexPoints[i].nNiZ);
        vectorAverageMeasure[nWhichMeasure]+=dMeasureValue;
      }
    }

    for (int nWhichMeasure=0; nWhichMeasure<nMaxMeasureNum; nWhichMeasure++)
    { // normalize by the number of points
      vectorAverageMeasure[nWhichMeasure]=(double)vectorAverageMeasure[nWhichMeasure]/(double)nVectorNativeIndexPointsSize;
    }

    return(vectorAverageMeasure);
  }

  bool updateMeasureDensities(MRI* pVolume)
  {
    TypeVectorFloat measureVector;
    // loop over locations in the labeledSlice
    TypeMapFractionalAreaAndPoints mapFractionalAreaAndPoints;
    TypeMapFractionalAreaAndPoints::iterator it;
    int nTiZ=nCurrSliceNum;

    for (int nTiX=0; nTiX<statisticsVolumeHeader.nXDIM; nTiX++)
      for (int nTiY=0; nTiY<statisticsVolumeHeader.nYDIM; nTiY++)
        //      for (int nTiX=76; nTiX<77; nTiX++)
        //        for (int nTiY=65; nTiY<66; nTiY++)
      {
        // [[]]
#if 0
        if ((nTiZ==77) && (nTiY==67) && (nTiX==50))
        {
          char cstrMsg[500];
          sprintf(cstrMsg, "neg var -- check out label 1");
          cout << cstrMsg << "\n";

        }
#endif
        if (statisticsVolumeHeader.AreasAndLabelsForVoxel(pVolume, nTiX, nTiY, nTiZ, mapFractionalAreaAndPoints))
        {
          for (it=mapFractionalAreaAndPoints.begin(); it!=mapFractionalAreaAndPoints.end(); it++)
          {
            int nLabel=it->first;
            measureVector=averageMeasureVectorAtPoints(it->second.vectorNativeIndexPoints);
            mapDensitySlices[nLabel].insertMeasure(measureVector,nTiX,nTiY,it->second.fFractionalArea);
          }
        }
      };

    return(true);

  }

  bool makeTestLabelConditionalDensitySlices()
  {
    unsigned char ucharLabel=1;
    TypeVectorFloat vectFloatMeasureMean;
    TypeMatrixFloat matrixFloatMeasureVariance;
    float floatNewNumMeasuresEntered=100.0;
    int nMaxMeasureNum=measureVolumes.size();

    for (int x=0; x<statisticsVolumeHeader.nXDIM; x++)
    {
      // create mean vector and cov matrix for insertion into the test 2D Spars matrix of gaussian distributions
      double dMean=fabs(64-x); // want the mean to be equal to the RAS coord to be equal so we can probe with tkmedit
      vectFloatMeasureMean.clear();
      matrixFloatMeasureVariance.clear();
      for (int nRow=0; nRow<nMaxMeasureNum; nRow++)
      {
        vectFloatMeasureMean.push_back(dMean*(nRow+1)); // higher vector elements are simply 2 or 3 etc times the value in the first vect component
      }
      for (int nCol=0; nCol<nMaxMeasureNum; nCol++)
      {
        matrixFloatMeasureVariance.push_back(vectFloatMeasureMean);
      }

      for (int y=0; y<statisticsVolumeHeader.nYDIM; y++)
      {
        mapDensitySlices[ucharLabel].setGausParams(x,y,vectFloatMeasureMean,matrixFloatMeasureVariance,floatNewNumMeasuresEntered);
      };
    }
    return(true);
  }

  unsigned char extractMeanComponent(ushort x, ushort y, int nWhichLabel, int nWhichComponent)
  {
    double dMeanComponent=mapDensitySlices[nWhichLabel].extractMeanComponent(x,y,nWhichComponent);
    return(( unsigned char)((int)dMeanComponent));
  }

  double extractVarianceComponent(ushort x, ushort y, int nWhichLabel, int nWhichComponent)
  {
    double dVarianceComponent=mapDensitySlices[nWhichLabel].extractVarianceComponent(x,y,nWhichComponent);
    return(dVarianceComponent);
  }



  bool freeMeasures()
  {
    // loop over all measure volumes and free the volumes
    int nMaxMeasureNum=measureVolumes.size();
    for (int nWhichMeasure=0; nWhichMeasure<nMaxMeasureNum; nWhichMeasure++)
    {
      if (measureVolumes[nWhichMeasure].pVolume!=NULL)
      {
        MRIfree(&(measureVolumes[nWhichMeasure].pVolume));
      }
    }

    return(true);
  }

private:
  CLabelLUT& labelLUT;
  CStatisticsVolumeHeader& statisticsVolumeHeader;

  TypeVectorMeasureVolume measureVolumes;



  // string strDensityDirectory;

  //  TypeMapGaussDist mapDensityDistributions;
  TypeMapSparse2DGausDistMatrix mapDensitySlices;
  int nCurrSliceNum;

  // --- used for labeling after priors are computed, holds one label's density volume
  TypeVectorSparse2DGausDistMatrix vectorDensitySlices;

};







#endif































