/**
 * @brief i have no idea
 *
 */
/*
 * Original Author: Kevin Teich
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

#ifndef StatisticsVolumeHeader_h
#define StatisticsVolumeHeader_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "mri_ca_configFile.h"
#include "mri_ca_util.h"
#include "mri_ca_trainingParameters.h"

typedef struct NativeIndexPointTag
{
  int nNiX;
  int nNiY;
  int nNiZ;
}
TypeNativeIndexPoint;

typedef vector<TypeNativeIndexPoint> TypeVectorNativeIndexPoints;

typedef struct FractionalAreaAndPointsTag
{
  float fFractionalArea;
  TypeVectorNativeIndexPoints vectorNativeIndexPoints;
}
TypeFractionalAreaAndPoints;


// Map of label -to-> fractional area
typedef map<unsigned char, TypeFractionalAreaAndPoints>  TypeMapFractionalAreaAndPoints;


class CStatisticsVolumeHeader
{

public:
  enum NMRAutoFixerLabelType { 
    AFLTUnknown=0,
    AFLTNonWhite,
    AFLTWhite,
    AFLTEditedWhite,
    AFLTEditedNonWhite };

  bool writeStatisticsVolumeInfoFile(string strStatisticsVolumePath,
                                     bool bEchoToStdOut=false)
  {
    CConfigFile statisticsVolumeInfoConfigFile;
    // erase what is there since we may 
    // be updating the contents that were there

    if (bEchoToStdOut)
    {
      statisticsVolumeInfoConfigFile.init(strStatisticsVolumePath);
    }
    else
    {
      statisticsVolumeInfoConfigFile.init(strStatisticsVolumePath, ios::trunc);
    }
    statisticsVolumeInfoConfigFile.writeSection("StatisticsVolumeInfo",
                                                bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(nNumSubjects,
                                         "numberOfSubjects",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(nNumSlices,
                                         "numberOfSlices",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(nNumMeasures,
                                         "numberOfMeasures",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(nXDIM,"xdim",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(nYDIM,"ydim",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(strLabelLUTPath,
                                         "LabelLUT",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fSliceThickness,
                                         "fSliceThickness",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fSquarePixelSize,
                                         "fSquarePixelSize",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fStartX, "fStartX",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fEndX, "fEndX",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fStartY, "fStartY",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fEndY, "fEndY",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fStartZ, "fStartZ",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(fEndZ, "fEndZ",bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(strLabeledVolumeSourceType,
                                         "strLabeledVolumeSourceType",
                                         bEchoToStdOut);
    statisticsVolumeInfoConfigFile.writeSection("TrainedSubjectsAndMeasures",
                                                bEchoToStdOut);
    statisticsVolumeInfoConfigFile.write(arr2DstrTrainedScans,
                                         "",
                                         bEchoToStdOut);

    return(true);
  }

  bool getStatisticsVolumeInfoFile
    (string strStatisticsVolumePath,
     string strDefaultLabelLUTPath="../LabelsAutoFixer.cfg")
  {
    bool bReturnValue=true;
    CConfigFile statisticsVolumeInfoConfigFile;
    try
    {
      // Define  defaults
      nNumSubjects=0;
      arr2DstrTrainedScans.clear();
      nNumMeasures=1;
      nXDIM=128;
      nYDIM=128;
      nNumSlices=128;
      strLabelLUTPath=strDefaultLabelLUTPath;
      fSliceThickness=0.002;
      fSquarePixelSize=0.002;
      fStartX=-0.128;
      fEndX=0.128;
      fStartY=-0.128;
      fEndY=0.128;
      fStartZ=-0.128;
      fEndZ=0.128;
      strLabeledVolumeSourceType="NMRAutoFixer";
      enumLabeledVolumeSourceType=enumNMRAutoFixer;

      statisticsVolumeInfoConfigFile.init(strStatisticsVolumePath,
                                          ios::nocreate);
      statisticsVolumeInfoConfigFile.get(nNumSubjects,
                                         "StatisticsVolumeInfo",
                                         "numberOfSubjects");
      statisticsVolumeInfoConfigFile.get(arr2DstrTrainedScans,
                                         "TrainedSubjectsAndMeasures", "");
      statisticsVolumeInfoConfigFile.get(nNumSlices,
                                         "StatisticsVolumeInfo",
                                         "numberOfSlices");
      statisticsVolumeInfoConfigFile.get(nNumMeasures,
                                         "StatisticsVolumeInfo",
                                         "numberOfMeasures");
      statisticsVolumeInfoConfigFile.get(nXDIM,"StatisticsVolumeInfo","xdim");
      statisticsVolumeInfoConfigFile.get(nYDIM,"StatisticsVolumeInfo","ydim");
      statisticsVolumeInfoConfigFile.get(strLabelLUTPath,
                                         "StatisticsVolumeInfo",
                                         "LabelLUT");
      statisticsVolumeInfoConfigFile.get(fSliceThickness,
                                         "StatisticsVolumeInfo",
                                         "fSliceThickness");
      statisticsVolumeInfoConfigFile.get(fSquarePixelSize,
                                         "StatisticsVolumeInfo",
                                         "fSquarePixelSize");
      statisticsVolumeInfoConfigFile.get(fStartX,
                                         "StatisticsVolumeInfo",
                                         "fStartX");
      statisticsVolumeInfoConfigFile.get(fEndX,"StatisticsVolumeInfo","fEndX");
      statisticsVolumeInfoConfigFile.get(fStartY,
                                         "StatisticsVolumeInfo",
                                         "fStartY");
      statisticsVolumeInfoConfigFile.get(fEndY,"StatisticsVolumeInfo","fEndY");
      statisticsVolumeInfoConfigFile.get(fStartZ,
                                         "StatisticsVolumeInfo",
                                         "fStartZ");
      statisticsVolumeInfoConfigFile.get(fEndZ,"StatisticsVolumeInfo","fEndZ");
      statisticsVolumeInfoConfigFile.get(strLabeledVolumeSourceType,
                                         "StatisticsVolumeInfo",
                                         "strLabeledVolumeSourceType");
      enumLabeledVolumeSourceType=stringToLabeledVolumeSourceType
        (strLabeledVolumeSourceType);
    }
    catch (char* s)
    {
      // Define and use these as defaults
      nNumSubjects=0;
      arr2DstrTrainedScans.clear();
      nNumMeasures=1;
      nXDIM=128;
      nYDIM=128;
      nNumSlices=128;
      strLabelLUTPath=strDefaultLabelLUTPath;
      fSliceThickness=0.002;
      fSquarePixelSize=0.002;
      fStartX=-0.128;
      fEndX=0.128;
      fStartY=-0.128;
      fEndY=0.128;
      fStartZ=-0.128;
      fEndZ=0.128;
      strLabeledVolumeSourceType="NMRAutoFixer";
      enumLabeledVolumeSourceType=enumNMRAutoFixer;
    }

    // extract the directory path from the strStatisticsVolumePath and use 
    // that for the directory for the statistics volume files
    char cstrTemp[500];
    sprintf(cstrTemp,"%s",strStatisticsVolumePath.c_str());
    char* cstrStatisticsVolumeDirPath = strrchr(cstrTemp, '/');
    if (cstrStatisticsVolumeDirPath == NULL)
    {
      strStatisticsVolumeDirPath=".";
    }
    else
    {
      *cstrStatisticsVolumeDirPath='\00';
      strStatisticsVolumeDirPath=cstrTemp;
    }

    // create an MRI structure so that we can map from 
    // the Talairach space Index coordinates to Talairach space RAS coordinates
    pTalairachMRI=MRIallocHeader(1, 1, 1, MRI_UCHAR);
    // contents are meaningless, except for the values 
    // needed by MRIvoxelToWorld
    if (pTalairachMRI!=NULL)
    {
      // Ripped this stuff out of mriio,c in corRead. 
      // Its ugly but dont know how else to get MRIworldToVoxel to work.
      // Dont want to have dummy
      // COR-.info file lying around by itself confusing things on disk.
      // pTalairachMRI->slice_direction=MRI_CORONAL;

      pTalairachMRI->xstart=fStartX * 1000;
      pTalairachMRI->xend=fEndX * 1000;
      pTalairachMRI->ystart=fStartY * 1000;
      pTalairachMRI->yend=fEndY * 1000;
      pTalairachMRI->zstart=fStartZ * 1000;
      pTalairachMRI->zend=fEndZ * 1000;

      pTalairachMRI->xsize=fSquarePixelSize * 1000;
      pTalairachMRI->ysize=fSquarePixelSize * 1000;
      pTalairachMRI->zsize=fSliceThickness * 1000;
    }
    else
    {
      bReturnValue=false;
    }

    return(bReturnValue);
  }



  int mapLabelToAFLT(int nLabel)
  { // AFLTNonWhite=1, AFLTWhite, AFLTEditedWhite, AFLTEditedNonWhite
    // Only this routine is specific to the NMR filled COR files:
    if (nLabel>=WM_MIN_VAL)
    {
      if (nLabel<WM_EDITED_ON_VAL)
      {
        return(AFLTWhite);
      }
      else
      {
        if (nLabel==WM_EDITED_ON_VAL)
          return(AFLTEditedWhite);
        else
          return(AFLTUnknown);
      }
    }
    else
    {
      if (nLabel!=WM_EDITED_OFF_VAL)
      {
        if (nLabel>=0)
        {
          return(AFLTNonWhite);
        }
        else
        {
          return(AFLTUnknown);
        }
      }
      else
      {
        return(AFLTEditedNonWhite);
      }
    }
  }


  // note this routine does NOT do bounds checking - this is done by the caller
  bool nativeIndexToTalairachIndex(MRI* pNativeVolume,
                                   int nNiX, int nNiY, int nNiZ,
                                   float& fTiX, float& fTiY, float& fTiZ)
  {
    double realTrX;
    double realTrY;
    double realTrZ;
    MRIvoxelToTalairach(pNativeVolume,
                        nNiX, nNiY, nNiZ,
                        &realTrX, &realTrY, &realTrZ); // Index->TalRAS
    double realTiX;
    double realTiY;
    double realTiZ;
    MRIworldToVoxel(pTalairachMRI,
                    realTrX, realTrY, realTrZ,
                    &realTiX, &realTiY, &realTiZ);  // RAS->Index
    fTiX=realTiX;
    fTiY=realTiY;
    fTiZ=realTiZ;

    return(true);
  }


  // note this routine does NOT do bounds checking - this is done by the caller
  bool talairachIndexToNativeIndex(MRI* pNativeVolume,
                                   float fTiX, float fTiY, float fTiZ,
                                   float& fNiX, float& fNiY, float& fNiZ)
  {
    /* Strategy for transforming talairach index space to native index space
       -------------
       1 Let  talairachMRI = an MRI* pointer to a structure with voxel
       dimensions to be 2x2x2 mm
       2 Transform TalIndex->TalRAS via  MRIVoxelToWorld(talairachMRI,...)
       3 Let NativeVolume = an MRI* pointer to an MRI structure read in from
       the WM COR volume by MRIread() which should also load the talairach.xfm
       from the transforms directory
       4. Transform TalRas->NativeIndex via
       MRItalairachToVoxel(NativeVolume,...)

       Naming conventions:  Type Abbreviations
       -------------------------
       real = Real
       n = int
       f = float

       Naming conventions:     Space and Coord sys abbreviations
       -------------------------------------------
       Tr = Talairach space, RAS coordinates
       Ti = Talairach space, Index coordinates
       Nr = Native (labeled vol) space, RAS coordinates
       Ni = Native (labeled vol) space, Index coordinates
    */


    // Make sure we use the same size as Real
    double realTrX;
    double realTrY;
    double realTrZ;
    MRIvoxelToWorld(pTalairachMRI,
                    fTiX,fTiY, fTiZ,
                    &realTrX, &realTrY, &realTrZ);   // Tal Index -> Tal RAS
    double realNiX;
    double realNiY;
    double realNiZ;
    MRItalairachToVoxel(pNativeVolume,
                        realTrX,realTrY,realTrZ,
                        &realNiX, &realNiY, &realNiZ);
                        // Tal RAS -> Native Index
    fNiX=realNiX;
    fNiY=realNiY;
    fNiZ=realNiZ;
    return(true);
  }



  void mapPrevNext(float fCoordinate,
                   int nMin, int nMax,
                   float& fPrev, float& fNext)
  {
    fPrev=max<int>(nMin,floor(fCoordinate));
    fNext=min<int>(nMax-1,ceil(fCoordinate));
    if (fPrev==fNext)
    {
      if (fNext==nMax-1)
      {
        fPrev=fNext-1;
      }
      else
      {
        fNext=fPrev+1;
      }
    }
  }

  float round(float fValue)
  {
    return( (floor(fValue - 0.5) == floor(fValue) ) ? 
            ceil(fValue) : floor(fValue));
  }

  bool inBounds(int nXMax, int nYMax, int nZMax, float fX, float fY, float fZ)
  {
    if  (fX<0)  return(false);
    if  (fY<0)  return(false);
    if  (fZ<0)  return(false);
    if  (fX>=nXMax)  return(false);
    if  (fY>=nYMax)  return(false);
    if  (fZ>=nZMax)  return(false);
    return(true);
  }

  // find8ClosestVoxels: The returned canonical points must have 
  // the order of the points preserved since InterpolateLinear depends
  // on the order to compute the
  // interpolated  value
  bool find8ClosestVoxels
    (float fTiX, float fTiY, float fTiZ,
     TypeMatrixFloat& vectorCanonicalPointsForInterpolation)
  {
    vectorCanonicalPointsForInterpolation.clear();
    float fTiXNext, fTiXPrev, fTiYNext, fTiYPrev,fTiZNext, fTiZPrev;
    mapPrevNext(fTiX,0, nXDIM, fTiXPrev,  fTiXNext );
    mapPrevNext(fTiY,0, nYDIM, fTiYPrev, fTiYNext);
    mapPrevNext(fTiZ,0, nNumSlices, fTiZPrev,  fTiZNext );

    TypeVectorFloat vectPoint;
    vectPoint.clear();
    vectPoint.push_back(fTiXPrev);
    vectPoint.push_back(fTiYPrev);
    vectPoint.push_back(fTiZPrev);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXPrev);
    vectPoint.push_back(fTiYPrev);
    vectPoint.push_back(fTiZNext);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXPrev);
    vectPoint.push_back(fTiYNext);
    vectPoint.push_back(fTiZPrev);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXPrev);
    vectPoint.push_back(fTiYNext);
    vectPoint.push_back(fTiZNext);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXNext);
    vectPoint.push_back(fTiYPrev);
    vectPoint.push_back(fTiZPrev);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXNext);
    vectPoint.push_back(fTiYPrev);
    vectPoint.push_back(fTiZNext);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXNext);
    vectPoint.push_back(fTiYNext);
    vectPoint.push_back(fTiZPrev);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    vectPoint.clear();
    vectPoint.push_back(fTiXNext);
    vectPoint.push_back(fTiYNext);
    vectPoint.push_back(fTiZNext);
    vectorCanonicalPointsForInterpolation.push_back(vectPoint);
    return(true);
  }

  float interpolateLinear
    (float fTiX, float fTiY, float fTiZ,
     const TypeMatrixFloat& vectorCanonicalPointsForInterpolation,
     const TypeVectorFloat& vectorValueAtTriplets, bool bDebug=false)
  {
    float fInterpolatedValue=0.0;

    // figure out the dimensions of the voxel whose vertices are these points
    float fTiXNext, fTiXPrev, fTiYNext, fTiYPrev,fTiZNext, fTiZPrev;
    mapPrevNext(fTiX,0, nXDIM, fTiXPrev,  fTiXNext );
    mapPrevNext(fTiY,0, nYDIM, fTiYPrev, fTiYNext);
    mapPrevNext(fTiZ,0, nNumSlices, fTiZPrev,  fTiZNext );

    // interpolate by computing the weighted sum of the 8 measurements
    float fXFractionOfPrev=1.0-fabs(fTiX-fTiXPrev);
    float fXFractionOfNext=1.0-fabs(fTiXNext-fTiX);
    float fYFractionOfPrev=1.0-fabs(fTiY-fTiYPrev);
    float fYFractionOfNext=1.0-fabs(fTiYNext-fTiY);
    float fZFractionOfPrev=1.0-fabs(fTiZ-fTiZPrev);
    float fZFractionOfNext=1.0-fabs(fTiZNext-fTiZ);

    fInterpolatedValue = fXFractionOfPrev * fYFractionOfPrev * fZFractionOfPrev * vectorValueAtTriplets[0] +
                         fXFractionOfPrev * fYFractionOfPrev * fZFractionOfNext * vectorValueAtTriplets[1] +
                         fXFractionOfPrev * fYFractionOfNext * fZFractionOfPrev * vectorValueAtTriplets[2] +
                         fXFractionOfPrev * fYFractionOfNext * fZFractionOfNext * vectorValueAtTriplets[3] +
                         fXFractionOfNext * fYFractionOfPrev * fZFractionOfPrev * vectorValueAtTriplets[4] +
                         fXFractionOfNext * fYFractionOfPrev * fZFractionOfNext * vectorValueAtTriplets[5] +
                         fXFractionOfNext * fYFractionOfNext * fZFractionOfPrev * vectorValueAtTriplets[6] +
                         fXFractionOfNext * fYFractionOfNext * fZFractionOfNext * vectorValueAtTriplets[7];

    if (bDebug)
    {
      char cstrMsg[500];
      sprintf(cstrMsg,"fTiXNext, fTiXPrev, fTiYNext, fTiYPrev,fTiZNext, fTiZPrev are\n %f,%f,%f,%f,%f,%f\n",fTiXNext, fTiXPrev, fTiYNext, fTiYPrev,fTiZNext, fTiZPrev);
      cout << cstrMsg;
      sprintf(cstrMsg,"fXFractionOfPrev,fXFractionOfNext,fYFractionOfPrev, fYFractionOfNext,fZFractionOfPrev,fZFractionOfNext are:\n %f,%f,%f,%f,%f,%f",fXFractionOfPrev,fXFractionOfNext,fYFractionOfPrev, fYFractionOfNext,fZFractionOfPrev,fZFractionOfNext);

      cout << cstrMsg;
      cout << "fInterpolatedValue=" << fInterpolatedValue;
    }





    return fInterpolatedValue;
  }

  // Returns a map from label (uchar) to fractional volume that the given canonical space voxel is labeled with the label
  //  (labels which do not label the voxel are not included in the map)
  // Along with each label is a list of the Native space index coordinates of which had that label. This is used when computing the measurement
  // vector which contributes to the conditional density
  bool AreasAndLabelsForVoxel(MRI* pLabeledVolume, float fTiX, float fTiY, float fTiZ, TypeMapFractionalAreaAndPoints& mapFractionalAreaAndPoints)
  {
    bool bReturnValue=true;
    mapFractionalAreaAndPoints.clear();
    /* Strategy for computing the fractionalIntersectionVolume
       -------------------------------------------------------

      Project the 8 verticies of the Talairach voxel into native space
        These verticies are found by adding and subtracting .4999 from the Talairach index point
        These are projected by mapping them to Native space index coordinates using the talairach transform and scale factors in the mri struct
      Next the indexes of the smallest native index space rectanguloid which bounds the talairach voxel is found
        These indexes are found by rounding the projected coordinates to the indexes of the nearest native space voxel
      Next boundary conditions are handled by moving points less than the minimum to the minimum native space index coordinate
                                       and by moving points more than the maximum to the maximum natice space index coordinate
      Next the largest and smallest native space index coordinate is found in ea dim. x,y,z
      Next these native space index coord bound are used to bound an integral over the labeled volume which is used to approximate the
        amount of volume with each label (from calc we have that as the label resolution goes to infinity the error goes to zero)
        Thus we 1) Count the number of Native space voxel-labels with each label
                2) Compute fractional volume by normalizing by the number of voxels in the rectanguloid

      Notes:
        For talairach voxel centroids which map to a location outside the bounds of the Native space simply return with no
        labels and no fractionalVolumes
     */

    float arr1DXPoints[8];
    float arr1DYPoints[8];
    float arr1DZPoints[8];
    const float AREA_FACTOR = 0.4999;

    float fNiX, fNiY, fNiZ;
    talairachIndexToNativeIndex(pLabeledVolume, fTiX, fTiY, fTiZ, fNiX, fNiY, fNiZ);
    if (!inBounds(pLabeledVolume->width, pLabeledVolume->height, pLabeledVolume->depth, fNiX, fNiY, fNiZ))
    {
      bReturnValue=false;
    }
    else  // corresponding Ni point is in bounds
    {

      talairachIndexToNativeIndex(pLabeledVolume, fTiX- AREA_FACTOR, fTiY- AREA_FACTOR, fTiZ- AREA_FACTOR, arr1DXPoints[0],arr1DYPoints[0],arr1DZPoints[0]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX- AREA_FACTOR, fTiY- AREA_FACTOR, fTiZ+ AREA_FACTOR, arr1DXPoints[1],arr1DYPoints[1],arr1DZPoints[1]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX- AREA_FACTOR, fTiY+ AREA_FACTOR, fTiZ- AREA_FACTOR, arr1DXPoints[2],arr1DYPoints[2],arr1DZPoints[2]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX- AREA_FACTOR, fTiY+ AREA_FACTOR, fTiZ+ AREA_FACTOR, arr1DXPoints[3],arr1DYPoints[3],arr1DZPoints[3]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX+ AREA_FACTOR, fTiY- AREA_FACTOR, fTiZ- AREA_FACTOR, arr1DXPoints[4],arr1DYPoints[4],arr1DZPoints[4]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX+ AREA_FACTOR, fTiY- AREA_FACTOR, fTiZ+ AREA_FACTOR, arr1DXPoints[5],arr1DYPoints[5],arr1DZPoints[5]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX+ AREA_FACTOR, fTiY+ AREA_FACTOR, fTiZ- AREA_FACTOR, arr1DXPoints[6],arr1DYPoints[6],arr1DZPoints[6]);
      talairachIndexToNativeIndex(pLabeledVolume, fTiX+ AREA_FACTOR, fTiY+ AREA_FACTOR, fTiZ+ AREA_FACTOR, arr1DXPoints[7],arr1DYPoints[7],arr1DZPoints[7]);

      // debug stuff
      int nDebug=0;
      char cstrMsg[500];
      if (nDebug==1)
      {

        sprintf(cstrMsg,"Tal (%f,%f,%f) has the 8 bounding native pts:",fTiX,fTiY, fTiZ);
        cout << cstrMsg;
        for (int ii=0; ii<8; ii++)
        {
          sprintf(cstrMsg," (%f,%f,%f) ",arr1DXPoints[ii],arr1DYPoints[ii],arr1DZPoints[ii]);
          cout << cstrMsg << "\n";
        }
      }

      const float MIN_NATIVE_INDEX_X=0;
      const float MAX_NATIVE_INDEX_X=pLabeledVolume->width-1;
      const float MIN_NATIVE_INDEX_Y=0;
      const float MAX_NATIVE_INDEX_Y=pLabeledVolume->height-1;
      const float MIN_NATIVE_INDEX_Z=0;
      const float MAX_NATIVE_INDEX_Z=pLabeledVolume->depth-1;

      float fMinX;
      float fMaxX;
      float fMinY;
      float fMaxY;
      float fMinZ;
      float fMaxZ;

      // the indexes of the smallest native index space rectanguloid which bounds the talairach voxel is found
      for (int i=0; i<8; i++)
      {
        arr1DXPoints[i]=round(arr1DXPoints[i]);
        arr1DYPoints[i]=round(arr1DYPoints[i]);
        arr1DZPoints[i]=round(arr1DZPoints[i]);

        if (arr1DXPoints[i]<MIN_NATIVE_INDEX_X) arr1DXPoints[i]=MIN_NATIVE_INDEX_X;
        if (arr1DXPoints[i]>MAX_NATIVE_INDEX_X) arr1DXPoints[i]=MAX_NATIVE_INDEX_X;
        if (arr1DYPoints[i]<MIN_NATIVE_INDEX_Y) arr1DYPoints[i]=MIN_NATIVE_INDEX_Y;
        if (arr1DYPoints[i]>MAX_NATIVE_INDEX_Y) arr1DYPoints[i]=MAX_NATIVE_INDEX_Y;
        if (arr1DZPoints[i]<MIN_NATIVE_INDEX_Z) arr1DZPoints[i]=MIN_NATIVE_INDEX_Z;
        if (arr1DZPoints[i]>MAX_NATIVE_INDEX_Z) arr1DZPoints[i]=MAX_NATIVE_INDEX_Z;

        if (i==0)
        {
          fMinX=arr1DXPoints[i];
          fMaxX=arr1DXPoints[i];
          fMinY=arr1DYPoints[i];
          fMaxY=arr1DYPoints[i];
          fMinZ=arr1DZPoints[i];
          fMaxZ=arr1DZPoints[i];
        }
        else
        {
          if (arr1DXPoints[i]>fMaxX) fMaxX=arr1DXPoints[i];
          if (arr1DXPoints[i]<fMinX) fMinX=arr1DXPoints[i];
          if (arr1DYPoints[i]>fMaxY) fMaxY=arr1DYPoints[i];
          if (arr1DYPoints[i]<fMinY) fMinY=arr1DYPoints[i];
          if (arr1DZPoints[i]>fMaxZ) fMaxZ=arr1DZPoints[i];
          if (arr1DZPoints[i]<fMinZ) fMinZ=arr1DZPoints[i];
        }
      }


      if (nDebug==1)
      {
        cout << "After adjusting them:\n";
        for (int ii=0; ii<8; ii++)
        {
          sprintf(cstrMsg," (%f,%f,%f) ",arr1DXPoints[ii],arr1DYPoints[ii],arr1DZPoints[ii]);
          cout << cstrMsg << "\n";
        }
        cout << "Tal voxel is labeled:\n";
      }

      int nNumNiVoxelsInTightestBoundingRectanguloid=(int)((fMaxX-fMinX+1)*(fMaxY-fMinY+1)*(fMaxZ-fMinZ+1));
      if  (nNumNiVoxelsInTightestBoundingRectanguloid<=0)   nNumNiVoxelsInTightestBoundingRectanguloid=1;

      // Count the number of voxels of each label
      //  and find the set of points which label the canonical voxel with that label
      unsigned char ucharLabel;
      TypeMapFractionalAreaAndPoints::iterator it;
      mapFractionalAreaAndPoints.clear();
      for (int nNiX=(int)fMinX; nNiX<=fMaxX; nNiX++)
        for (int nNiY=(int)fMinY; nNiY<=fMaxY; nNiY++)
          for (int nNiZ=(int)fMinZ; nNiZ<=fMaxZ; nNiZ++)
          {
            TypeNativeIndexPoint aPoint;
            aPoint.nNiX=nNiX;
            aPoint.nNiY=nNiY;
            aPoint.nNiZ=nNiZ;
            ucharLabel=MRIvox(pLabeledVolume,nNiX,nNiY, nNiZ);
            if (enumLabeledVolumeSourceType==enumNMRAutoFixer)
            {
              ucharLabel=mapLabelToAFLT(ucharLabel);
            }

            if ((it=mapFractionalAreaAndPoints.find(ucharLabel))==mapFractionalAreaAndPoints.end())
            {
              TypeFractionalAreaAndPoints fractionalAreaAndPoints;
              fractionalAreaAndPoints.vectorNativeIndexPoints.push_back(aPoint);
              fractionalAreaAndPoints.fFractionalArea=1.0;
              mapFractionalAreaAndPoints.insert(TypeMapFractionalAreaAndPoints::value_type(ucharLabel,fractionalAreaAndPoints));
            }
            else
            {
              it->second.fFractionalArea=it->second.fFractionalArea+1;  // use fractional area to count the number of voxels with the same label
              it->second.vectorNativeIndexPoints.push_back(aPoint);
            }
          }

      // for ea label compute an approx fraction of volume intersection with the canonical voxel
      for (it=mapFractionalAreaAndPoints.begin(); it!=mapFractionalAreaAndPoints.end(); it++)
      {
        it->second.fFractionalArea=(float)it->second.fFractionalArea/(float)nNumNiVoxelsInTightestBoundingRectanguloid;

        if (nDebug==1)
        {
          int nLabel=it->first;
          float fFractionalVolume= it->second.fFractionalArea;
          cout << "label " << nLabel << " w frac vol " << fFractionalVolume;
          cout << "\n\n";
        }

      }
    }
    return(true);
  }




  ~CStatisticsVolumeHeader()
  {
    if (pTalairachMRI!=NULL)
    {
      MRIfree(&pTalairachMRI);
    }
  }


  bool setLabeledVolumeSourceType(LabeledVolumeSourceType enumNewLabeledVolumeSourceType)
  {
    enumLabeledVolumeSourceType=enumNewLabeledVolumeSourceType;
    return(true);
  }


public:
  int nNumSlices;
  int nXDIM;
  int nYDIM;
  int nNumSubjects;
  int nNumMeasures;
  float fSliceThickness;
  ;
  float fSquarePixelSize;
  float fStartX;
  float fEndX;
  float fStartY;
  float fEndY;
  float fStartZ;
  float fEndZ;
  string strLabelLUTPath;
  TypeMatrixString   arr2DstrTrainedScans;
  string strStatisticsVolumeDirPath;
  MRI* pTalairachMRI;
  string strLabeledVolumeSourceType;
  LabeledVolumeSourceType enumLabeledVolumeSourceType;

};

#endif
