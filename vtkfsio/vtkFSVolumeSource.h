/**
 * @file  vtkFSVolumeSource.h
 * @brief Source for FreeSurfer MRI volumes
 *
 * This will read any FS MRI volume with the MRIRead command, and
 * provide a VTK StructuredPointsSource output port. The output is in
 * 'index' space, with the corner at 0,0,0. There are also functions
 * for getting the 'RAS' space transforms, which should be used to
 * display the output properly in RAS space.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:56 $
 *    $Revision: 1.8 $
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

#ifndef __vtkFSVolumeSource_h
#define __vtkFSVolumeSource_h

#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsSource.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
extern "C" {
#include "mri.h"
}

// These are defined in the mri headers somewhere. You will probably
// have to do this in files that include this file.
#ifdef X
#undef X
#endif
#ifdef Y
#undef Y
#endif


class vtkFSVolumeSource : public vtkStructuredPointsSource {
public:

  static vtkFSVolumeSource *New();
  vtkTypeRevisionMacro( vtkFSVolumeSource, vtkStructuredPointsSource );

  // Description:
  // This will call the MRIread function from the FS library.
  void MRIRead ( char const* ifnVolume );

  // Description:
  // This writes the volume. The volume type is determined by the
  // filename extension or format. If no filename is given, it will
  // use the orignal filename, overwriting the original volume.
  void MRIWrite ( char const* ifnVolume );
  void MRIWrite ();

  // Description:
  // If set to true, the output image will have the proper spacing
  // from the MRI {x,y,z}size values. If false, spacing will just be
  // 1,1,1.
  void ActualSpacingOn ();
  void ActualSpacingOff ();

  // Description:
  // Coordinate conversion. RAS space is defined by various header
  // metadata and used by the Layer to display things in the right
  // space.
  int ConvertIndexToRAS ( float iIdxX, float iIdxY, float iIdxZ,
                          float& oRASX, float& oRASY, float& oRASZ );
  int ConvertRASToIndex ( float iRASX, float iRASY, float iRASZ,
                          float& oIdxX, float& oIdxY, float& oIdxZ );
  int ConvertRASToIndex ( float iRASX, float iRASY, float iRASZ,
                          int& oIdxX, int& oIdxY, int& oIdxZ );

  void GetRASBounds ( float oRASBounds[6] );

  void GetUnscaledRASBounds ( float oUnscaledRASBounds[6] );

  double* GetVoxelToRASMatrix () {
    return mVoxelToRASMatrix;
  }
 
  // Description:
  // Returns the RAS matrix, but without pixel scaling
  void GetUnscaledVoxelToRASMatrix ( double oMatrix[16] );
  
  double* GetRASToVoxelMatrix () {
    return mRASToVoxelMatrix;
  }

  // Description:
  // Returns the RAS matrix, but without pixel scaling
  void GetUnscaledRASToVoxelMatrix ( double* oMatrix );

  void SetVoxelToRASMatrix ( vtkMatrix4x4& iMatrix );

  float GetRASCenterX ();
  float GetRASCenterY ();
  float GetRASCenterZ ();

  // Description:
  // Get and set the value at the given index.
  float GetValueAtIndex ( float iIdxX, float iIdxY, float iIdxZ );
  void  SetValueAtIndex ( float iIdxX, float iIdxY, float iIdxZ,
                          float iValue );
  float GetValueAtIndex (float iIdxX, float iIdxY, float iIdxZ, float iIdxFrame );


  // Description:
  // Get the min and max value in the volume.
  float GetMinValue ();
  float GetMaxValue ();

  // Description:
  // Get the pixel size in RAS space from the metadata.
  float GetPixelSizeX ();
  float GetPixelSizeY ();
  float GetPixelSizeZ ();
  float GetPixelSize ( const int iDimension );

  // Description:
  // Returns the number of elements in a dimension. This is in index
  // space.
  float GetXDimension ();
  float GetYDimension ();
  float GetZDimension ();
  int*  GetDimensions ();

  // Description:
  // Returns a best guess value increment for a GUI.
  float GetPreferredValueIncrement ();

protected:

  vtkFSVolumeSource ();
  ~vtkFSVolumeSource () {};

  // Overriding to output the volume data.
  void Execute ();
  void ExecuteInformation ();

  // If we have a new MRI, we copy the data to our local storage.
  void CopyMRIToImage ();
  void CopyMatricesFromMRI ();

  // Pointer to the MRI.
  MRI* mMRI;

  // Whether or not to use actual pixel spacing in the output image.
  bool mbUseActualPixelSpacing;

  // The local copy of the SP data.
  vtkImageData* mImageData;

  //   [  0  1  2  3 ]
  //   [  4  5  6  7 ]
  //   [  8  9 10 11 ]
  //   [ 12 13 14 15 ]
  double mVoxelToRASMatrix[16];
  double mRASToVoxelMatrix[16];

  // RAS bounds.
  bool  mbBoundsCacheDirty;
  float mRASBounds[6];

private:
  vtkFSVolumeSource(const vtkFSVolumeSource&);  // Not implemented.
  void operator=(const vtkFSVolumeSource&);  // Not implemented.
};

#endif


