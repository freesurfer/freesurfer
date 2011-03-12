/**
 * @file  LivewireTool.h
 * @brief LivewireTool data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:52 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef LivewireTool_h
#define LivewireTool_h

#include <vtkSmartPointer.h>

#include <string>
#include <vector>

class vtkImageClip;
class vtkImageData;
class vtkDijkstraImageGeodesicPath;
class vtkPoints;
class vtkImageChangeInformation;

class LivewireTool
{
public:
  LivewireTool( );
  virtual ~LivewireTool();

  void UpdateImageDataInfo( vtkImageData* image, int nPlane, int nSlice );

  void GetLivewirePoints( double* pt1_in, double* pt2_in, vtkPoints* pts_out );

  void GetLivewirePoints( vtkImageData* image, int nPlane, int nSlice,
                          double* pt1_in, double* pt2_in, vtkPoints* pts_out );

protected:
  int  m_nPlane;
  int  m_nSlice;
  vtkSmartPointer<vtkDijkstraImageGeodesicPath>  m_path;
  vtkImageData*         m_imageData;
  vtkImageData*         m_imageSlice;
};

#endif


