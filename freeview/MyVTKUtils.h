/**
 * @file  MyVTKUtils.h
 * @brief Misc utility class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:52 $
 *    $Revision: 1.1 $
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

#ifndef MyVTKUtils_h
#define MyVTKUtils_h

#include <math.h>
#include <vector>
#include <QStringList>

class vtkRenderer;
class vtkRenderWindow;
class vtkImageData;
class vtkActor;
class vtkVolume;
class vtkPoints;

class MyVTKUtils
{
public:
  MyVTKUtils()
  {}
  static bool VTKScreenCapture( vtkRenderWindow* renderWindow, vtkRenderer* renderer,
                                const char* filename, bool bAntiAliasing = false, int nMag = 1 );
  static void ViewportToWorld( vtkRenderer* renderer, double x, double y,
                               double& world_x, double& world_y, double& world_z );
  
  static void ViewportToWorld( vtkRenderer* renderer, double x, double y, double z, 
                               double& world_x, double& world_y, double& world_z );
  
  static void NormalizedViewportToWorld( vtkRenderer* renderer, double x, double y,
                                         double& world_x, double& world_y, double& world_z );
  
  static void NormalizedViewportToWorld( vtkRenderer* renderer, double x, double y, double z,
                                         double& world_x, double& world_y, double& world_z );

  static void WorldToViewport( vtkRenderer* renderer, double world_x, double world_y,
                               double world_z, double& x, double& y, double& z );

  static bool BuildContourActor( vtkImageData* data_in, double dTh1, double dTh2, vtkActor* actor_out,
                                 int nSmoothIterations = 0, int* ext = NULL, bool bAllRegion = false );

  static bool BuildVolume( vtkImageData* data_in, double dTh1, double dTh2, vtkVolume* vol_out );

  static void GetLivewirePoints( vtkImageData* image_in, int nPlane_in, int nSlice_in,
                                 double* pt1_in, double* pt2_in, vtkPoints* pts_out );
};

#endif
