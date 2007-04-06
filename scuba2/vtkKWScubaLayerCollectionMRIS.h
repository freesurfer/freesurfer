/**
 * @file  vtkKWScubaLayerCollectionMRIS.h
 * @brief Implementation for MRIS viewers.
 *
 * In 2D, the MRIS is viewed as the intersection of the surface with a
 * slice. In 3D, the MRIS is viewed as a fully 3D mesh surface.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:05 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2002-2007,
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

#ifndef vtkKWScubaLayerCollectionMRIS_h
#define vtkKWScubaLayerCollectionMRIS_h

#include "vtkKWScubaLayerCollection.h"
#include "ScubaCollectionPropertiesMRIS.h"
#include "vtkFSSurfaceSource.h"

class vtkTransformPolyDataFilter;
class vtkAlgorithmOutput;
class vtkDecimatePro;

class vtkKWScubaLayerCollectionMRIS : public vtkKWScubaLayerCollection
                    //BTX
                    , public ScubaCollectionPropertiesMRIS
		    //ETX
{

 public:

  static vtkKWScubaLayerCollectionMRIS* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayerCollectionMRIS, 
			vtkKWScubaLayerCollection );

  // Description:
  // Set the surface file name.
  void SetSurfaceFileName ( const char* ifnVolume );

  // Layers will get the source objects with these calls. These
  // impelement ScubaCollectionPropertiesMRIS.
  vtkFSSurfaceSource* GetSource () const;
  vtkAlgorithmOutput* GetNormalModeOutputPort () const;
  vtkAlgorithmOutput* GetFastModeOutputPort () const;
  vtkPolyData* GetNormalModeOutput () const;
  vtkPolyData* GetFastModeOutput () const;

 protected:

  vtkKWScubaLayerCollectionMRIS ();
  ~vtkKWScubaLayerCollectionMRIS ();

  //BTX
  virtual vtkKWScubaLayer*
    MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
  //ETX

  void LoadSurfaceFromFileName ();

 private:

  //BTX
  vtkFSSurfaceSource* mSource;
  vtkDecimatePro* mDecimator;
  std::string mfnSurface;
  //ETX

};

#endif
