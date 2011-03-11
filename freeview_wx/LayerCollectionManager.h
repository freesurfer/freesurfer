/**
 * @file  LayerCollectionManager.h
 * @brief Manage the collections of layers, such as volumes, surfaces and way points, etc.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
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

#ifndef LayerCollectionManager_h
#define LayerCollectionManager_h

#include <string>
#include <vector>
#include "Listener.h"
#include "Broadcaster.h"

class LayerMRI;
class LayerCollection;
class vtkRenderer;
class Layer;

class LayerCollectionManager : public Listener, public Broadcaster
{
public:
  LayerCollectionManager();
  virtual ~LayerCollectionManager();

  void Append2DProps( vtkRenderer* renderer, int nImagePlane );
  void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  virtual void DoListenToMessage( std::string const iMsg, void* iData, void* sender );

  LayerCollection* GetLayerCollection( std::string strType );

  bool SetSlicePosition( int nPlane, double dPos, bool bRoundToGrid = true );
  bool SetSlicePosition( double* pos );
  bool OffsetSlicePosition( int nPlane, double dPosDiff, bool bRoundToGrid = true );

  void RefreshSlices();

  bool HasAnyLayer();
  bool HasLayer( std::string type );

  std::vector<Layer*> GetAllLayers();

protected:
  std::vector<LayerCollection*> m_layerCollections;
  LayerCollection*    m_layerCollectionMRI;
};

#endif


