/**
 * @file  LayerROI.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2017/02/09 17:20:12 $
 *    $Revision: 1.21.2.5 $
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
 *
 */

#ifndef LayerROI_h
#define LayerROI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"

class FSLabel;
class vtkImageReslice;
class vtkImageMapToColors;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkImageActor;
class vtkImageData;
class vtkProp;
class LayerMRI;
class LayerPropertyROI;
class LayerSurface;

class LayerROI : public LayerVolumeBase
{
  Q_OBJECT
public:
  LayerROI( LayerMRI* layerMRI, QObject* parent = NULL  );
  virtual ~LayerROI();

  bool LoadROIFromFile( const QString& filename );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL );

  bool HasProp( vtkProp* prop );

  void SetVisible( bool bVisible = true );
  bool IsVisible();

  virtual bool HasUndo();
  virtual bool HasRedo();

  virtual void Undo();
  virtual void Redo();

  virtual void SaveForUndo(int nPlane = 0);

  inline LayerPropertyROI* GetProperty()
  {
    return (LayerPropertyROI*)mProperty;
  }

  bool SaveROI();

  void UpdateLabelData();

  virtual void SetModified();

  bool GetCentroidPosition(double* pos);

  void GetStats(int nPlane, int *count_out, float *area_out,
                LayerMRI *underlying_mri, double *mean_out, double *sd_out);

  LayerSurface* GetMappedSurface()
  {
      return m_layerMappedSurface;
  }

  void MapLabelColorData( unsigned char* colordata, int nVertexCount);

public slots:
  void UpdateOpacity();
  void UpdateColorMap();
  void UpdateThreshold();
  void SetMappedSurface(LayerSurface* s);
  void OnUpdateLabelRequested();
  void EditVertex(int nvo, bool bAdd);
  void EditVertex(const QList<int> list_nvo, bool bAdd);
  void Dilate(int nTimes = 1);
  void Erode(int nTimes = 1);
  void Open(int nTimes = 1);
  void Close(int nTimes = 1);

protected slots:
  void OnBaseVoxelEdited(const QList<int> voxel_list, bool bAdd);

protected:
  bool DoRotate( std::vector<RotationElement>& rotations );
  void DoRestore();
  void InitializeActors();
  void UpdateProperties();
  void OnLabelDataUpdated();

  virtual void OnSlicePositionChanged( int nPlane );

  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>   mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];

  LayerMRI*  m_layerSource;
  FSLabel*   m_label;
  LayerSurface* m_layerMappedSurface;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];
  int*      m_nVertexCache;
};

#endif


