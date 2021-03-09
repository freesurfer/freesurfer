/**
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
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
 *
 */

#ifndef LayerROI_h
#define LayerROI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"
#include <QVector>

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

 
#include "label.h"


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

  virtual void SaveForUndo(int nPlane = -1, bool bAllFrames = false );

  inline LayerPropertyROI* GetProperty()
  {
    return (LayerPropertyROI*)mProperty;
  }

  bool SaveROI();

  void UpdateLabelData();

  bool GetCentroidPosition(double* pos);

  void GetStats(int nPlane, int *count_out, float *area_out,
                LayerMRI *underlying_mri, double *mean_out, double *sd_out);

  LayerSurface* GetMappedSurface()
  {
    return m_layerMappedSurface;
  }

  void MapLabelColorData( unsigned char* colordata, int nVertexCount);

  LABEL*  GetRawLabel();

  LayerMRI* GetRefMRI()
  {
    return m_layerSource;
  }

public slots:

  virtual void SetModified();
  void UpdateOpacity();
  void UpdateColorMap();
  void UpdateThreshold();
  void SetMappedSurface(LayerSurface* s);
  void OnUpdateLabelRequested();
  void EditVertex(int nvo, bool bAdd);
  void EditVertex(const QVector<int> list_nvo, bool bAdd);
  void Dilate(int nTimes = 1);
  void Erode(int nTimes = 1);
  void Open(int nTimes = 1);
  void Close(int nTimes = 1);
  void Resample();
  void Clear();
  void OnSurfaceDestroyed(QObject* obj);

protected slots:
  void OnBaseVoxelEdited(const QVector<int>& voxel_list, bool bAdd);

protected:
  bool DoRotate( std::vector<RotationElement>& rotations );
  void DoRestore();
  void InitializeActors();
  void UpdateProperties();
  void OnLabelDataUpdated();
  void UpdateFilteredImage(vtkImageData* mask_before, vtkImageData* mask_after);
  vtkSmartPointer<vtkImageData> GetThresholdedMaskImage();

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


