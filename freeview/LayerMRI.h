/**
 * @file  LayerMRI.h
 * @brief Layer class for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/24 19:59:46 $
 *    $Revision: 1.81 $
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

#ifndef LayerMRI_h
#define LayerMRI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"
#include <QString>
#include <QList>

extern "C"
{
#include "colortab.h"
}

class vtkImageReslice;
class vtkImageMapToColors;
class vtkSimpleLabelEdgeFilter;
class vtkTransform;
class vtkTexture;
class vtkPolyDataMapper;
class vtkActor;
class vtkImageActor;
class vtkImageData;
class vtkProp;
class vtkVolume;
class vtkPolyData;
class vtkPolyDataAlgorithm;
class vtkImageResample;
class vtkUnsignedCharArray;
class LayerPropertyMRI;
class FSVolume;
class BuildContourThread;
class Contour2D;
class SurfaceRegion;
class SurfaceRegionGroups;
class LayerMRIWorkerThread;

#ifndef IntList
typedef QList<int> IntList;
#endif

class LayerMRI : public LayerVolumeBase
{
  friend class ThreadBuildContour;
  friend class VolumeCropper;
  friend class SurfaceRegionGroups;

  Q_OBJECT
public:
  LayerMRI( LayerMRI* ref, QObject* parent = NULL );
  virtual ~LayerMRI();

  inline LayerPropertyMRI* GetProperty()
  {
    return (LayerPropertyMRI*)mProperty;
  }

  bool LoadVolumeFromFile();
  bool Create( LayerMRI* mri, bool bCopyVoxel, int data_type = -1, int voxel_option = -1 );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );
  virtual bool HasProp( vtkProp* prop );

  void Remove2DProps( vtkRenderer* render, int nPlane );

//  void SetSliceNumber( int* sliceNumber );
  void SetSlicePositionToWorldCenter();

  virtual double GetVoxelValue( double* pos );
  double GetVoxelValueByOriginalIndex( int i, int j, int k, int frame = -1 );
  QList<double> GetVoxelValueByOriginalIndexAllFrames(int i, int j, int k);
  double GetSampledVoxelValueByRAS(double* ras, int frame = -1);

  std::vector<double> GetSampledVoxelValues(std::vector< std::vector<double> >& line3d, int frame = -1);
  std::vector<double> GetMeanSegmentValues(std::vector< std::vector<double> >& line3d, int frame = -1);

  virtual QString GetLabelName( double value );

  void RASToOriginalIndex( const double* pos, int* n );
  void OriginalIndexToRAS( const int* n, double* pos );

  virtual void SetVisible( bool bVisible = true );
  virtual bool IsVisible();

  virtual void UpdateVoxelValueRange( double dValue );

  FSVolume* GetSourceVolume()
  {
    return m_volumeSource;
  }

  FSVolume* GetRefVolume()
  {
    return m_volumeRef;
  }

  void SetRefVolume( FSVolume* ref )
  {
    m_volumeRef = ref;
  }

  bool SaveVolume();

  void SetResampleToRAS( bool bResample );

  bool GetResampleToRAS()
  {
    return m_bResampleToRAS;
  }

  void RemapPositionToRealRAS( const double* pos_in, double* pos_out );

  void RemapPositionToRealRAS( double x_in, double y_in, double z_in,
                               double& x_out, double& y_out, double& z_out );

  void TargetToRAS( const double* pos_in, double* pos_out )
  {
    RemapPositionToRealRAS( pos_in, pos_out );
  }

  void RASToTarget( const double* pos_in, double* pos_out );

  void NativeRASToTkReg( const double* pos_in, double* pos_out );

  void TkRegToNativeRAS( const double* pos_in, double* pos_out );

  int GetNumberOfFrames();

  void GetRASCenter( double* pt );

  bool IsTransformed();

  void SetReorient( bool bReorient );

  void SetSampleMethod( int nSampleMethod )
  {
    m_nSampleMethod = nSampleMethod;
  }

  void SetConform( bool bConform );

  bool GetVoxelValueRange( const double* pt0, const double* pt1,
                           int nPlane, double* range_out );

  bool GetVoxelStatsRectangle( const double* pt0, const double* pt1,
                               int nPlane, double* mean_out, double* sd_out = NULL, int* cnt_out = NULL );

  bool GetVoxelsOnLine( const double* pt0, const double* pt1, int nPlane, int*& indice_out, double*& value_out, int* cnt_out );

  bool GetVoxelStats(QList<int>& indices, double* mean_out, double* sd_out = NULL);

  void ResetWindowLevel();

  int GetDataType();

  virtual COLOR_TABLE* GetEmbeddedColorTable();

  void SnapToVoxelCenter( const double* pt_in, double* pt_out );

  void GetCurrentLabelStats( int nPlane, float* label_out, int* count_out, float* area_out,
                             LayerMRI* underlying_mri = NULL, double* mean_out = NULL, double* sd_out = NULL );

  vtkImageData* GetSliceImageData( int nPlane );

  bool FloodFillByContour2D( double* ras, Contour2D* c2d );

  virtual void SetModified();

  bool SaveContourToFile(const QString& fn);

  SurfaceRegion* CreateNewSurfaceRegion( double* pt );

  void AddSurfaceRegionLoopPoint( double* pt );

  void CloseSurfaceRegion();

  SurfaceRegion* SelectSurfaceRegion( double* pos );
  SurfaceRegion* SelectSurfaceRegion( int nId );

  SurfaceRegion* GetCurrentSurfaceRegion()
  {
    return m_currentSurfaceRegion;
  }

  bool DeleteCurrentSurfaceRegion();

  int GetNumberOfSurfaceRegions()
  {
    return m_surfaceRegions.size();
  }

  bool SaveAllSurfaceRegions( const QString& fn );

  bool LoadSurfaceRegions( const QString& fn );

  QString GetOrientationString();

  void SetCropToOriginal(bool bCropToOriginal);

  void SetCroppingBounds( double* bounds );

  virtual void GetDisplayBounds( double* bounds );

  bool SaveRegistration( const QString& filename );

  void GetLabelStats( LayerMRI* label, int nPlane,
                      std::vector<int>& id,
                      std::vector<int>& number,
                      std::vector<double>& mean,
                      std::vector<double>& std );

  SurfaceRegionGroups* GetSurfaceRegionGroups()
  {
    return m_surfaceRegionGroups;
  }

  void SetWriteResampled(bool bResample)
  {
    m_bWriteResampled = bResample;
  }

  bool GetWriteResampled()
  {
    return m_bWriteResampled;
  }

  int GoToLabel(int orientation, const QString& label_name);

  void ReplaceVoxelValue(double orig_value, double new_value, int nPlane = -1);

  void SetGotoLabel(int nOrientation, const QString& name)
  {
    m_nGotoLabelOrientation = nOrientation;
    m_strGotoLabelName = name;
  }

  int GetGotoLabelSlice()
  {
    return m_nGotoLabelSlice;
  }

  double GetTR();

  QList<int> GetAvailableLabels()
  {
    return m_nAvailableLabels;
  }

  bool SaveIsoSurface(const QString& fn);

  bool HasReg();

public slots:
  void SetActiveFrame( int nFrame );
  void SetActiveFrameOneBase( int nFrame )
  {
    SetActiveFrame( nFrame-1 );
  }

Q_SIGNALS:
  void ResampleFactorChanged();
//  void ColorMapChanged();
  void ActiveFrameChanged( int nFrame );
  void SurfaceRegionAdded();
  void SurfaceRegionUpdated();
  void SurfaceRegionRemoved();
  void IsoSurfaceUpdated();
  void LabelStatsReady();

protected slots:
  void UpdateDisplayMode();
  virtual void UpdateOpacity();
  void UpdateResliceInterpolation();
  void UpdateTextureSmoothing();
  void UpdateContour( int nSegIndex = -1 );
  void UpdateContourActor( int nSegIndex );
  void UpdateContourColor();
  void ShowContour();
  void UpdateVolumeRendering();
  void UpdateVectorActor();
  void UpdateVectorActor( int nPlane, vtkImageData* imagedata );
  virtual void UpdateVectorActor( int nPlane );

  void ResetSurfaceRegionIds();

  void UpdateLabelOutline();
  void UpdateUpSampleMethod();
  void UpdateProjectionMap();

  void UpdateTensorActor();
  virtual void UpdateColorMap();

  void OnContourThreadFinished(int thread_id);

  void OnAvailableLabels(const IntList& vals);

protected:
  virtual void DoTransform(double *mat, int sample_method);
  virtual bool DoRotate( std::vector<RotationElement>& rotations );
  virtual void DoTranslate( double* offset );
  virtual void DoScale( double* rscale, int nSampleMethod );
  virtual void DoRestore();
  virtual void DoTransform(int sample_method);

  void InitializeVolume();
  void InitializeActors();
  void ConnectProperty();
  void UpdateTensorActor( int nPlane, vtkImageData* imagedata = NULL );

  std::vector<int> GetVoxelIndicesBetweenPoints( int* n0, int* n1 );
  void BuildTensorGlyph( vtkImageData* imagedata,
                         int i, int j, int k,
                         double* pt, double scale,
                         vtkPolyData* sourcepolydata,
                         vtkUnsignedCharArray* scalars,
                         vtkPolyDataAlgorithm* a);


  virtual void OnSlicePositionChanged( int nPlane );


  // Pipeline ------------------------------------------------------------
  vtkSmartPointer<vtkImageReslice>      mReslice[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMap[3];
  vtkSmartPointer<vtkImageMapToColors>  mColorMapMaxProjection[3];
  vtkSmartPointer<vtkSimpleLabelEdgeFilter>   mEdgeFilter[3];
  vtkSmartPointer<vtkImageResample>     mResample[3];

  FSVolume*   m_volumeSource;
  FSVolume*   m_volumeRef;
  bool    m_bResampleToRAS;
  bool    m_bReorient;
  int     m_nSampleMethod;
  bool    m_bConform;
  bool    m_bWriteResampled;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];

  vtkActor*       m_glyphActor2D[3];
  vtkActor*       m_glyphActor3D[3];

  vtkImageActor*  m_projectionMapActor[3];

  struct SegmentationActor
  {
    int id;
    vtkActor* actor;
  };

  QList<SegmentationActor>    m_segActors;

  vtkSmartPointer<vtkActor>   m_actorContour;
  vtkSmartPointer<vtkVolume>  m_propVolume;

  int         m_nThreadID;
  vtkSmartPointer<vtkActor>       m_actorContourTemp;

  QList<SurfaceRegion*>           m_surfaceRegions;
  SurfaceRegion*                  m_currentSurfaceRegion;
  SurfaceRegionGroups*            m_surfaceRegionGroups;

  int         m_nOrientationIndex[3];

  int         m_nGotoLabelSlice;
  int         m_nGotoLabelOrientation;
  QString     m_strGotoLabelName;

private:
  double**    private_buf1_3x3;
  double**    private_buf2_3x3;

  LayerMRIWorkerThread* m_worker;
  QList<int>  m_nAvailableLabels;
};


#endif


