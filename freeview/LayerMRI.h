/**
 * @brief Layer class for MRI volume.
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

#ifndef LayerMRI_h
#define LayerMRI_h

#include "LayerVolumeBase.h"
#include "vtkSmartPointer.h"
#include <QString>
#include <QList>



#include "colortab.h"
#include "nifti1.h"


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
class LayerSurface;
class LayerROI;
class GeoSWorker;
class Region3D;
class LayerPointSet;

#ifndef IntList
typedef QList<int> IntList;
#endif

class LayerMRI : public LayerVolumeBase
{
  friend class ThreadBuildContour;
  friend class VolumeCropper;
  friend class SurfaceRegionGroups;
  friend class LayerMRIWorkerThread;

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
  bool CreateFromMRIData(void* mri);  // hack
  bool LoadVolumeTransform();
  void UnloadVolumeTransform();

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane );
  virtual void Append3DProps( vtkRenderer* renderer, bool* bPlaneVisibility = NULL );
  virtual bool HasProp( vtkProp* prop );

  void Remove2DProps( vtkRenderer* render, int nPlane );

  //  void SetSliceNumber( int* sliceNumber );
  void SetSlicePositionToWorldCenter();

  virtual double GetVoxelValue( double* pos );
  double GetVoxelValueByOriginalIndex( int i, int j, int k, int frame = -1 );
  QList<double> GetVoxelValueByOriginalIndexAllFrames(int i, int j, int k);
  void GetVoxelValueByOriginalIndexAllFrames(int i, int j, int k, float* buffer);
  double GetSampledVoxelValueByRAS(double* ras, int frame = -1);

  std::vector<double> GetSampledVoxelValues(std::vector< std::vector<double> >& line3d, int frame = -1);
  std::vector<double> GetMeanSegmentValues(std::vector< std::vector<double> >& line3d, int frame = -1);

  virtual QString GetLabelName( double value );

  void RASToOriginalIndex( const double* pos, int* n );
  void RASToOriginalIndex(const double *pos, double *n_out);
  void OriginalIndexToRAS( const int* n, double* pos );
  void OriginalVoxelToRAS( const double* vcoord, double* pos);
  void TargetIndexToOriginalIndex(const int* n_in, int* n_out);

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

  virtual int GetNumberOfFrames();

  void GetRASCenter( double* pt );

  bool IsTransformed();

  void SetReorient( bool bReorient );

  void SetSampleMethod( int nSampleMethod );

  void SetConform( bool bConform );

  bool GetVoxelValueRange( const double* pt0, const double* pt1,
                           int nPlane, double* range_out );

  bool GetVoxelStatsRectangle( const double* pt0, const double* pt1,
                               int nPlane, double* mean_out, double* sd_out = NULL, int* cnt_out = NULL );

  bool GetVoxelsOnLine( const double* pt0, const double* pt1, int nPlane, int*& indice_out, double*& value_out, int* cnt_out );

  bool GetVoxelStats(QVector<int>& indices, double* mean_out, double* sd_out = NULL);

  bool GetVoxelStatsByTargetRAS(QVector<float> &coords, double* mean_out, double *sd_out = NULL);

  void ResetWindowLevel();

  int GetDataType();

  virtual COLOR_TABLE* GetEmbeddedColorTable();

  void SnapToVoxelCenter( const double* pt_in, double* pt_out );

  void GetCurrentLabelStats( int nPlane, float* label_out, int* count_out, float* area_out,
                             LayerMRI* underlying_mri = NULL, double* mean_out = NULL, double* sd_out = NULL );

  vtkImageData* GetSliceImageData( int nPlane );

  bool FloodFillByContour2D( double* ras, Contour2D* c2d );

  bool SaveContourToFile(const QString& fn);

  SurfaceRegion* CreateNewSurfaceRegion( double* pt, vtkProp* prop );

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

  void SetMaskLayer(LayerMRI* layer_mask);

  double GetMaskThreshold()
  {
    return m_dMaskThreshold;
  }

  LayerMRI* GetMaskLayer()
  {
    return m_layerMask;
  }

  void SetCorrelationSurface(LayerSurface* surf);

  LayerSurface* GetCorrelationSurface()
  {
    return m_correlationSurface;
  }

  double GetHistoValueFromPercentile(double percentile);

  double GetHistoPercentileFromValue(double value);

  bool HasValidHistogram();

  void Threshold(int frame, LayerMRI* src, int src_frame, double th_low, double th_high, bool replace_in, double in_value, bool replace_out, double out_value);

  bool Segment(int min_label_index, int max_label_index, int min_num_of_voxels);

  void RestoreFromBackup();

  bool GetLayerLabelCenter(double val, double* pos_out);
  
  bool IsWindowAdjustable();
  
  bool IsObscuring();

  bool GeodesicSegmentation(LayerMRI* seeds, double lambda, int wsize, double max_dist, double smoothing_std, LayerMRI* mask, double max_foreground_dist = 0);

  void GeodesicSegmentationAbort();

  void GeodesicSegmentationApply(LayerMRI* filled);

  void GetVolumeInfo(int* dim, double* voxel_size);

  void SetIgnoreHeader(bool b)
  {
    m_bIgnoreHeader = b;
  }

  QVector<double> GetVoxelList(int nVal, bool bForce = false);

  int GetLabelCount(int nVal);

  QVariantMap GetTimeSeriesInfo();

  QString GetGeoSegErrorMessage();

  bool ExportLabelStats(const QString& fn);

  QList<vtkActor*> GetContourActors(bool bVisibleOnly = false);

  Region3D* GetCurrent3DRegion()
  {
    return m_current3DRegion;
  }
  
  Region3D* CreateNew3DRegion( double* pt, vtkProp* prop );

  bool DeleteCurrent3DRegion();

  void DeleteAll3DRegions();

  void Add3DRegionPoint( double* pt );

  int GetNumberOf3DRegions()
  {
    return m_3DRegions.size();
  }

  Region3D* Select3DRegion( double* pos, double dist);

  bool SaveAll3DRegions(const QString& fn);

  bool Load3DRegions(const QString& fn);

  void Close3DRegion();

  void UpdateVoxelsByPointSet(LayerPointSet* ps, int nPlane);

  void LocateLocalMaximumAtRAS(double* ras_in, double dx, double dy, double dz, double* ras_out, double dist_in_vox = 3);

public slots:
  virtual void SetModified();
  void SetActiveFrame( int nFrame );
  void SetActiveFrameOneBase( int nFrame )
  {
    SetActiveFrame( nFrame-1 );
  }
  void UpdateResliceInterpolation();
  void UpdateMRIToImage();
  void SetMaskThreshold(double val);

Q_SIGNALS:
  void ResampleFactorChanged();
  //  void ColorMapChanged();
  void ActiveFrameChanged( int nFrame );
  void SurfaceRegionAdded();
  void SurfaceRegionUpdated();
  void SurfaceRegionRemoved();
  void IsoSurfaceUpdating();
  void IsoSurfaceUpdated();
  void LabelStatsReady();
  void CorrelationSurfaceChanged(LayerSurface*);
  void GeodesicSegmentationApplied();
  void GeodesicSegmentationFinished(double time_in_secs);
  void GeodesicSegmentationProgress(double percentage);
  void Region3DAdded();
  void Region3DRemoved();

protected slots:
  virtual void UpdateDisplayMode();
  virtual void UpdateOpacity();
  void UpdateTextureSmoothing();
  void UpdateContour( int nSegIndex = -1 );
  void UpdateContourActor( int nSegIndex );
  void UpdateContourColor();
  void ShowContour();
  void UpdateVolumeRendering();
  void UpdateVectorActor();
  void UpdateVectorActor( int nPlane, vtkImageData* imagedata, vtkImageData* scaledata = NULL );
  virtual void UpdateVectorActor( int nPlane );

  void ResetSurfaceRegionIds();

  void UpdateLabelOutline();
  void UpdateUpSampleMethod();
  void UpdateProjectionMap();

  void UpdateTensorActor();
  virtual void UpdateColorMap();

  void OnContourThreadFinished(int thread_id);
  void UpdateSurfaceCorrelationData();

  void ResetRef();

  void OnLabelContourChanged(int n = -1);
  void RebuildContour();

  void UpdateLabelInformation();
  void OnLabelInformationReady();

  void UpdateVectorLineWidth(double val);

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
  void GetColorWheelColor(double* v, int plane, unsigned char* c_out);
  void UpdateNiftiHeader();

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
  vtkSmartPointer<vtkImageReslice>     mResample[3];

  FSVolume*   m_volumeSource;
  FSVolume*   m_volumeRef;
  bool    m_bResampleToRAS;
  bool    m_bReorient;
  int     m_nSampleMethod;
  bool    m_bConform;
  bool    m_bWriteResampled;
  bool    m_bIgnoreHeader;

  vtkImageActor*  m_sliceActor2D[3];
  vtkImageActor*  m_sliceActor3D[3];

  vtkActor*       m_glyphActor2D[3];
  vtkActor*       m_vectorDotActor2D[3];
  vtkActor*       m_glyphActor3D[3];

  vtkImageActor*  m_projectionMapActor[3];

  vtkSmartPointer<vtkActor>   m_actorContour;
  vtkSmartPointer<vtkVolume>  m_propVolume;
  QMap<int, vtkActor*>            m_labelActors;
  vtkActor*                   m_actorCurrentContour;

  int         m_nThreadID;
  vtkSmartPointer<vtkActor>       m_actorContourTemp;
  QMap<int, vtkActor*>            m_labelActorsTemp;

  QList<SurfaceRegion*>           m_surfaceRegions;
  SurfaceRegion*                  m_currentSurfaceRegion;
  SurfaceRegionGroups*            m_surfaceRegionGroups;

  QList<Region3D*>            m_3DRegions;
  Region3D*                   m_current3DRegion;

  int         m_nOrientationIndex[3];

  int         m_nGotoLabelSlice;
  int         m_nGotoLabelOrientation;
  QString     m_strGotoLabelName;

  vtkSmartPointer<vtkImageData> m_imageDataBackup;
  LayerMRI*   m_layerMask;

  vtkSmartPointer<vtkImageData> m_imageRawDisplay;
  LayerSurface* m_correlationSurface;

private:
  double**    private_buf1_3x3;
  double**    private_buf2_3x3;

  LayerMRIWorkerThread* m_worker;
  QList<int>  m_nAvailableLabels;
  QMap<int, QList<double> > m_listLabelCenters;
  QMap<int, QVector<double> > m_voxelLists;
  QMap<int, int> m_labelVoxelCounts;

  QMap<QObject*, double>  m_mapMaskThresholds;
  double      m_dMaskThreshold;

  nifti_1_header    m_niftiHeader;

  GeoSWorker* m_geos;
};


#endif


