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

#ifndef LayerVolumeBase_h
#define LayerVolumeBase_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include <vector>
#include <QFile>
#include <QVariantMap>
#include <QVector>

class vtkImageData;
class BrushProperty;
class LivewireTool;

class LayerVolumeBase : public LayerEditable
{
  Q_OBJECT
public:
  LayerVolumeBase( QObject* parent = NULL );
  virtual ~LayerVolumeBase();

  void SetVoxelByRAS( double* ras, int nPlane, bool bAdd = true, bool ignore_brush_size = false );
  void SetVoxelByRAS( double* ras1, double* ras2, int nPlane, bool bAdd = true, bool ignore_brush_size = false, bool draw_as_box = false  );
  bool FloodFillByRAS( double* ras, int nPlane, bool bAdd = true, bool b3D = false, char* mask_out = 0, bool ignore_exclusion = false );
  void CloneVoxelByRAS( double* ras, int nPlane );
  void CloneVoxelByRAS( double* ras1, double* ras2, int nPlane );
  void ShiftVoxelsByRAS( double* ras_offset, int nPlane);
  void ShiftVoxels(int *nOffset, int nPlane);

  bool BorderFillByRAS( double* ras, int nPlane, bool b3D = false);

  void SetLiveWireByRAS( double* ras1, double* raw2, int nPlane );
  std::vector<double> GetLiveWirePointsByRAS( double* pt1, double* pt2, int nPlane );

  virtual bool HasUndo();
  virtual bool HasRedo();

  virtual void Undo();
  virtual void Redo();

  void Copy( int nPlane );
  void Paste( int nPlane );
  bool CopyStructure( int nPlane, double* ras );

  virtual void SaveForUndo( int nPlane = -1, bool bAllFrames = false );

  double GetFillValue();

  double GetBlankValue();

  int GetBrushRadius();

  virtual void UpdateVoxelValueRange( double fValue ) { Q_UNUSED(fValue); }

  virtual int GetActiveFrame()
  {
    return m_nActiveFrame;
  }

  virtual void SetActiveFrame( int nFrame )
  {
    m_nActiveFrame = nFrame;
  }

  vtkImageData* GetImageData()
  {
    return m_imageData;
  }

  vtkImageData* GetImageDataRef()
  {
    return m_imageDataRef;
  }

  void SetBrushProperty( BrushProperty* brush )
  {
    m_propertyBrush = brush;
  }

  bool IsValidToPaste( int nPlane );

  double GetMinimumVoxelSize();

  virtual void GetDisplayBounds( double* bounds );

signals:
  void FillValueChanged( double );
  void EraseValueChanged( double );
  void BrushRadiusChanged( int );
  void BaseVoxelEdited(const QVector<int>, bool bAdd);

public slots:
  void SetFillValue( double fFill );
  void SetBlankValue( double fBlank );
  void SetBrushRadius( int nRadius );
  void ClearVoxels();
  void PrepareShifting(int nPlane);
  void DoneShifting();

protected:
  QVector<int> SetVoxelByIndex( int* n, int nPlane, bool bAdd = true, bool ignore_brush_size = false ); // true is to add, false is to remove
  QVector<int> SetVoxelByIndex( int* n1, int* n2, int nPlane, bool bAdd = true, bool ignore_brush_size = false, bool draw_as_box = false  );
  QVector<int> FloodFillByIndex( int* n, int nPlane, bool bAdd = true, bool ignore_overflow = true, char* mask_out = NULL, bool ignore_exclusion = false );
  QVector<int> BorderFillByRAS(int *n, int nPlane);
  void FillBoxByIndex( int* n1, int* n2, int nPlane);
  bool SetLiveWireByIndex( int* n1, int* n2, int nPlane );
  bool CloneVoxelByIndex( int* n, int nPlane );
  bool CloneVoxelByIndex( int* n1, int* n2, int nPlane );

  bool GetConnectedToOld( vtkImageData* img, int nFrame, int* n, int nPlane );

  struct UndoRedoBufferItem
  {
    UndoRedoBufferItem()
    {
      data = 0;
    }
    void Clear()
    {
      if (data)
      {
        delete[] data;
      }
      data = 0;

      if (!cache_filename.isEmpty())
      {
        if (QFile::exists(cache_filename))
        {
          QFile::remove(cache_filename);
        }
      }
    }

    int  plane;                 // -1 means whole 3d volume
    int  slice;
    char* data;
    int  frame;
    QString cache_filename;     // if not empty, ignore data and read from cache file.
    QVariantMap mri_settings;
  };

  void SaveBufferItem( UndoRedoBufferItem& item, int nPlane = -1, int nSlice = 0, int nFrame = 0, const char* mask = NULL );
  void LoadBufferItem( UndoRedoBufferItem& item, bool bIgnoreZeros = false );
  QString GenerateCacheFileName();

  vtkSmartPointer<vtkImageData> m_imageData;
  vtkSmartPointer<vtkImageData> m_imageDataRef;

  double   m_fFillValue;
  double   m_fBlankValue;

  std::vector<UndoRedoBufferItem>  m_bufferUndo;
  std::vector<UndoRedoBufferItem>  m_bufferRedo;
  UndoRedoBufferItem    m_bufferClipboard;

  int   m_nBrushRadius;

  BrushProperty*  m_propertyBrush;

  LivewireTool*  m_livewire;

  int     m_nActiveFrame;

  char*   m_shiftBackgroundData;
  char*   m_shiftForegroundData;
};

#endif


