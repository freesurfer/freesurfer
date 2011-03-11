/**
 * @file  LayerVolumeBase.h
 * @brief Layer data object for MRI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:39 $
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

#ifndef LayerVolumeBase_h
#define LayerVolumeBase_h

#include "LayerEditable.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include <string>
#include <vector>

class vtkImageData;
class BrushProperty;
class LivewireTool;
class Contour2D;

class LayerVolumeBase : public LayerEditable
{
public:
  LayerVolumeBase();
  virtual ~LayerVolumeBase();

  void SetVoxelByRAS( double* ras, int nPlane, bool bAdd = true );
  void SetVoxelByRAS( double* ras1, double* ras2, int nPlane, bool bAdd = true );
  bool FloodFillByRAS( double* ras, int nPlane, bool bAdd = true, char* mask_out = false );

  void SetLiveWireByRAS( double* ras1, double* raw2, int nPlane );
  std::vector<double> GetLiveWirePointsByRAS( double* pt1, double* pt2, int nPlane );

  bool HasUndo();
  bool HasRedo();

  void Undo();
  void Redo();

  void Copy( int nPlane );
  void Paste( int nPlane );
  bool CopyStructure( int nPlane, double* ras );

  void SaveForUndo( int nPlane );

  float GetFillValue();
  void SetFillValue( float fFill );

  float GetBlankValue();
  void SetBlankValue( float fBlank );

  int GetBrushRadius();
  void SetBrushRadius( int nRadius );

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane ) = 0;
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL ) = 0;

  virtual void SetVisible( bool bVisible = true ) = 0;
  virtual bool IsVisible() = 0;
  
  virtual void UpdateVoxelValueRange( double fValue ) {}

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
  
protected:
  bool SetVoxelByIndex( int* n, int nPlane, bool bAdd = true ); // true is to add, false is to remove
  bool SetVoxelByIndex( int* n1, int* n2, int nPlane, bool bAdd = true );
  bool FloodFillByIndex( int* n, int nPlane, bool bAdd = true, char* mask_out = false );
  bool SetLiveWireByIndex( int* n1, int* n2, int nPlane );

  bool GetConnectedToOld( vtkImageData* img, int nFrame, int* n, int nPlane );

  struct UndoRedoBufferItem
  {
    int  plane;
    int  slice;
    char* data;
  };

  void SaveBufferItem( UndoRedoBufferItem& item, int nPlane, int nSlice, const char* mask = NULL );
  void LoadBufferItem( UndoRedoBufferItem& item, bool bIgnoreZeros = false );

  vtkSmartPointer<vtkImageData> m_imageData;
  vtkSmartPointer<vtkImageData> m_imageDataRef;

  float   m_fFillValue;
  float   m_fBlankValue;

  std::vector<UndoRedoBufferItem>  m_bufferUndo;
  std::vector<UndoRedoBufferItem>  m_bufferRedo;
  UndoRedoBufferItem    m_bufferClipboard;

  int   m_nBrushRadius;

  BrushProperty*  m_propertyBrush;

  LivewireTool*  m_livewire;

  int     m_nActiveFrame;
};

#endif


