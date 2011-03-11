/**
 * @file  BrushProperty.h
 * @brief Class to hold brush properties for voxel editing
 *
 * Simpleclass for use with the Listener class so text
 * messages with a pointer data can be sent to a list of listeners.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:35 $
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


#ifndef BrushProperty_h
#define BrushProperty_h

class LayerVolumeBase;

class BrushProperty
{
public:

  BrushProperty ();
  virtual ~BrushProperty ();

  int  GetBrushSize();
  void SetBrushSize( int nSize );

  int  GetBrushTolerance();
  void  SetBrushTolerance( int nTolerance );

  LayerVolumeBase*  GetReferenceLayer();
  void     SetReferenceLayer( LayerVolumeBase* layer );

  double* GetDrawRange();
  void  SetDrawRange( double* range );
  void  SetDrawRange( double low, double high );

  bool GetDrawRangeEnabled();
  void SetDrawRangeEnabled( bool bEnable );

  double* GetExcludeRange();
  void SetExcludeRange( double* range );
  void  SetExcludeRange( double low, double high );

  bool GetExcludeRangeEnabled();
  void SetExcludeRangeEnabled( bool bEnable );

  bool GetDrawConnectedOnly();
  void SetDrawConnectedOnly( bool bEnable );

protected:
  int  m_nBrushSize;
  int  m_nBrushTolerance;
  double m_dDrawRange[2];
  bool m_bEnableDrawRange;
  double m_dExcludeRange[2];
  bool m_bEnableExcludeRange;
  bool m_bDrawConnectedOnly;

  LayerVolumeBase* m_layerRef;
};


#endif
