/**
 * @file  LayerEditable.h
 * @brief Base Layer class for editable volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.13 $
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

#ifndef LayerEditable_h
#define LayerEditable_h

#include "Layer.h"
#include "vtkSmartPointer.h"
#include <string>

class LayerEditable : public Layer
{
public:
  LayerEditable();
  virtual ~LayerEditable();

  virtual bool HasUndo()
  {
    return false;
  }
  virtual bool HasRedo()
  {
    return false;
  }

  virtual void Undo()
  {}
  virtual void Redo()
  {}

  virtual void Append2DProps( vtkRenderer* renderer, int nPlane ) = 0;
  virtual void Append3DProps( vtkRenderer* renderer, bool* bSliceVisibility = NULL ) = 0;

  virtual void SetVisible( bool bVisible = true ) = 0;
  virtual bool IsVisible() = 0;

  void ResetModified()
  {
    m_bModified = false;
  }

  bool IsModified()
  {
    return m_bModified;
  }

  void SetEditable( bool bEditable = true )
  {
    m_bEditable = bEditable;
  }

  bool IsEditable()
  {
    return m_bEditable;
  }

  const char* GetFileName()
  {
    return m_sFilename.c_str();
  }

  void SetFileName( const char* fn )
  {
    m_sFilename = fn;
  }

  const char* GetRegFileName()
  {
    return m_sRegFilename.c_str();
  }

  void SetRegFileName( const char* fn )
  {
    m_sRegFilename = fn;
  }

  virtual void SetModified();
  
protected:

  int   m_nMaxUndoSteps;
  bool  m_bModified;

  bool   m_bEditable;

  std::string m_sFilename;
  std::string m_sRegFilename;
};

#endif


