/**
 * @file  ConnectivityData.h
 * @brief Holder for connectivity data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:36 $
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

#ifndef ConnectivityData_h
#define ConnectivityData_h

#include "Listener.h"
#include "Broadcaster.h"
#include <vtkSmartPointer.h>
#include <wx/string.h>

extern "C"
{
#include "mri.h"
#include "colortab.h"
}


class LayerSurface;
class vtkPolyData;
class vtkActor;
class vtkRGBAColorTransferFunction;
class vtkRenderer;

class ConnectivityData : public Listener, public Broadcaster
{
  public:
    enum DISPLAY_MODE { DM_Single, DM_Both, DM_All };
    
    ConnectivityData();
    virtual ~ConnectivityData();
    
    bool IsValid();
    
    bool Load( const char* filename, const char* lut_fn, LayerSurface* surf1, LayerSurface* surf2 );

    void SetVaryRadius  ( bool bVary );
    bool GetVaryRadius  ()
    {
      return m_bVaryRadius;
    }    
    
    void SetRadius      ( double dRadius );
    void SetRadiusMax   ( double dRadius );
    void SetRadiusRange ( double dMin, double dMax );
    void GetRadiusRange ( double* range )
    {
      range[0] = m_dRadius;
      range[1] = m_dRadiusMax;
    }
    
    void SetThresholdMin  ( double dMin );
    void SetThresholdMax  ( double dMax );
    void SetThresholdRange( double dMin, double dMax );
    void GetThresholdRange( double* range )
    {
      range[0] = m_dThresholdMin;
      range[1] = m_dThresholdMax;
    }
    
    void SetDisplayMode ( int nMode );
    int  GetDisplayMode ()
    {
      return m_nDisplayMode;
    }
    
    void SetIncrementalDisplay( bool bIncre );
    bool GetIncrementalDisplay()
    {
      return m_bIncrementalDisplay;
    }
    
    void GetDataRange( double* range )
    {
      range[0] = m_dMinConnectivity;
      range[1] = m_dMaxConnectivity;
    }
    
    void BuildConnectivityActors();
    
    void AppendProps( vtkRenderer* renderer );
    
    bool Export( const char* filename );
    
    bool HasAnySeeds();
    
  protected:
    int GetColorTableIndex( LayerSurface* surf, int annotation_index );
    virtual void DoListenToMessage ( std::string const iMsg, void* iData, void* sender );
    void BuildWireActor( int seed_index, double* dStartPos, double* dPos, double val, vtkPolyData* polydata ); 
    
    // properties
    bool          m_bVaryRadius;
    double        m_dRadius;
    double        m_dRadiusMax;
    double        m_dThresholdMin;
    double        m_dThresholdMax;
    int           m_nDisplayMode;
    int           m_nLastDisplayMode;
    bool          m_bIncrementalDisplay;
    
    // data
    MRI*          m_MRI;
    COLOR_TABLE*  m_lut;
    LayerSurface* m_surfLeft;
    LayerSurface* m_surfRight;
    
    vtkSmartPointer<vtkRGBAColorTransferFunction> m_vtkLut;
    vtkSmartPointer<vtkActor>     m_actorStartPoint;
    vtkSmartPointer<vtkActor>     m_actorTube;

    double**      m_dMatConnectivity;
    double        m_dMinConnectivity;
    double        m_dMaxConnectivity;
    int*          m_nSeedLabels;
    int           m_nMatSize;
    
    wxString      m_name;
};

#endif
