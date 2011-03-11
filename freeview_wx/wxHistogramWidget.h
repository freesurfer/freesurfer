/**
 * @file  wxHistogramWidget.h
 * @brief Widget to display histogram data.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:43 $
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

#ifndef wxHistogramWidget_h
#define wxHistogramWidget_h

#include <wx/panel.h>
#include <wx/colour.h>
#include "Broadcaster.h"
#include <vector>

class LineMarker;
typedef std::vector<LineMarker> MarkerArray;

struct LineMarker
{
  enum { Solid = 0, Dash, Dot };
  
  double    position;
  wxColour  color;
  int       style;
};

class WXDLLEXPORT wxHistogramWidget : public wxPanel, public Broadcaster
{
public:
	wxHistogramWidget();
	wxHistogramWidget( wxWindow *parent, 
						wxWindowID id = wxID_ANY,
						const wxPoint& pos = wxDefaultPosition,
						const wxSize& size = wxDefaultSize, 
						long style = wxTAB_TRAVERSAL,
						const wxString& name = wxT( "HistogramWidget" ) );

	virtual ~wxHistogramWidget();

	bool Create( wxWindow *parent, 
				 wxWindowID id = wxID_ANY,
				 const wxPoint& pos = wxDefaultPosition,
				 const wxSize& size = wxDefaultSize, 
				 long style = wxTAB_TRAVERSAL,
				 const wxString& name = wxT( "HistogramWidget" ) );

	template <class T> void SetInputData( T* data, long size );	
    
	void GetOutputRange( double* dRange );
	
	void SetOutputRange( double* dRange );
	
	bool GetAutoRange();
	
	void SetAutoRange( bool bRange );
	
	int GetNumberOfBins();
	
	void SetNumberOfBins( int nBins );
	 
	int GetOutputSize();
	
	void GetOutputData( int* buffer_out );
	
  int GetMaximumCount()
  { 
    return m_nMaxCount; 
  }
  
  void SetForegroundColor( const wxColour& color );
    
  void SetMarkers( const MarkerArray& markers, bool bRefresh = true );
  
  void SetColorTableData( unsigned char* colortable, bool bRefresh = true );
  
protected:
  void Initialize();
	void OnPaint            ( wxPaintEvent& event );
  void OnMouseButtonDown  ( wxMouseEvent& event );
	
	void UpdateData( bool bRepaint = true );
  void UpdateColorTable( );
	
	double*	    m_dInputData;
	long		    m_nInputSize;
	double	    m_dInputRange[2];
	
	int*		    m_nOutputData;
	int			    m_nOutputSize;	
	double		  m_dOutputRange[2];
	bool		    m_bAutoRange;
	int			    m_nNumberOfBins;
  unsigned char* m_nColorTable;       // color table for histogram drawing as RGBA
  
	double      m_dBinWidth;
  int         m_nMaxCount;
    
	wxColour	  m_colorBackground;
  wxColour    m_colorForeground;
  
  wxRect      m_rectGraph;
  
  MarkerArray m_markers;
	
	DECLARE_DYNAMIC_CLASS( wxHistogramWidget )
	DECLARE_EVENT_TABLE()
};

template <class T> void wxHistogramWidget::SetInputData( T* data, long size )
{
	if ( m_dInputData )
		delete[] m_dInputData;
	
	m_dInputData = new double[size];
	if ( m_dInputData )
	{
		m_dInputRange[0] = m_dInputRange[1] = data[0];
		for ( long i = 0; i < size; i++ )
		{
			m_dInputData[i] = data[i];
			if ( m_dInputRange[0] > m_dInputData[i] )
				m_dInputRange[0] = m_dInputData[i];
			else if ( m_dInputRange[1] < m_dInputData[i] )
				m_dInputRange[1] = m_dInputData[i];
		}
		m_nInputSize = size;
    m_dOutputRange[0] = m_dInputRange[0];
    m_dOutputRange[1] = m_dInputRange[1];
	}
    
  UpdateData();
}

#endif


