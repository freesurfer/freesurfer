/**
 * @file  FSWayPoints.h
 * @brief Base way points class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:37 $
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


#include <wx/wx.h>
#include <wx/ffile.h>
#include "FSWayPoints.h"
#include "FSVolume.h"
#include <stdexcept>
#include "vtkImageData.h"
#include "MyUtils.h"
#include <vector>

using namespace std;

FSWayPoints::FSWayPoints() :
    m_label( NULL )
{}

FSWayPoints::~FSWayPoints()
{
  if ( m_label )
    ::LabelFree( &m_label );
}

bool FSWayPoints::IsLabelFormat( const char* filename )
{
  char* fn = strdup( filename );
  LABEL* label = ::LabelRead( NULL, fn );
  free( fn );

  if ( label == NULL )
  {
    return false;
  }
  else
  {
    ::LabelFree( &label );
    return true;
  }
}

bool FSWayPoints::ReadAsLabel( const char* filename )
{
  if ( m_label )
    ::LabelFree( &m_label );

  char* fn = strdup( filename );
  m_label = ::LabelRead( NULL, fn );
  free( fn );

  if ( m_label == NULL )
  {
    cerr << "LabelRead failed";
    return false;
  }

  return true;
}

bool FSWayPoints::ReadAsControlPoints( const char* filename )
{
  if ( m_label )
    ::LabelFree( &m_label );
  
  wxFFile file( wxString::FromAscii(filename) );
  wxString strg;
  if ( !file.ReadAll( &strg ) )
    return false;
  
  wxArrayString ar = MyUtils::SplitString( strg, _("\n") );
  int nCount = 0;
  double dvalue;
  vector<float> values;
  for ( size_t i = 0; i < ar.size(); i++ )
  {
    wxArrayString subs = MyUtils::SplitString( ar[i], _(" ") );
    if ( subs.size() > 2 )
    {
      for ( int j = 0; j < 3; j++ )
      {
        subs[j].ToDouble( &dvalue );
        values.push_back( dvalue );
      }
      nCount ++;
    }
  }

  m_label = ::LabelAlloc( nCount, NULL, (char*)"" );
  m_label->n_points = nCount;
  for ( int i = 0; i < nCount; i++ )
  {
    m_label->lv[i].x = values[i*3];
    m_label->lv[i].y = values[i*3+1];
    m_label->lv[i].z = values[i*3+2];
    m_label->lv[i].vno = -1;
    m_label->lv[i].deleted = false;
    m_label->lv[i].stat = 0;
  }

  return true;
}


bool FSWayPoints::WriteAsLabel( const char* filename )
{
  char* fn = strdup( filename );
  int err = ::LabelWrite( m_label, fn );
  free( fn );

  if ( err != 0 )
  {
    cerr << "Way Points Write failed" << endl;
  }

  return err == 0;
}


bool FSWayPoints::WriteAsControlPoints( const char* filename )
{
  wxFFile file( wxString::FromAscii(filename), "w" );
  wxString strg;
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    strg.Printf( wxT("%f %f %f\n"), m_label->lv[i].x, m_label->lv[i].y, m_label->lv[i].z );
    if ( !file.Write( strg ) )
    {
      cerr << "Control Points Write failed" << endl;
      return false;
    }
  }
  strg.Printf( wxT("info\nnumpoints %d\nuseRealRAS 1\n"), m_label->n_points );
  file.Write( strg );
  
  return true;
}


void FSWayPoints::UpdateLabel( std::vector<WayPoint>& points_in, FSVolume* ref_vol )
{
  if ( m_label )
    ::LabelFree( &m_label );

  int nCount = 0;
  double pos[3];
  vector<float> values;
  for ( size_t i = 0; i < points_in.size(); i++ )
  {
    ref_vol->TargetToRAS( points_in[i].pt, pos );
    ref_vol->RASToNativeRAS( pos, pos );
    values.push_back( pos[0] );
    values.push_back( pos[1] );
    values.push_back( pos[2] );
    values.push_back( points_in[i].value );
    nCount ++;
  }

  m_label = ::LabelAlloc( nCount, NULL, (char*)"" );
  m_label->n_points = nCount;
  for ( int i = 0; i < nCount; i++ )
  {
    m_label->lv[i].x = values[i*4];
    m_label->lv[i].y = values[i*4+1];
    m_label->lv[i].z = values[i*4+2];
    m_label->lv[i].vno = -1;
    m_label->lv[i].deleted = false;
    m_label->lv[i].stat = values[i*4+3];
  }
}

void FSWayPoints::LabelToWayPoints( std::vector<WayPoint>& points_out, FSVolume* ref_vol )
{
  if ( !m_label )
  {
    cerr << "Label is empty" << endl;
    return;
  }

  WayPoint wp;
  for ( int i = 0; i < m_label->n_points; i++ )
  {
    wp.pt[0] = m_label->lv[i].x;
    wp.pt[1] = m_label->lv[i].y;
    wp.pt[2] = m_label->lv[i].z;
    wp.value = m_label->lv[i].stat;
    ref_vol->NativeRASToRAS( wp.pt, wp.pt );
    ref_vol->RASToTarget( wp.pt, wp.pt );
    points_out.push_back( wp );
  }
}
