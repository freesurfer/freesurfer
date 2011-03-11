/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
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
#include <wx/filename.h>
#include <iostream>
#include "LUTDataHolder.h"
#include <iostream>
#include <wx/xrc/xmlres.h>
#include <wx/file.h>
#include "res/FreeSurferColorLUT.h"

LUTDataHolder::LUTDataHolder()
{
  ColorTableData ctd;  
  wxFileName fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/FreeSurferColorLUT.txt") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
  ctd.name = "FreeSurferColorLUT";
  if ( ctd.table )
    m_tables.push_back( ctd );

  fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/tkmeditParcColorsCMA") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
  ctd.name = "tkmeditParcColorsCMA";
  if ( ctd.table )
    m_tables.push_back( ctd );

  fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/Simple_surface_labels2009.txt") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
  ctd.name = "Simple_surface_labels2009";
  if ( ctd.table )
  m_tables.push_back( ctd );
  
  // did not find any stock tables, load build-in one
  if ( m_tables.size() == 0)
  {

    wxString tempfn = wxFileName::GetTempDir() + "/FreeSurferColorLUT.txt";
    wxFile file;
    if ( file.Open( tempfn, wxFile::write) )
    {
      file.Write( FreeSurferColorLUT_binary, FreeSurferColorLUT_binary_LEN );
      file.Flush();
      file.Close();
      ctd.filename = "FreeSurferColorLUT.txt";
      ctd.table = CTABreadASCII( tempfn.c_str() );
      ctd.name = "FreeSurferColorLUT";
      if ( ctd.table )
        m_tables.push_back( ctd );
    }
    else
    {
      std::cerr << "Can not load any stock look up tables." << std::endl;
    }
  }
}

LUTDataHolder::~LUTDataHolder()
{
  for ( int i = 0; i < GetCount(); i++ )
  {
    CTABfree( &m_tables[i].table );
  }
}

const char* LUTDataHolder::GetName( int i )
{
  if ( i < GetCount() )
    return m_tables[i].name.c_str();
  else
    return NULL;
}

COLOR_TABLE* LUTDataHolder::GetColorTable( int i )
{
  if ( i < GetCount() )
    return m_tables[i].table;
  else
    return NULL;
}

COLOR_TABLE* LUTDataHolder::GetColorTable( const char* name )
{
  for ( int i = 0; i < GetCount(); i++ )
  {
    if ( m_tables[i].name == name || m_tables[i].filename == name )
      return m_tables[i].table;
  }
  
  // if input is an index number instead of name, may still be valid
  wxString str = wxString::FromAscii(name);
  long n;
  if ( str.ToLong( &n ) && n < GetCount() )
    return m_tables[n].table;
  
  return NULL;
}

COLOR_TABLE* LUTDataHolder::LoadColorTable( const char* filename )
{
  // first trying to see if we've already loaded the lut from the same file
  wxFileName fn = wxFileName::FileName( wxString::FromAscii(filename) );
  fn.Normalize();
  COLOR_TABLE* ct = GetColorTable( fn.GetFullPath().ToAscii() );
  std::string filename_full = fn.GetFullPath().ToAscii();
  if ( !ct )
  {
    if ( !fn.FileExists() )
      fn = wxFileName::FileName( wxString( _("$FREESURFER_HOME/") ) + wxString::FromAscii(filename) );
    fn.Normalize();
    filename_full = fn.GetFullPath().char_str();
    
    for ( int i = 0; i < GetCount(); i++ )
    {
      if ( m_tables[i].filename == filename_full )
      {
        ct = m_tables[i].table;
        break;
      }
    }
  }
  
  // if so, reload the lut file and update the lut data. do not create new entry in the lut list
  if ( ct )
  {
    int nId = -1;
    for ( int i = 0; i < GetCount(); i++ )
    {
      if ( m_tables[i].table == ct )
      {
        nId = i;
        break;
      }
    }
    
    ct = CTABreadASCII( m_tables[nId].filename.c_str() );
    if ( ct )
    {
      CTABfree( &m_tables[nId].table );
      m_tables[nId].table = ct;
    }
    else
    {
      std::cerr << "Can not load color table from '" << filename << "'." << std::endl;
    }
  }
  // otherwise, load and create a new lut entry
  else
  {  
    ct = CTABreadASCII( filename_full.c_str() );
    if ( ct )
    {
      ColorTableData ctd;
      ctd.table = ct;
      ctd.filename = filename_full;
      wxFileName fn = wxFileName::FileName( wxString::FromAscii( filename_full.c_str() ) );
      ctd.name = fn.GetName().char_str();
      int n = 2;
      while ( GetColorTable( ctd.name.c_str() ) )
      {
        ctd.name = (fn.GetName() << _("_") << n).char_str();
        n++;
      }
      m_tables.push_back( ctd );
    }
    else
    {
      std::cerr << "Can not load color table from '" << filename << "'." << std::endl;
    }
  }
  
  return ct;
}

int LUTDataHolder::GetCount()
{
  return m_tables.size();
}

int LUTDataHolder::GetIndex( COLOR_TABLE* ct )
{
  for ( int i = 0; i < GetCount(); i++ )
  {
    if ( ct == m_tables[i].table )
      return i;
  }

  return -1;
}
