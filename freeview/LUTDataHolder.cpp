/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2009/04/29 22:53:51 $
 *    $Revision: 1.2.2.3 $
 *
 * Copyright (C) 2008-2009,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <wx/wx.h>
#include <wx/filename.h>
#include "LUTDataHolder.h"

LUTDataHolder::LUTDataHolder()
{
  ColorTableData ctd;
  wxFileName fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/FreeSurferColorLUT.txt") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( (char*)ctd.filename.c_str() );
  ctd.name = "FreeSurferColorLUT";
  ctd.filename = fn.GetFullPath().char_str();
  if ( ctd.table )
    m_tables.push_back( ctd );

  fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/tkmeditParcColorsCMA") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( (char*)ctd.filename.c_str() );
  ctd.name = "tkmeditParcColorsCMA";
  ctd.filename = fn.GetFullPath().char_str();
  if ( ctd.table )
    m_tables.push_back( ctd );

  fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/surface_labels.txt") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( (char*)ctd.filename.c_str() );
  ctd.name = "surface_labels";
  ctd.filename = fn.GetFullPath().char_str();
  if ( ctd.table )
    m_tables.push_back( ctd );
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
  COLOR_TABLE* ct = GetColorTable( filename );
  std::string filename_full = filename;
  if ( !ct )
  {
    wxFileName fn = wxFileName::FileName( wxString::FromAscii(filename) );
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
    
    ct = CTABreadASCII( (char*)m_tables[nId].filename.c_str() );
    if ( ct )
    {
      CTABfree( &m_tables[nId].table );
      m_tables[nId].table = ct;
    }
    else
    {
      std::cerr << "Cannot load color table from '" 
                << filename << "'." << std::endl;
    }
  }
  // otherwise, load and create a new lut entry
  else
  {  
    ct = CTABreadASCII( (char*)filename_full.c_str() );
    if ( ct )
    {
      ColorTableData ctd;
      ctd.table = ct;
      ctd.filename = filename_full;
      wxFileName fn = wxFileName::FileName
        ( wxString::FromAscii( filename_full.c_str() ) );
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
      std::cerr << "Cannot load color table from '" 
                << filename << "'." << std::endl;
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
