/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2009/03/31 22:00:12 $
 *    $Revision: 1.6 $
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
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
  ctd.name = "FreeSurferColorLUT";
  ctd.filename = fn.GetFullPath().char_str();
  if ( ctd.table )
    m_tables.push_back( ctd );

  fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/tkmeditParcColorsCMA") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
  ctd.name = "tkmeditParcColorsCMA";
  ctd.filename = fn.GetFullPath().char_str();
  if ( ctd.table )
    m_tables.push_back( ctd );

  fn = wxFileName::FileName
    ( _("$FREESURFER_HOME/surface_labels.txt") );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
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
  wxString str = name;
  long n;
  if ( str.ToLong( &n ) && n < GetCount() )
    return m_tables[n].table;
  
  return NULL;
}


COLOR_TABLE* LUTDataHolder::LoadColorTable( const char* filename )
{
  if ( GetColorTable( filename ) )
    return GetColorTable( filename );
  
  ColorTableData ctd;
  wxFileName fn = wxFileName::FileName( filename );
  if ( !fn.FileExists() )
    fn = wxFileName::FileName( wxString("$FREESURFER_HOME/") + filename );
  fn.Normalize();
  ctd.filename = fn.GetFullPath().char_str();
  // if the same file already loaded, return the table
  for ( int i = 0; i < GetCount(); i++ )
  {
    if ( m_tables[i].filename == ctd.filename )
      return m_tables[i].table;
  }
  
  ctd.table = CTABreadASCII( ctd.filename.c_str() );
  if ( ctd.table )
  {
    ctd.name = fn.GetName().c_str();
    int n = 2;
    while ( GetColorTable( ctd.name.c_str() ) )
    {
      ctd.name = (fn.GetName() << "_" << n).c_str();
      n++;
    }
    m_tables.push_back( ctd );
  }
  else
  {
    std::cerr << "Can not load color table from '" << ctd.filename.c_str() << "'." << std::endl;
  }
  
  return ctd.table;
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
