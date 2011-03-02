/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:02 $
 *    $Revision: 1.6 $
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

#ifndef LUTDataHolder_h
#define LUTDataHolder_h

#include <string>
#include <vector>

extern "C"
{
#include "colortab.h"
}

class LUTDataHolder
{
public:
  LUTDataHolder();
  virtual ~LUTDataHolder();

  const char* GetName( int i );

  COLOR_TABLE* GetColorTable( int i );
  
  COLOR_TABLE* GetColorTable( const char* name );

  int GetIndex( COLOR_TABLE* ct );

  int GetCount();
  
  COLOR_TABLE* LoadColorTable( const char* fn );

protected:
  struct ColorTableData
  {
    COLOR_TABLE* table;
    std::string  name;
    std::string  filename;
  };

  std::vector<ColorTableData> m_tables;
};

#endif


