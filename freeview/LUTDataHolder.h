/**
 * @file  LUTDataHolder.h
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
 *    $Revision: 1.9 $
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

#ifndef LUTDataHolder_h
#define LUTDataHolder_h

#include <QList>
#include <QString>

extern "C"
{
#include "colortab.h"
}

class LUTDataHolder
{
public:
  LUTDataHolder();
  virtual ~LUTDataHolder();

  QString GetName( int i );

  COLOR_TABLE* GetColorTable( int i );
  
  COLOR_TABLE* GetColorTable( const QString& name );

  int GetIndex( COLOR_TABLE* ct );

  int GetCount();
  
  COLOR_TABLE* LoadColorTable( const QString& fn );

protected:
  struct ColorTableData
  {
    COLOR_TABLE* table;
    QString  name;
    QString  filename;
  };

  QList<ColorTableData> m_tables;
};

#endif


