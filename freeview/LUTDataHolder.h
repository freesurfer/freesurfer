/**
 * @brief LUT data object.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#ifndef LUTDataHolder_h
#define LUTDataHolder_h

#include <QList>
#include <QString>



#include "colortab.h"


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


