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

#ifndef FSWayPoints_h
#define FSWayPoints_h


#include "vtkMatrix4x4.h"
#include "CommonDataStruct.h"


extern "C"
{
#include "label.h"
}

#include <string>
#include <vector>

class FSVolume;
class wxWindow;
class wxCommandEvent;

struct WayPoint
{
  double pt[3];
  double value;
};

class FSWayPoints
{
public:
  FSWayPoints();
  virtual ~FSWayPoints();

  bool ReadAsLabel( const char* filename );
  bool ReadAsControlPoints( const char* filename );
  bool WriteAsLabel( const char* filename );
  bool WriteAsControlPoints( const char* filename );
  
  static bool IsLabelFormat( const char* filename );

  void UpdateLabel( std::vector<WayPoint>& points_in, FSVolume* vol_ref );
  void LabelToWayPoints( std::vector<WayPoint>& points_out, FSVolume* vol_ref );

protected:
  // use label to save way points
  LABEL*   m_label;

};

#endif


