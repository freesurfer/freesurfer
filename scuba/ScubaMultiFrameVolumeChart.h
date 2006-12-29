/**
 * @file  ScubaMultiFrameVolumeChart.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2002-2007,
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


#ifndef ScubaMultiFrameVolumeChart_h
#define ScubaMultiFrameVolumeChart_h

#include "string_fixed.h"
#include <map>
#include "ChartWindow.h"
#include "Listener.h"
#include "TclCommandManager.h"
#include "VolumeCollection.h"
#include "ScubaROIVolume.h"

// This allows you to create charts based on volumes that will get the
// current cursor and draw the volume's values at that cursor over
// multiple frames. If there is time metadata available, it will use
// that to format the graph into conditions and times.

class ScubaMultiFrameVolumeChart : public Listener, public IDTracker<ScubaMultiFrameVolumeChart> {

public:

  ScubaMultiFrameVolumeChart( VolumeCollection& iVolume );
  ScubaMultiFrameVolumeChart( VolumeCollection& iVolume,
                              ScubaROIVolume& iROI );
  ~ScubaMultiFrameVolumeChart();

  virtual void
  DoListenToMessage ( std::string isMessage, void* iData );

protected:

  void Draw();

private:

  VolumeCollection& mVolume;
  ScubaROIVolume*   mROI;
  ChartWindow*      mChart;
};

class ScubaMultiFrameVolumeChartFactory : public TclCommandManager {

public:

  static ScubaMultiFrameVolumeChartFactory& GetFactory();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

protected:
  static bool mbAddedTclCommands;

};


#endif
