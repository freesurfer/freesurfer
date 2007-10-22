/**
 * @file  ScubaMultiFrameVolumeChart.h
 * @brief Creates charts for values in multiframe volumes
 *
 * Uses the ChartWindow class to create charts based on volumes that
 * will get the current cursor and draw the volume's values at that
 * cursor over multiple frames. If there is time metadata available,
 * it will use that to format the graph into conditions and times.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/22 04:39:29 $
 *    $Revision: 1.4 $
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
