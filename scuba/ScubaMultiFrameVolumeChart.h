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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:38 $
 *    $Revision: 1.5 $
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
