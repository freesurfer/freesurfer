/**
 * @file  ScubaVolumeROIIntensityChart.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/16 22:25:38 $
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


#ifndef ScubaVolumeROIIntensityChart_h
#define ScubaVolumeROIIntensityChart_h

#include "string_fixed.h"
#include "ChartWindow.h"
#include "Listener.h"
#include "TclCommandManager.h"
#include "VolumeCollection.h"
#include "ScubaROIVolume.h"

class ScubaVolumeROIIntensityChart : public Listener, public IDTracker<ScubaVolumeROIIntensityChart> {

public:

  // Order in which to sort the RAS points before sending them to the
  // chart.
  enum SortOrder { x, y, z };

  ScubaVolumeROIIntensityChart( VolumeCollection& iVolume,
                                ScubaROIVolume& iROI, SortOrder iOrder );
  ~ScubaVolumeROIIntensityChart();

  virtual void
  DoListenToMessage ( std::string isMessage, void* iData );

protected:

  // Builds the list of points and sends the data to the chart. The
  // first function switches on the sort type and calls the templated
  // function.
  void Draw();
  template <typename T> void DrawWithSort();


  VolumeCollection& mVolume;
  ScubaROIVolume&   mROI;
  ChartWindow*      mChart;
  SortOrder         mSortOrder;

  // These compare the first VolumeLocation to the second and return
  // true if the given coordinate in the first is greater than the
  // second. If they're equal, it will move onto the next coordinate.
  class VolumeLocationRASXComparatorGT {
  public:
    bool operator()(VolumeLocation const& v1, VolumeLocation const& v2) const;
  };
  class VolumeLocationRASYComparatorGT {
  public:
    bool operator()(VolumeLocation const& v1, VolumeLocation const& v2) const;
  };
  class VolumeLocationRASZComparatorGT {
  public:
    bool operator()(VolumeLocation const& v1, VolumeLocation const& v2) const;
  };
};

class ScubaVolumeROIIntensityChartFactory : public DebugReporter, public TclCommandListener {

public:

  static ScubaVolumeROIIntensityChartFactory& GetFactory();

  virtual TclCommandResult
  DoListenToTclCommand ( char* isCommand, int iArgc, char** iArgv );

protected:
  static bool mbAddedTclCommands;

};


#endif
