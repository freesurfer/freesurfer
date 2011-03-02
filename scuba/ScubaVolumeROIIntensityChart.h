/**
 * @file  ScubaVolumeROIIntensityChart.h
 * @brief Given a volume and an ROI, uses a ChartWindow to plot intenity data
 *
 * Uses a ChartWindow to plot the voxel index on the x axes and the
 * intensity on the y axis. Voxel index can be sorted by x, y, or z.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
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
