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
  class VolumeLocationRASXComparatorGT:
    public std::binary_function<VolumeLocation*,VolumeLocation*,bool> {
    public:
    bool operator()(const VolumeLocation* v1, const VolumeLocation* v2) const;
  };
  class VolumeLocationRASYComparatorGT:
    public std::binary_function<VolumeLocation*,VolumeLocation*,bool> {
    public:
    bool operator()(const VolumeLocation* v1, const VolumeLocation* v2) const;
  };
  class VolumeLocationRASZComparatorGT:
    public std::binary_function<VolumeLocation*,VolumeLocation*,bool> {
    public:
    bool operator()(const VolumeLocation* v1, const VolumeLocation* v2) const;
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
