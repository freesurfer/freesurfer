#ifndef ScubaMultiFrameVolumeChart_h
#define ScubaMultiFrameVolumeChart_h

#include "string_fixed.h"
#include <map>
#include "ChartWindow.h"
#include "Listener.h"
#include "TclCommandManager.h"
#include "VolumeCollection.h"

// This allows you to create charts based on volumes that will get the
// current cursor and draw the volume's values at that cursor over
// multiple frames. If there is time metadata available, it will use
// that to format the graph into conditions and times.

class ScubaMultiFrameVolumeChart : public Listener, public IDTracker<ScubaMultiFrameVolumeChart> {
  
 public:

  ScubaMultiFrameVolumeChart( VolumeCollection& iVolume );
  ~ScubaMultiFrameVolumeChart();

  virtual void
    DoListenToMessage ( std::string isMessage, void* iData );

 protected:

  void Draw();

 private:

  VolumeCollection& mVolume;
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
