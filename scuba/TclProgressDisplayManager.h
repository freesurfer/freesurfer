#ifndef TclProgressDisplayManager_h
#define TclProgressDisplayManager_h

#include "ProgressDisplayManager.h"
#include "string_fixed.h"
#include <list>

class TclProgressDisplayManager : public ProgressDisplayManager{

 public:

  TclProgressDisplayManager() {}
  ~TclProgressDisplayManager() {}

  void NewTask ( std::string isTitle,
		 std::string isText,
		 bool ibUseMeter,
		 std::list<std::string> ilsButtons );
  
  void UpdateTask ( std::string isText,
		    float iPercent );
  
  int CheckTaskForButton ();
  
  void EndTask ();

 protected:

  std::list<std::string> mlButtons;

};


#endif
