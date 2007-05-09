#ifndef __FieldLineInitPathStrategy_h
#define __FieldLineInitPathStrategy_h

#include "InitializePath.h"

class FieldLineInitPathStrategy : public InitializePath
{
public:

  FieldLineInitPathStrategy();
  FieldLineInitPathStrategy( PoistatsModel* model );
  ~FieldLineInitPathStrategy();  

  /**
   * This does the work of finding the initial path.
   */
  void CalculateInitialPath();

};

#endif 

