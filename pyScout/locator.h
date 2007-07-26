
#ifndef _h_locator_h
#define _h_locator_h

#include <vtkCellLocator.h>

#include "WrapperPolyData.h"

class Locator
{
 public:
  Locator();
  ~Locator();
  void SetInputData(WrapperPolyData* input);

  int FindClosestPoint(double x, double y, double z);
 private:
  vtkCellLocator* locator;
};

#endif
