
#include <iostream>

#include "locator.h"

Locator::Locator()
{
  locator = 0;
}

Locator::~Locator()
{}

void
Locator::SetInputData(WrapperPolyData* input)
{
  vtkPolyData* data = GetPolyData(input);

  locator = vtkCellLocator::New();
  locator->SetDataSet(data);
  locator->SetAutomatic(1);
  locator->BuildLocator();
}

int
Locator::FindClosestPoint(double x,
			  double y,
			  double z)
{
  if (!locator) return -1;

  double input[3];
  input[0] = x;
  input[1] = y;
  input[2] = z;
  double returnedPoint[3];
  double dist2;
  vtkIdType id;
  int subId;
  
  locator->FindClosestPoint(input, returnedPoint,
			    id, subId, dist2);
  return id;
}

