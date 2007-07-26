
#include <iostream>

#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "WrapperPolyData.h"

WrapperPolyData::WrapperPolyData()
{
  points = new PointContainerType;
  lines = new CellArrayType;
}

WrapperPolyData::~WrapperPolyData()
{
  delete points;
  delete lines;
}

double
WrapperPolyData::GetPointX(int idx) const
{ return (*points)[idx].x; }

double
WrapperPolyData::GetPointY(int idx) const
{ return (*points)[idx].y; }

double
WrapperPolyData::GetPointZ(int idx) const
{ return (*points)[idx].z; }

void 
WrapperPolyData::InsertPoint(int idx,
			  double x, double y, double z)
{
  PointType point;
  point.x = x;
  point.y = y;
  point.z = z;
  points->push_back(point);
}
void
WrapperPolyData::InsertNextCell(int numPoints)
{
  CellType cell;
  lines->push_back( cell );
}
void
WrapperPolyData::InsertCellPoint(int counter)
{
  lines->back().push_back(counter);
}

int
WrapperPolyData::GetNumberOfCells() const
{
  return lines->size();
}

void
WrapperPolyData::InitTraversal() const
{ citCells = lines->begin(); }

bool
WrapperPolyData::GoToNextCell() const
{ ++citCells; return (citCells!=lines->end()); }

int 
WrapperPolyData::GetCellNumberOfPoints() const
{ return citCells->size(); }

int 
WrapperPolyData::GetCellPoint(int idx) const
{ return (*citCells)[idx]; }

int
WrapperPolyData::GetNumberOfPoints() const
{ return points->size(); }

/////////////////////////////////

vtkPolyData*
GetPolyData(const WrapperPolyData* wrapper) 
{
  vtkPolyData* outData = vtkPolyData::New();
  vtkPoints* outPoints = vtkPoints::New();
  vtkCellArray* outLines = vtkCellArray::New();

  for(unsigned int ui(0), noPts(wrapper->GetNumberOfPoints());
      ui < noPts; ++ui)
    outPoints->InsertNextPoint(wrapper->GetPointX(ui),
			       wrapper->GetPointY(ui),
			       wrapper->GetPointZ(ui));

  if ( wrapper->GetNumberOfCells() )
    {
      wrapper->InitTraversal();
      do
	{
	  unsigned int noPts = wrapper->GetCellNumberOfPoints();
	  outLines->InsertNextCell( noPts );
	  for(unsigned int ui(0);
	      ui < noPts; ++ui )
	    outLines->InsertCellPoint( wrapper->GetCellPoint(ui) );
	}
      while ( wrapper->GoToNextCell() );
    }

  outData->SetPoints(outPoints);
  outData->SetLines(outLines);
  return outData;
}


//////////////////////

Line::Line()
  : data()
{}

Line::~Line()
{}

void
Line::AddPoint(double x, double y)
{
  data.push_back( std::make_pair(x,y) );
}

double
Line::GetPoint(int idx,
	       int dim) const
{
  if (dim==0)
    return data[idx].first;
  else if (dim==1)
    return data[idx].second;

  return .0;
}

int
Line::GetNumberOfPoints() const
{
  std::cout << " num pts = " << data.size() << std::endl;
  return data.size();
}

