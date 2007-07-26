
#ifndef _wrapper_polyData_h
#define _wrapper_polyData_h

#include <algorithm>
#include <vector>

#include <vtkPolyData.h>


/*

for passing data from Python to C

*/
class WrapperPolyData
{
 public:
  WrapperPolyData();
  ~WrapperPolyData();
  
  void InsertPoint(int idx, double, double, double);
  void InsertNextCell(int);
  void InsertCellPoint(int);

  int GetNumberOfCells() const;
  void InitTraversal() const;
  bool GoToNextCell() const;
  int  GetCellNumberOfPoints() const;
  int  GetCellPoint(int idx) const;
  double GetPointX(int idx) const;
  double GetPointY(int idx) const;
  double GetPointZ(int idx) const;
  int GetNumberOfPoints() const;
  
 private:
  typedef struct { double x,y,z; } PointType;
  typedef std::vector<PointType>  PointContainerType;
  typedef std::vector<int> CellType;
  typedef std::vector<CellType> CellArrayType;

  PointContainerType *points;
  CellArrayType *lines;

  mutable CellArrayType::const_iterator citCells;
};

vtkPolyData* GetPolyData(const WrapperPolyData*);


/*

for passing objects from C to Python
mainly lists

*/
class Line
{
 public:
  typedef std::pair<double,double> Point;

  Line();
  ~Line();
  void AddPoint(double x, double y);
  int GetNumberOfPoints() const;
  double GetPoint(int idx, int dim) const; // simplify interface
 private:
  typedef std::vector<Point> Container;
  Container data;
};



#endif
