#ifndef POINTSET_H
#define POINTSET_H

#include <string>
#include <vector>


/// \class fsPointSet
/// \brief A 3d point collection
///
/// fsPointSet is a collection of 3d points used mostly by the user
/// within freeview. Freeview has it's own fsPointSet class, but we don't
/// want to rely on Qt for such a simple class. This class is a
/// work-in-progress - and right now, pointsets can only be
/// saved as json files

class fsPointSet
{
public:

  // simple 3d point
  struct Point {
    double x, y, z, value;
    int index=0, count=0;
    std::string labelname;
    Point() : x(0), y(0), z(0), value(1) {}
    Point(double _x, double _y, double _z, double _value = 1.0)
      : x(_x), y(_y), z(_z), value(_value) {}
  };
  int IncludeZeroCount = 1;
  int IncludeZeroIndex = 1;

  bool save(std::string filename);
  // bool load(std::string filename);
  bool save_as_ctrlpoint(std::string filename);

  int add(Point point) {
    points.push_back(point);
    return points.size() - 1;
  }
  int add(double x, double y, double z, double value) {
    points.push_back({x, y, z, value});
    return points.size() - 1;
  }
  int set(int i, double x, double y, double z, double value) {
    if(i < 0 || i >= points.size()) return(1);
    points[i].x = x;
    points[i].y = y;
    points[i].z = z;
    points[i].value = value;
    return(0);
  }
  int set(int i, Point p){
    if(i < 0 || i >= points.size()) return(1);
    points[i].x = p.x;
    points[i].y = p.y;
    points[i].z = p.z;
    points[i].value = p.value;
    points[i].count = p.count;
    points[i].index = p.index;
    return(0);
  }
  Point get(int i){ 
    Point p;
    if(i < 0 || i >= points.size()) {
      p.count = -1;
      return(p);
    }
    return(points[i]);
  }

  int print(FILE *fp){
    int i=0;
    for (const Point& point : points) {
      if(point.count == 0 && !IncludeZeroCount) continue;
      if(point.index == 0 && !IncludeZeroIndex) continue;
      fprintf(fp,"%3d %4d %6d %8.4lf %8.4lf %8.4lf %8.4lf\n",
	      i,point.index,point.count,point.x,point.y,point.z,point.value);
    }
    fflush(fp);
    return(0);
  }
  int printCoords(FILE *fp){
    for (const Point& point : points) {
      if(point.count == 0 && !IncludeZeroCount) continue;
      if(point.index == 0 && !IncludeZeroIndex) continue;
      fprintf(fp,"%4d %8.4lf %8.4lf %8.4lf\n",point.index,point.x,point.y,point.z);
    }
    fflush(fp);
    return(0);
  }
  int writeCoords(std::string fname){
    FILE *fp = fopen(fname.c_str(),"w");
    if(fp==NULL) return(1);
    printCoords(fp);
    fclose(fp);
    return(0);
  }

  // let's make fsPointSet an extension of std::vector for ease of use
  void clear() { points.clear(); }
  void remove(int index) { points.erase(points.begin() + index); }
  Point& operator[](int index) { return points[index]; }

  int writeCentroidTable(std::string outfile,std::string ltafile);

  // let's support iteration as well
  typedef std::vector<Point>::iterator iterator;
  typedef std::vector<Point>::const_iterator const_iterator;
  iterator begin() { return points.begin(); }
  iterator end() { return points.end(); }

  std::string vox2ras = "scanner_ras";

 private:
  std::vector<Point> points;

};

fsPointSet loadfsPointSet(std::string filename);

#endif
