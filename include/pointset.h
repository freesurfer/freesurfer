#ifndef POINTSET_H
#define POINTSET_H

#include <string>
#include <vector>


/// \class PointSet
/// \brief A 3d point collection
///
/// PointSet is a collection of 3d points used mostly by the user
/// within freeview. Freeview has it's own PointSet class, but we don't
/// want to rely on Qt for such a simple class. This class is a
/// work-in-progress - and right now, pointsets can only be
/// saved as json files

class fsPointSet
{
public:

  // simple 3d point
  struct Point {
    double x, y, z, value;
    Point() : x(0), y(0), z(0), value(1) {}
    Point(double _x, double _y, double _z, double _value = 1.0)
      : x(_x), y(_y), z(_z), value(_value) {}
  };

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

  // let's make PointSet an extension of std::vector for ease of use
  void clear() { points.clear(); }
  void remove(int index) { points.erase(points.begin() + index); }
  Point& operator[](int index) { return points[index]; }

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
