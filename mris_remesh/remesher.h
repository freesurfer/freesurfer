#pragma once

#include <array>
#include <vector>
#include <list>
#include <iterator>

#include "mrisurf.h"


// -------- 3D Vector --------


class Vector
{
public:

  Vector() {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
  }

  Vector(const double x, const double y, const double z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }

  Vector(const Vector& vect) {
    v[0] = vect.v[0];
    v[1] = vect.v[1];
    v[2] = vect.v[2];
  }

  inline double& operator[](int);
  inline const double& operator[](int) const;

  friend inline std::ostream& operator<<(std::ostream&, const Vector&);

  inline Vector& operator=(const Vector&);
  inline bool operator==(const Vector&) const;
  inline bool operator!=(const Vector&) const;

  inline Vector& operator+=(const Vector&);
  inline Vector& operator-=(const Vector&);
  inline Vector& operator*=(const double);
  inline Vector operator+(const Vector&) const;
  inline Vector operator-(const Vector&) const;
  inline double operator*(const Vector&) const;
  inline Vector operator*(const double) const;

  inline double norm(void) const;
  inline double normSquared(void) const;

private:
  double v[3];
};

inline std::ostream& operator<<(std::ostream &os, const Vector& vect) {
  os << "[" << vect.v[0] << ", " << vect.v[1] << ", " << vect.v[2] << "]";
  return os;
}

inline double& Vector::operator[](int index) {
  return v[index];
}

inline const double& Vector::operator[](int index) const {
  return v[index];
}

inline Vector& Vector::operator=(const Vector &vect) {
  v[0] = vect[0];
  v[1] = vect[1];
  v[2] = vect[2];
  return *this;
}

inline bool Vector::operator==(const Vector &vect) const {
  return (v[0] == vect[0] && v[1] == vect[1] && v[2] == vect[2]) ? true : false;
}

inline bool Vector::operator!=(const Vector &vect) const {
  return (v[0] == vect[0] && v[1] == vect[1] && v[2] == vect[2]) ? false : true;
}

inline Vector& Vector::operator+=(const Vector &vect) {
  v[0] += vect[0];
  v[1] += vect[1];
  v[2] += vect[2];
  return *this;
}

inline Vector& Vector::operator-=(const Vector &vect) {
  v[0] -= vect[0];
  v[1] -= vect[1];
  v[2] -= vect[2];
  return *this;
}

inline Vector& Vector::operator*=(const double scalar) {
  v[0] *= scalar;
  v[1] *= scalar;
  v[2] *= scalar;
  return *this;
}

inline Vector Vector::operator+(const Vector &vect) const {
  return Vector(v[0] + vect[0], v[1] + vect[1], v[2] + vect[2]);
}

inline Vector Vector::operator-(const Vector &vect) const {
  return Vector(v[0] - vect[0], v[1] - vect[1], v[2] - vect[2]);
}

inline double Vector::operator*(const Vector &vect) const {
  return v[0] * vect[0] + v[1] * vect[1] + v[2] * vect[2];
}

inline Vector Vector::operator*(const double scalar) const {
  return Vector(v[0] * scalar, v[1] * scalar, v[2] * scalar);
}

inline double Vector::norm() const {
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

inline double Vector::normSquared() const {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}       

inline Vector operator*(const double scalar, const Vector vect) {
  return vect * scalar;
}

inline Vector cross(const Vector u, const Vector v) {
   return Vector(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]);
}

// support arbitrary vector printing
template <typename T> std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if (!v.empty()) {
    out << '[';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "]";
  }
  return out;
}


// -------- 3D Matrix --------


class Matrix
{
public:
  Matrix() {}
  Matrix(
    const double m11, const double m12, const double m13,
    const double m21, const double m22, const double m23,
    const double m31, const double m32, const double m33) {
    m[0] = m11; m[1] = m12; m[2] = m13;
    m[3] = m21; m[4] = m22; m[5] = m23;
    m[6] = m31; m[7] = m32; m[8] = m33;
  }

  inline Matrix operator-(const Matrix&) const;
  inline Matrix operator*(const double) const;   
  inline Matrix operator*(const Matrix&) const;
  inline Vector operator*(const Vector&) const;

  inline double& operator()(int, int);

private:
  double m[9];
};


inline double& Matrix::operator()(int i, int j) {
  return m[3 * i + j];
}


inline Matrix Matrix::operator-(const Matrix &mat) const {
  return Matrix(
    m[0] - mat.m[0], m[1] - mat.m[1], m[2] - mat.m[2],
    m[3] - mat.m[3], m[4] - mat.m[4], m[5] - mat.m[5],
    m[6] - mat.m[6], m[7] - mat.m[7], m[8] - mat.m[8]);
}

inline Matrix Matrix::operator*(const Matrix &mat) const {
  return Matrix(
    m[0] * mat.m[0] + m[1] * mat.m[3] + m[2] * mat.m[6], 
    m[0] * mat.m[1] + m[1] * mat.m[4] + m[2] * mat.m[7],
    m[0] * mat.m[2] + m[1] * mat.m[5] + m[2] * mat.m[8], 
    m[3] * mat.m[0] + m[4] * mat.m[3] + m[5] * mat.m[6], 
    m[3] * mat.m[1] + m[4] * mat.m[4] + m[5] * mat.m[7], 
    m[3] * mat.m[2] + m[4] * mat.m[5] + m[5] * mat.m[8], 
    m[6] * mat.m[0] + m[7] * mat.m[3] + m[8] * mat.m[6], 
    m[6] * mat.m[1] + m[7] * mat.m[4] + m[8] * mat.m[7], 
    m[6] * mat.m[2] + m[7] * mat.m[5] + m[8] * mat.m[8]);  
}

inline Matrix Matrix::operator*(const double scalar) const {
  return Matrix(
    m[0] * scalar, m[1] * scalar, m[2] * scalar,
    m[3] * scalar, m[4] * scalar, m[5] * scalar, 
    m[6] * scalar, m[7] * scalar, m[8] * scalar);
}

inline Vector Matrix::operator*(const Vector &v) const {
  return Vector(
    m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
    m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
    m[6] * v[0] + m[7] * v[1] + m[8] * v[2]);
}


// -------- Sparse Matrix --------

class SparseMatrix
{
public:

  struct entry
  {
    int pos;
    double value;
  };

  void redim(int m, int n);

  void fillzero();
  bool erase(int i, int j);
  double eraseGetVal(int i, int j);

  bool isSet(int i, int j) const;
  double getVal(int i, int j) const;
  void setVal(int i, int j, double d);

  int getNoRow() { return row; }
  const std::list<entry>& getRow(unsigned int r) const { return entries[r]; };
  void setRow(unsigned int r,const std::list < entry>& rr) { entries[r] = rr; };

  const double& getRef(int, int) const;
  double& getRef(int, int);

  inline double& operator()(int, int);
  inline const double& operator()(int, int) const;

  double getMax() const;
  double getMin() const;

private:
  int nonzero = 0;
  int row = 0;
  int col = 0;
  std::vector<std::list<entry>> entries;
};

inline double& SparseMatrix::operator()(int i, int j) { return getRef(i, j); }
inline const double& SparseMatrix::operator()(int i, int j) const { return getRef(i, i); }


// -------- Remeshing Class --------


class Remesher
{
public:

  Remesher(const MRIS *surf);
  MRIS * toSurface();

  void remeshBK(unsigned int it, double l = -1, bool ridge = false);
  
  int verbose = 0;

private:
  std::vector<Vector> points3d;
  std::vector<std::vector<int>> tria;

  void init();

  void createVtoT();
  void createEdges();
  void createVBoundary();
  void createQualities(double tmin = 0.04) const;

  bool contractEdge(int eidx, bool cleanup=true);
  bool rmTrias(const std::vector< int > & tnums, bool fix = true);
  bool rmTrias(bool fix = true);
  bool rmTriaInVtoT(int vidx, int tidx);
  bool removeTria(unsigned int tidx);
  bool removeHangingEdge(unsigned int eidx);

  int rmFreeVertices(bool fix);

  int insertVertexOnEdge(const Vector& pos);
  void createElementN(bool force = false) const;
  void createVertexNSimple(bool force = false) const;

  std::vector<int> get1Ring(unsigned int i, bool ordered = true) const;

  SparseMatrix v_to_e;
  std::vector<std::vector<int>> t_to_e;
  std::vector<std::vector<int>> e_to_t;
  std::vector<std::vector<int>> e_to_v;
  std::vector<std::vector<int>> v_to_t;
  std::vector<std::vector<short>> v_to_lv;

  int bedgecount;
  bool ismanifold;

  std::vector<bool> onridge;
  std::vector<bool> onboundary;
  std::vector<int> indexwoboundary;

  int innercount;
  int pointcount;

  mutable std::vector<Vector> vertexn;
  mutable std::vector<Vector> elementn;

  mutable std::vector<double> qualities;
  mutable double worstqual;
  mutable double averagequal;

  bool mergeEdges(unsigned int etarg, unsigned int edel);
  bool rmEdgeInVtoE(int eidx);
  double getAverageEdgeLength() const;
  int replaceVertexInEdges(int vold, int vnew, bool simulate = false, int vtip0=-1, int vtip1=-1);
  int contractShortEdges(double l);
  void contractEdgeInTria(int eidx, int tidx);
  int createVertexOnEdge(const Vector& pos);
  int insertVerticesOnLongEdges(double l);
  int makeRidge(double angle =-1);

  int removeCollapsedTrias(int vidx, int mididx);

  double computeArea() const;
  double computeArea(unsigned int tidx) const;

  double computeQuality(int tInd) const;
  int checkStructure();

  void tangentialSmoothing(int it = 1,bool gravity = true, bool tangent = true);
};
