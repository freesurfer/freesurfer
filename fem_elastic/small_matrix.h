
#ifndef H_SMALL_MATRIX_H
#define H_SMALL_MATRIX_H

#include <iostream>

class SmallMatrix
{
public:
  SmallMatrix(); // generates an empty matrix
  SmallMatrix(const SmallMatrix& sm);
  SmallMatrix& operator=(const SmallMatrix& sm);
  SmallMatrix(int rows, int cols);

  ~SmallMatrix();

  void set(double val);

  void identity(int size);

  double  operator()(int i, int j) const;
  double& operator()(int i, int j);

  int rows() const
  {
    return m_rows;
  }
  int cols() const
  {
    return m_cols;
  }

  void  set_block(const SmallMatrix& m,
                  size_t start_row,
                  size_t start_col);

  SmallMatrix operator*(const SmallMatrix& m);
  const SmallMatrix& operator*=(const double val);

  SmallMatrix operator+(const SmallMatrix& m) const;
  SmallMatrix operator-(const SmallMatrix& m) const;

  void operator+=(const SmallMatrix& m);
  void operator-=(const SmallMatrix& m);

  SmallMatrix transpose() const; // will generate a new matrix
  void inplace_transpose(); // will perform the operation in-place

  double norm() const; // returns the 2-norm of a matrix (as a vector)

private:
  double* m_pdata;
  int     m_rows, m_cols;

  void clone(const SmallMatrix& sm);
};

std::ostream& operator<<(std::ostream& os,
                         const SmallMatrix& sm);

#endif // H_SMALL_MATRIX_H
