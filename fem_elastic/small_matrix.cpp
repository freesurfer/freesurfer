
#include <algorithm>
#include <math.h>

#include "small_matrix.h"

SmallMatrix::SmallMatrix()
{
  m_pdata = NULL;
  m_rows = -1;
  m_cols = -1;
}

SmallMatrix::SmallMatrix(int _rows,
                         int _cols)
{
  m_rows = _rows;
  m_cols = _cols;

  m_pdata = new double[m_rows*m_cols];
}

SmallMatrix::SmallMatrix(const SmallMatrix& sm)
{
  m_pdata = NULL;
  m_rows = -1;
  m_cols = -1;

  clone(sm);
}

void
SmallMatrix::set(double val)
{
  std::fill_n(m_pdata, m_rows*m_cols, val);
}

SmallMatrix&
SmallMatrix::operator=(const SmallMatrix& sm)
{
  clone(sm);
  return *this;
}

SmallMatrix::~SmallMatrix()
{
  if ( m_pdata ) delete[] m_pdata;
}

void
SmallMatrix::identity(int size)
{
  if ( m_pdata ) delete[] m_pdata;
  m_rows = size;
  m_cols = size;

  m_pdata = new double[size*size];
  std::fill_n(m_pdata, size*size,
              0.0);

  for (int i=0; i<size; ++i)
    m_pdata[ i+ i*size ] = 1.0;
}

double
SmallMatrix::operator()(int i,
                        int j) const
{
  return m_pdata[ i*m_cols + j ];
}

double&
SmallMatrix::operator()(int i,
                        int j)
{
  return m_pdata[ i*m_cols + j ];
}


void
SmallMatrix::clone(const SmallMatrix& sm)
{
  if ( m_pdata ) delete[] m_pdata;

  m_rows = sm.rows();
  m_cols = sm.cols();

  m_pdata = new double[rows()*cols()];
  std::copy(sm.m_pdata, sm.m_pdata+rows()*cols(), m_pdata);
}

void
SmallMatrix::set_block(const SmallMatrix& m,
                       size_t start_row,
                       size_t start_col)
{
  size_t max_row, max_col;
  max_row = std::min( rows(), int(start_row+m.rows()) );
  max_col = std::min( cols(), int(start_col+m.cols()) );
  for ( size_t i=start_row;
        i < max_row; ++i)
    for (size_t j=start_col;
         j< max_col; ++j)
      (*this)(i,j) = m(i-start_row, j-start_col);
}

SmallMatrix
SmallMatrix::operator*(const SmallMatrix& m)
{
  if ( this->cols() != m.rows() )
  {
    std::cerr << " SmallMatrix::operator* -> wrong matrix sizes\n"
    << " this size = " << this->rows() << " , " << this->cols() << std::endl
    << " m size = " << m.rows() << " , " << m.cols() << std::endl;
    return SmallMatrix();
  }
  SmallMatrix retVal( this->rows(), m.cols() );

  double dval = 0.0;
  for (int i=0; i<retVal.rows(); ++i)
    for (int j=0; j<retVal.cols(); ++j)
    {
      dval = 0;
      for (int k=0; k<this->cols(); ++k)
        dval += (*this)(i,k) * m(k,j);
      retVal(i,j) = dval;
    } // next j, i
  return retVal;
}

const SmallMatrix&
SmallMatrix::operator*=(const double val)
{
  double *pos = m_pdata;
  int no = m_rows * m_cols;

  for (int i=0; i< no; ++i, ++pos)
    *pos *= val;

  return *this;
}

SmallMatrix
SmallMatrix::transpose() const
{
  SmallMatrix retVal(this->cols(), this->rows());

  for (int i=0; i<this->rows(); ++i)
    for (int j=0; j<this->cols(); ++j)
      retVal(j,i) = (*this)(i,j);

  return retVal;
}

void
SmallMatrix::inplace_transpose()
{
  double *pdata = new double[m_rows*m_cols];

  for (int i=0; i<m_rows; ++i)
    for (int j=0; j<m_cols; ++j)
      pdata[ j*m_rows + i ] = (*this)(i,j);

  delete m_pdata;
  m_pdata = pdata;

  int ibuf = m_rows;

  m_rows = m_cols;
  m_cols = ibuf;
}

SmallMatrix
SmallMatrix::operator+(const SmallMatrix& m) const
{
  // error if the two matrices are not the same size
  if ( this->rows() != m.rows() ||
       this->cols() != m.cols() )
  {
    std::cerr << " SmallMatrix::operator+ -> matrices have different sizes\n";
    return SmallMatrix();
  }

  SmallMatrix retVal(this->rows(), this->cols() );

  double* pthis = m_pdata;
  double* pm    = m.m_pdata;
  double* pret  = retVal.m_pdata;

  int no = this->rows() * this->cols();

  for (int i=0; i<no; ++i, ++pthis, ++pm, ++pret)
  {
    *pret = *pm + *pthis;
  } // next i, pthis, pm, pret

  return retVal;
}

SmallMatrix
SmallMatrix::operator-(const SmallMatrix& m) const
{
  // error if the two matrices do not have the same size
  if ( this->rows() != m.rows() ||
       this->cols() != m.cols() )
  {
    std::cerr << " SmallMatrix::operator- -> matrices have different sizes\n";
    return SmallMatrix();
  }

  SmallMatrix retVal(this->rows(), this->cols());

  double *pthis = m_pdata;
  double *pm    = m.m_pdata;
  double *pret  = retVal.m_pdata;

  int no = this->rows() * this->cols();

  for (int i=0; i<no; ++i, ++pthis, ++pm, ++pret)
  {
    *pret = *pthis - *pm;
  } // next i, pthis, pm, pret

  return retVal;
}

void
SmallMatrix::operator+=(const SmallMatrix& m)
{
  // error if the two matrices do not have the same size
  if ( this->rows() != m.rows() ||
       this->cols() != m.cols() )
  {
    std::cerr << " SmallMatrix::operator+= -> matrices have different sizes\n";
    exit(1);
  }

  double *pthis = m_pdata;
  double *pm    = m.m_pdata;

  int no = this->rows() * this->cols();

  for (int i=0; i<no; ++i, ++pthis, ++pm)
  {
    *pthis += *pm;
  } // next i, pthis, pm
}

void
SmallMatrix::operator-=(const SmallMatrix& m)
{
  // error if the two matrices do not have the same size
  if ( this->rows() != m.rows() ||
       this->cols() != m.cols() )
  {
    std::cerr << " SmallMatrix::operator-= -> matrices have different sizes\n";
    return;
  }

  double *pthis = m_pdata;
  double *pm    = m.m_pdata;

  int no = this->rows() * this->cols();

  for (int i=0; i<no; ++i, ++pthis, ++pm)
  {
    *pthis -= *pm;
  } // next i, pthis, pm
}

double
SmallMatrix::norm() const
{
  double dret(0);

  double *pthis = m_pdata;
  int no = this->rows() * this->cols();

  for (int i=0; i<no; ++i, ++pthis)
  {
    dret += *pthis * *pthis;
  }

  return std::sqrt(dret);
}


//-------------------------

std::ostream&
operator<<(std::ostream& os,
           const SmallMatrix& sm)
{
  for (int i=0; i<sm.rows(); ++i)
  {
    for (int j=0; j<sm.cols(); ++j)
      os << sm(i,j) << " ";
    os << std::endl;
  }
  return os;
}
