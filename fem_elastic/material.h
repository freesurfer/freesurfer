
#ifndef H_MATERIAL_H
#define H_MATERIAL_H

#include "small_matrix.h"

//
// this class is an abstract shell
// to get the matrix, create a derived class
//
// this is because the formulae for 2D and 3D are different
class VMaterial
{
public:
  VMaterial()
  {}
  VMaterial(double E, double nu) : m_strLabel()
  {
    m_E=E;
    m_nu=nu;
  }
  virtual ~VMaterial()
  {}

  virtual SmallMatrix get_matrix() const =0;

  double get_E() const
  {
    return m_E;
  }
  double get_nu() const
  {
    return m_nu;
  }

  void set_label(std::string _label)
  {
    m_strLabel = _label;
  }
  std::string label() const
  {
    return m_strLabel;
  }

private:
  double m_E; // Young's constant
  double m_nu; // Poisson

  std::string m_strLabel;
};

#endif // H_MATERIAL_H
