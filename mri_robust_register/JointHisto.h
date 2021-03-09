/**
 * @brief A class for a joint histogram of two images
 *
 */

/*
 * Original Author: Martin Reuter
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
//
// written by Martin Reuter
// Sep. 1st ,2011
//
#ifndef JointHisto_H
#define JointHisto_H

#include "mri.h"
#include "mriBSpline.h"

#define export // obsolete feature 'export template' used in these headers 
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_matlab_print.h>
#undef export

#include <vcl_iostream.h>

class JointHisto
{
public:

  JointHisto() :
      n(256), sum(0.0), histo(256, 256), haseps(false)
  {
  }
  ;
  JointHisto(MRI *mri1, MRI * mri2, int d1 = 1, int d2 = 1, int d3 = 1) :
      n(256), sum(0.0), histo(256, 256), haseps(false)
  {
    create(mri1, mri2, d1, d2, d3);
  }
  ;
  JointHisto(MRI *mri1, MRI * mri2, const vnl_matrix_fixed<double, 4, 4>& M1,
      const vnl_matrix_fixed<double, 4, 4>& M2, int d1 = 1, int d2 = 1, int d3 =
          1) :
      n(256), sum(0.0), histo(256, 256), haseps(false)
  {
    create(mri1, mri2, M1, M2, d1, d2, d3);
  }
  ;
  void create(MRI *mri1, MRI * mri2, int d1 = 1, int d2 = 1, int d3 = 1);
  void create(MRI *mri1, MRI * mri2, const vnl_matrix_fixed<double, 4, 4>& M1,
      const vnl_matrix_fixed<double, 4, 4>& M2, int d1 = 1, int d2 = 1, int d3 =
          1);
  void create(MRI_BSPLINE *bspline1, MRI_BSPLINE * bspline2,
      const vnl_matrix_fixed<double, 4, 4>& M1,
      const vnl_matrix_fixed<double, 4, 4>& M2, int d1 = 1, int d2 = 1, int d3 =
          1);
  MRI * locate(MRI *mri1, MRI * mri2, const vnl_matrix_fixed<double, 4, 4>& M1,
      const vnl_matrix_fixed<double, 4, 4>& M2, int d1 = 1, int d2 = 1, int d3 =
          1);
  void set(const vnl_matrix<double> & histo);
  void smooth(double fwhm1 = 7.0);
  void print(const std::string & n = "H")
  {
    vnl_matlab_print(vcl_cout,histo,n.c_str());std::cout << std::endl;};
  void save(const std::string & fname,const std::string & n = "H");
  double clip (double thres);
  void normalize()
  { if (sum == 0.0 || sum == 1.0) return;
    else
    { histo /= sum; sum = 1.0;}};
  std::pair < double , double> getMinMax();

  double computeMI();  // mutual information
        double computeNMI();// normalized mutual information
        double computeECC();
        double computeNCC();
        double computeLS();
        double computeSCR();// symmetric correlation ratio

      protected:

        void computeRCsums()
        {
          // compute row and column sums
          vnl_matrix < double > v1(n,1,1);
          vnl_matrix < double > v2(1,n,1);
          if ((int)rowsum.size() != n) rowsum.set_size(n);
          if ((int)colsum.size() != n) colsum.set_size(n);
          rowsum = (histo *v1).get_column(0);
          colsum = (v2 * histo).get_row(0);
        };

        void addeps(double eps=2.2204E-16)
        { if (haseps) return;
          else
          { histo += eps;haseps=true;}};

        int n;
        double sum;
        vnl_matrix < double > histo;
        vnl_vector < double > rowsum;
        vnl_vector < double > colsum;
        bool haseps;
      };

#endif
