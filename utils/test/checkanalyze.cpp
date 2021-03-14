/**
 * @brief testing analyze orient handling
 *
 */
/*
 * Original Author: Yasunari Tosa
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

#include <iomanip>
#include <iostream>
// just the hack for RH7.3
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

#include "affine.h"
#include "matrix.h"
#include "mri.h"

const char *Progname = "checkanalyze";

#define V4_LOAD(v, x, y, z, r)                                                 \
  (VECTOR_ELT(v, 1) = x, VECTOR_ELT(v, 2) = y, VECTOR_ELT(v, 3) = z,           \
   VECTOR_ELT(v, 4) = r);

bool compareBad(const MRI *base, const MRI *tmp) {
  VECTOR *vb  = VectorAlloc(4, MATRIX_REAL);
  VECTOR *vt  = VectorAlloc(4, MATRIX_REAL);
  bool    bad = false;

  MATRIX *base_i_to_r = MatrixAlloc(4, 4, MATRIX_REAL);
  GetAffineMatrix(base_i_to_r, base->i_to_r__);

  // RAS voxel position should be the same
  // get the base->tmp voxel map
  MATRIX *b2t = MatrixMultiply(tmp->r_to_i__, base_i_to_r, NULL);
  std::cout << std::endl;
  for (int k = 0; k < base->depth; k++)
    for (int j = 0; j < base->height; ++j)
      for (int i = 0; i < base->width; ++i) {
        V4_LOAD(vb, i, j, k, 1.);
        MatrixMultiply(b2t, vb, vt);
        int it, jt, kt;
        it = int(floor(V3_X(vt) + 0.5));
        jt = int(floor(V3_Y(vt) + 0.5));
        kt = int(floor(V3_Z(vt) + 0.5));
        // now compare
        if (it >= 0 && it < tmp->width && jt >= 0 && jt < tmp->height &&
            kt >= 0 && kt < tmp->depth) {
          if (MRIvox(base, i, j, k) != MRIvox(tmp, it, jt, kt)) {
            if (bad == false)
              bad = true;
            std::cout << std::endl;
            std::cout << "base: (" << i << "," << j << "," << k
                      << ") =" << MRIvox(base, i, j, k) << "    tmp: (" << it
                      << "," << jt << "," << kt
                      << ") = " << MRIvox(tmp, i, j, k) << std::endl;
          }
        }
      }
  std::cout << std::endl;
  return bad;
}

int main(int argc, char *argv[]) {
  std::string file[7];
  file[0] = "./testLAS.img";
  file[1] = "./testLPS.img";
  file[2] = "./testLSA.img";
  file[3] = "./testLIA.img";
  file[4] = "./testASL.img";
  file[5] = "./testAIL.img";
  file[6] = "./testCOR.mgh";
  // first read COR file
  MRI *base = MRIread(const_cast<char *>(file[6].c_str()));
  bool bad  = false;
  for (int i = 0; i < 6; ++i) {
    MRI *tmp = MRIread(const_cast<char *>(file[i].c_str()));
    // verify
    if (compareBad(base, tmp) == true) {
      bad = true;
      std::cout << "Failed comparison for " << file[6].c_str() << " and "
                << file[i].c_str() << std::endl;
    } else
      std::cout << "Good comparison for " << file[6].c_str() << " and "
                << file[i].c_str() << std::endl;
  }
  return (bad == true) ? 1 : 0;
}
