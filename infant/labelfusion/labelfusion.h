#ifndef LABELFUSION_H
#define LABELFUSION_H

#include "numpy.h"

typedef py::array_t<float>  pyarrayf;
typedef py::array_t<double> pyarrayd;

pyarrayf CMF3D_ML(const pyarrayf &bound, int iNbIters, float fError, float cc, float steps);
pyarrayd performFrontPropagation3D(const pyarrayd &weights, const pyarrayd &start, int max_iters, const pyarrayd &values);

#endif
