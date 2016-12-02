// Created by Eugenio on Nov 7th 2013
// You can set the following two variables to your favorite optimizer for each section of the code
// The first section corresponds to the first iteration, when (most of) the co-registration happens
// The second section corresponds to subsequent iterations, corresponding to (mainly) refinement of the mesh
// The options for each section are CONJUGATE_GRADIENT, LEVENBERG_MARQUARDT and GRADIENT_DESCENT 

#include "kvlOptimizerConstants.h"

#define OPTIMIZER_SECTION_1 CONJUGATE_GRADIENT
#define OPTIMIZER_SECTION_2 GRADIENT_DESCENT

