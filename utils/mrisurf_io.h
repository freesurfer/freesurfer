#pragma once

#include "mrisurf_deform.h"

int mrisReadTransform(MRIS *mris, const char *mris_fname);

// STL
MRIS* mrisReadSTL(const std::string &filename);
bool mrisWriteSTL(const MRIS* mris, const std::string &filename);
