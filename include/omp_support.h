#pragma once

#ifndef HAVE_OMP

#define OMP_FOR_BEGIN
#define OMP_FOR_END

#else

#include "romp_support.h"

#define OMP_PFOR_BEGIN
#define OMP_PFOR_END

#endif
