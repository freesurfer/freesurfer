
#ifndef __h_untangler_h
#define __h_untangler_h

#include "fem_3d.h"

void
solve_topology_problems( TMesh3d& mesh,
                         unsigned int neighborhood=3,
                         bool dbgMode = false);

#endif
