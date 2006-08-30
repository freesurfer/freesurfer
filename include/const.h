/*----------------------------------------------------------------------

      File Name:     const.h

      Description:   miscellaneous constants.

----------------------------------------------------------------------*/

#include <math.h>

#ifndef CONST_H
#define CONST_H

#ifndef True
#define True 1
#define False 0
#endif

#define USE_SC_GSL_REPLACEMENT 1


#define FSMALL          0.00001f
#define MAX_LINE_LEN    4096
#define UNDEFINED       255
#define DEFINED(r, s)   ((r != UNDEFINED) && (s != UNDEFINED))

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif

#ifndef PI
#define PI              M_PI
#endif

#define INV_SQRTPI      0.56419

#define STR_LEN         1024   /* misc. string length */
#define STRLEN          STR_LEN
#define CMD_LINE_LEN    4096

/* predefined file descriptors */
#define FD_STDIN       0
#define FD_STDOUT      1
#define FD_STDERR      2

#endif
