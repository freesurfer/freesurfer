/*----------------------------------------------------------------------

      File Name:     const.h

      Description:   miscellaneous constants.

----------------------------------------------------------------------*/

#ifndef CONST_H
#define CONST_H

#ifndef True
#define True 1
#define False 0
#endif

#define MAX_LINE_LEN    2000
#define UNDEFINED       255
#define DEFINED(r, s)   ((r != UNDEFINED) && (s != UNDEFINED))

#ifndef PI
#define PI             3.1415926
#endif

#define INV_SQRTPI     0.56419

#define STR_LEN        100   /* misc. string length */

#endif
