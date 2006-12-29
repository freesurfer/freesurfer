/* @(#)pr_util.h 1.18 88/02/08 SMI */

/*
 * Copyright 1983, 1986, 1987 by Sun Microsystems, Inc.
 */

#ifndef pr_util_DEFINED
#define pr_util_DEFINED

/*
 * WARNING:  This include file is obsolete and may disappear in the future.
 *
 * pr_product has been moved to pixrect.h
 * struct pr_devdata etc. has been moved to pr_impl_make.h
 * new loop macros are in pr_impl_util.h
 */

/*
 * Utilities for implementing pixrect operations.
 */

/*
 * Aids to handling overlapping source and destination.
 * Given the from and to pr_pos's, rop_direction tells
 * whether the rasterop is up or down and left or right,
 * encoded as the ROP_UP and ROP_LEFT bits or their absence.
 * The macros rop_is(up|down|left|right) can then be used.
 */
#define ROP_UP  0x1
#define ROP_LEFT 0x2

#define rop_direction(src, so, dst, do) \
     (   ( (((dst).x+(do).x) < ((src).x+(so).x)) << 1) | \
           (((dst).y+(do).y) < ((src).y+(so).y))      )
#define rop_isleft(dir)  ((dir)&ROP_LEFT)
#define rop_isup(dir)  ((dir)&ROP_UP)
#define rop_isright(dir) (((dir)&ROP_LEFT)==0)
#define rop_isdown(dir)  (((dir)&ROP_UP)==0)

/*
 * Aids to producing fast loops, either unrolled or very tight:
 *
 * Cases8(n, op) produces the dense case part of a case statement
 * for the cases [n+1..n+8), repeating ``op'' 1-8 times respectively.
 *
 * Rop_slowloop(n, op) produces a loop to do ``op'' n times, in little space.
 *
 * Rop_fastloop(n, op) produces a loop to do ``op'' n times, in little time.
 *
 * Loop_d6(label, op) produces a dbra loop to do ``op'' the number of times
 * in register d6 (second non-pointer register variable).
 *
 * Loop_d6 is only possible on a 68000 family processor, and rop_fastloop
 * generates an unrolled loop only on a 68010 (assumes other processors
 * will have some kind of I-cache).
 */
#ifdef mc68000
/* generates a nice dbra loop */
#define rop_slowloop(n, op) \
 { register int _loop = (n); \
  if (--_loop >= 0) do { op; } while (--_loop != -1); }

#define loopd6(label, op)      \
 if (0) {       \
  asm("label:");      \
  op;       \
 };        \
 asm("dbra d6,label");
#else /* mc68000*/
#define rop_slowloop(n, op) \
 { register int _loop = (n); \
  while (--_loop >= 0) { op; } }
#endif /* mc68000*/

#ifdef mc68010
#define cases8(n, op)       \
     case (n)+8: op; case (n)+7: op; case (n)+6: op;  \
     case (n)+5: op; case (n)+4: op; case (n)+3: op;  \
     case (n)+2: op; case (n)+1: op;    \

#define rop_fastloop(n, op) \
 { register int _loop ; \
  for (_loop = (n); _loop > 15; _loop -= 16) \
       { op; op; op; op; op; op; op; op; \
           op; op; op; op; op; op; op; op; } \
    switch (_loop) { \
   cases8(8, op); \
   cases8(0, op); \
   case 0: break; \
  } }
#else /*mc68010*/
#define rop_fastloop rop_slowloop
#endif /*mc68010*/

/*
 * Alloctype(datatype) allocates a datatype structure using calloc
 * with the appropriate type cast.
 */
#define alloctype(datatype)      \
     (datatype *)calloc(1, sizeof (datatype))

/*
 * Pr_product is used when doing multiplications involving pixrects,
 * and casts its arguments to that the compiler will use 16 by 16 multiplies.
 */
#ifndef pr_product
#ifdef sun
#define pr_product(a, b) ((short)(a) * (short)(b))
#else
#define pr_product(a, b) ((a) * (b))
#endif
#endif

/*
 * Pr_area is the area of a rectangle.
 */
#define pr_area(size) pr_product((size).x, (size).y)

/*
 * Pr_devdata is used to keep track of the valloced/mmapped virtual
 * address of a device to prevent doing it more than necessary.
 */
struct pr_devdata
{
  struct pr_devdata *next; /* link to next device of this type */
  dev_t rdev;  /* device type */
  int count;  /* reference count */
  int fd;  /* fd of frame buffer, -1 if unused */
  short  *va;   /* virtual address */
  int bytes;  /* size of va, 0 for no munmap */
  caddr_t va2;  /* second virtual address, 0 if unused */
  int bytes2;  /* second size */
};

#ifndef KERNEL
Pixrect *pr_makefromfd();
#endif /*!KERNEL*/

#endif /*pr_util_DEFINED*/
