/* @(#)mem_rop_impl_ops.h 1.4 88/02/08 SMI */

/*
 * Copyright 1986 by Sun Microsystems,  Inc.
 */

#ifndef	mem_rop_impl_ops_DEFINED
#define	mem_rop_impl_ops_DEFINED

/*
 * rasterop function macros
 */

/*
 * unmasked rop
 */
#define	OP_urop0(d,s)	(0)
#define	OP_urop1(d,s)	(~((d) | (s)))
#define	OP_urop2(d,s)	((d) & ~(s))
#define	OP_urop3(d,s)	(~(s))
#define	OP_urop4(d,s)	(~(d) & (s))
#define	OP_urop5(d,s)	(~(d))
#define	OP_urop6(d,s)	((d) ^ (s))
#define	OP_urop7(d,s)	(~((d) & (s)))
#define	OP_urop8(d,s)	((d) & (s))
#define	OP_urop9(d,s)	((d) ^ ~(s))
#define	OP_uropA(d,s)	(d)
#define	OP_uropB(d,s)	((d) | ~(s))
#define	OP_uropC(d,s)	(s)
#define	OP_uropD(d,s)	(~(d) | (s))
#define	OP_uropE(d,s)	((d) | (s))
#define	OP_uropF(d,s)	(~0)

/*
 * masked rop
 *
 * correct mask:	m = OP_mmsk(mask)
 * perform rop:		d = OP_mrop(dst, mask, src)
 */
#define	OP_mmsk0(m)	(~(m))
#define	OP_mrop0(d,m,s)	((d) & (m))

#define	OP_mmsk1(m)	(m)
#define	OP_mrop1(d,m,s)	((d) ^ ((m) & ((d) | ~(s))))

#define	OP_mmsk2(m)	(~(m))
#define	OP_mrop2(d,m,s)	((d) & ((m) | ~(s)))

#define	OP_mmsk3(m)	(m)
#define	OP_mrop3(d,m,s)	((d) ^ ((m) & ((d) ^ ~(s))))

#define	OP_mmsk4(m)	(m)
#define	OP_mrop4(d,m,s)	((d) ^ ((m) & ((d) | (s))))

#define	OP_mmsk5(m)	(m)
#define	OP_mrop5(d,m,s)	((d) ^ (m))

#define	OP_mmsk6(m)	(m)
#define	OP_mrop6(d,m,s)	((d) ^ ((m) & (s)))

#define	OP_mmsk7(m)	(m)
#define	OP_mrop7(d,m,s)	((d) ^ ((m) & (~(d) | (s))))

#define	OP_mmsk8(m)	(~(m))
#define	OP_mrop8(d,m,s)	((d) & ((m) | (s)))

#define	OP_mmsk9(m)	(m)
#define	OP_mrop9(d,m,s)	((d) ^ ((m) & ~(s)))

#define	OP_mmskA(m)	(m)
#define	OP_mropA(d,m,s)	(d)

#define	OP_mmskB(m)	(m)
#define	OP_mropB(d,m,s)	((d) | ((m) & ~(s)))

#define	OP_mmskC(m)	(m)
#define	OP_mropC(d,m,s)	((d) ^ ((m) & ((d) ^ (s))))

#define	OP_mmskD(m)	(m)
#define	OP_mropD(d,m,s)	((d) ^ ((m) & ~((d) & (s))))

#define	OP_mmskE(m)	(m)
#define	OP_mropE(d,m,s)	((d) | ((m) & (s)))

#define	OP_mmskF(m)	(m)
#define	OP_mropF(d,m,s)	((d) | (m))


/*
 * unmasked fill
 *
 * generate fill constant:	k = OP_ufgen(color)
 * perform unmasked fill:	d = OP_ufill(dst, constant)
 */
#define	OP_ufgen0(c)	(0)
#define	OP_ufill0(d,k)	(k)

#define	OP_ufgen1(c)	(c)
#define	OP_ufill1(d,k)	(~((d) | (k)))

#define	OP_ufgen2(c)	(~(c))
#define	OP_ufill2(d,k)	((d) & (k))

#define	OP_ufgen3(c)	(~(c))
#define	OP_ufill3(d,k)	(k)

#define	OP_ufgen4(c)	(c)
#define	OP_ufill4(d,k)	(~(d) & (k))

#define	OP_ufgen5(c)	(c)
#define	OP_ufill5(d,k)	(~(d))

#define	OP_ufgen6(c)	(c)
#define	OP_ufill6(d,k)	((d) ^ (k))

#define	OP_ufgen7(c)	(c)
#define	OP_ufill7(d,k)	(~((d) & (k)))

#define	OP_ufgen8(c)	(c)
#define	OP_ufill8(d,k)	((d) & (k))

#define	OP_ufgen9(c)	(~(c))
#define	OP_ufill9(d,k)	((d) ^ (k))

#define	OP_ufgenA(c)	(c)
#define	OP_ufillA(d,k)	(d)

#define	OP_ufgenB(c)	(~(c))
#define	OP_ufillB(d,k)	((d) | (k))

#define	OP_ufgenC(c)	(c)
#define	OP_ufillC(d,k)	(k)

#define	OP_ufgenD(c)	(c)
#define	OP_ufillD(d,k)	(~(d) | (k))

#define	OP_ufgenE(c)	(c)
#define	OP_ufillE(d,k)	((d) | (k))

#define	OP_ufgenF(c)	(~0)
#define	OP_ufillF(d,k)	(k)


/*
 * masked fill
 *
 * generate fill constant:	k = OP_mfgen(mask, color)
 * generate mask for fill:	m = OP_mfmsk(mask)
 * perform masked fill:		d = OP_mfill(dst, mask, constant)
 */
#define	OP_mfgen0(m,c)	(c)
#define	OP_mfmsk0(m)	(~(m))
#define	OP_mfill0(d,m,k) ((d) & (m))

#define	OP_mfgen1(m,c)	(~(c))
#define	OP_mfmsk1(m)	(m)
#define OP_mfill1(d,m,k) ((d) ^ ((m) & ((d) | (k))))

#define	OP_mfgen2(m,c)	(~(m) | ~(c))
#define	OP_mfmsk2(m)	(m)
#define	OP_mfill2(d,m,k) ((d) & (k))

#define	OP_mfgen3(m,c)	((m) & ~(c))
#define	OP_mfmsk3(m)	(~(m))
#define	OP_mfill3(d,m,k) (((d) & (m)) | (k))	/* needs work */

#define	OP_mfgen4(m,c)	(c)
#define	OP_mfmsk4(m)	(m)
#define	OP_mfill4(d,m,k) ((d) ^ ((m) & ((d) | (k))))

#define	OP_mfgen5(m,c)	(c)
#define	OP_mfmsk5(m)	(m)
#define	OP_mfill5(d,m,k) ((d) ^ (m))

#define	OP_mfgen6(m,c)	((m) & (c))
#define	OP_mfmsk6(m)	(m)
#define	OP_mfill6(d,m,k) ((d) ^ (k))

#define	OP_mfgen7(m,c)	(c)
#define	OP_mfmsk7(m)	(m)
#define	OP_mfill7(d,m,k) ((d) ^ ((m) & (~(d) | (k))))

#define	OP_mfgen8(m,c)	(~(m) | (c))
#define	OP_mfmsk8(m)	(m)
#define	OP_mfill8(d,m,k) ((d) & (k))

#define	OP_mfgen9(m,c)	((m) & ~(c))
#define	OP_mfmsk9(m)	(m)
#define	OP_mfill9(d,m,k) ((d) ^ (k))

#define	OP_mfgenA(m,c)	(c)
#define	OP_mfmskA(m)	(m)
#define	OP_mfillA(d,m,k) (d)

#define	OP_mfgenB(m,c)	((m) & ~(c))
#define	OP_mfmskB(m)	(m)
#define	OP_mfillB(d,m,k) ((d) | (k))

#define	OP_mfgenC(m,c)	((m) & (c))
#define	OP_mfmskC(m)	(~(m))
#define	OP_mfillC(d,m,k) (((d) & (m)) | (k))	/* needs work */

#define	OP_mfgenD(m,c)	(~(c))
#define	OP_mfmskD(m)	(m)
#define	OP_mfillD(d,m,k) ((d) ^ ((m) & (~(d) | (k))))

#define	OP_mfgenE(m,c)	((m) & (c))
#define	OP_mfmskE(m)	(m)
#define	OP_mfillE(d,m,k) ((d) | (k))

#define	OP_mfgenF(m,c)	(c)
#define	OP_mfmskF(m)	(m)
#define	OP_mfillF(d,m,k) ((d) | (m))


/*
 * unmasked fill, color = 0
 * (pcc is too dumb for this)
 */
#define	OP_uzero0(d)	(0)
#define	OP_uzero1(d)	(~(d))
#define	OP_uzero2(d)	(d)
#define	OP_uzero3(d)	(~0)
#define	OP_uzero4(d)	(0)
#define	OP_uzero5(d)	(~(d))
#define	OP_uzero6(d)	(d)
#define	OP_uzero7(d)	(~0)
#define	OP_uzero8(d)	(0)
#define	OP_uzero9(d)	(~(d))
#define	OP_uzeroA(d)	(d)
#define	OP_uzeroB(d)	(~0)
#define	OP_uzeroC(d)	(0)
#define	OP_uzeroD(d)	(~(d))
#define	OP_uzeroE(d)	(d)
#define	OP_uzeroF(d)	(~0)

/*
 * masked fill, color = 0
 *
 * generate mask for fill:	m = OP_mfmsk(mask)
 * perform masked zero fill:	d = OP_mzero(dst, mask)
 */
#define	OP_mzero0(d,m)	((d) & (m))
#define OP_mzero1(d,m)	((d) ^ (m))
#define	OP_mzero2(d,m)	(d)
#define	OP_mzero3(d,m)	((d) | ~(m))	/* problem */
#define	OP_mzero4(d,m)	((d) & ~(m))	/* problem */
#define	OP_mzero5(d,m)	((d) ^ (m))
#define	OP_mzero6(d,m)	(d)
#define	OP_mzero7(d,m)	((d) | (m))
#define	OP_mzero8(d,m)	((d) & (m))
#define	OP_mzero9(d,m)	((d) ^ (m))
#define	OP_mzeroA(d,m)	(d)
#define	OP_mzeroB(d,m)	((d) | (m))
#define	OP_mzeroC(d,m)	((d) & (m))
#define	OP_mzeroD(d,m)	((d) ^ (m))
#define	OP_mzeroE(d,m)	(d)
#define	OP_mzeroF(d,m)	((d) | (k))

#endif
