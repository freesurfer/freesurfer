/* @(#)mem_rop_impl_util.h 1.7 89/05/19 SMI */

/*
 * Copyright 1986, 1987 by Sun Microsystems,  Inc.
 */

#ifndef	mem_rop_impl_util_DEFINED
#define	mem_rop_impl_util_DEFINED

/*
 * Utility macros for memory pixrect code
 */

/* 
 * Code selection macros 
 */
#define	IFTRUET(t, a, b)	_STMT(if (t) {a;} else {b;})
#define	IFFALSET(t, a, b)	_STMT(b;)

/*
 * Option control macros
 */
#define ALWAYS(option,a,b)	IFTRUE(a,b)
#define NEVER(option,a,b)	IFFALSE(a,b)
#define OPTION(option,a,b)	option(a,b)

#define	ALWAYST(option,t,a,b)	IFTRUET(t,a,b)
#define	NEVERT(option,t,a,b)	IFFALSET(t,a,b)
#define	OPTIONT(option,t,a,b)	option(t,a,b)

/*
 * Speed choice macros
 *
 * The T form is used when the "never" half of the fast code
 * is the same as the slow code.
 *
 * Define NEVER_SLOW to get all fast code.
 * Define NEVER_FAST to get all slow code.
 */
 
#if defined(NEVER_FAST) && !defined(NEVER_SLOW)
#define	FAST	IFFALSE
#define	FASTT	IFFALSET
#else
#define	FAST	IFTRUE
#define	FASTT	IFTRUET
#endif

#ifdef NEVER_SLOW
#define	SLOW	IFTRUE
#define	SLOWT	IFTRUET
#else
#define	SLOW	IFFALSE
#define	SLOWT	IFFALSET
#endif

/* fast in user, slow in kernel */
#define UFAST IFKERNEL(SLOW,FAST)
#define UFASTT IFKERNEL(SLOWT,FASTT)

#endif	mem_rop_impl_util_DEFINED
