/* @(#)pr_impl_make.h 1.2 88/02/08 SMI */

#ifndef	pr_impl_make_DEFINED
#define	pr_impl_make_DEFINED

/*
 * Copyright 1987 by Sun Microsystems,  Inc.
 */

/*
 * Pr_devdata is used to keep track of the mmapped virtual
 * address of a device to prevent mapping more than once.
 */
struct pr_devdata {
	struct pr_devdata *next; /* link to next device of this type */
	dev_t	rdev;		/* device type */
	int	count;		/* reference count */
	int	fd;		/* fd of frame buffer, -1 if unused */
	caddr_t	va; 		/* virtual address */
	int	bytes;		/* size of va, 0 for no munmap */
	caddr_t	va2;		/* second virtual address, 0 if unused */
	int	bytes2;		/* second size */
};

#ifndef KERNEL
Pixrect *pr_makefromfd();
Pixrect *pr_makefromfd_2();
#endif !KERNEL

#endif	pr_impl_make_DEFINED
