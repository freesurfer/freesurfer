/* @(#)gp1reg.h 1.11 89/03/14 SMI */

/*
 * Copyright 1985, 1987 by Sun Microsystems, Inc.
 */

#ifndef		_sundev_gp1reg_h
#define	_sundev_gp1reg_h

#define	GP1_SIZE		65536
#define	GP1_SHMEM_SIZE		32768
#define	GP1_SHMEM_OFFSET	32768
#define	GP2_SHMEM_SIZE		0x20000

#define	GP2_SIZE		262144

#define	GP1_BLOCK_SHORTS	512
#define	GP1_BLOCKS		32
#define	GP1_CONTROL_BLOCKS	1
#define	GP1_POST_BLOCKS		23
#define	GP1_STATIC_BLOCKS	8

#define	GP1_POST_SLOTS		256

/* GP1 board layout */
struct gp1 {
	union gp1_reg {
		struct gp1reg {
			short gpr_ident;
			union {
				short status;
				short control;
			} gpr_csr;
			short gpr_ucode_addr;
			short gpr_ucode_data;
		} reg;
		struct gp2reg {
			short ident;
			short :16;
			u_int xp_addr;
			u_int xp_data_h;
			u_int xp_data_l;
			u_int rp_addr;
			u_int rp_data;
			u_int pp_addr;
			u_int pp_data_h;
			u_int pp_data_l;
			u_int status;
			u_int control;
			u_int pp_addr2;
		} reg2;
		u_short ureg[16384];
	} reg;
	struct gp1_shmem {
		u_short	host_count;
		u_short gp_count;
		u_char ver_flag;
		u_char alloc_sem;
		u_short host_alloc_h;
		u_short host_alloc_l;
		u_short gp_alloc_h;
		u_short gp_alloc_l;
		u_char flag2;
		u_char post_sem;
		u_short post_buf[GP1_POST_SLOTS];
		u_char block_owners[GP1_BLOCKS];
		u_short fill280[228];
		u_short ver_release;
		u_short ver_serial;
		u_short fill510[2];
		u_short post[GP1_POST_BLOCKS][GP1_BLOCK_SHORTS];
		u_short stat[GP1_STATIC_BLOCKS][GP1_BLOCK_SHORTS];
	} shmem;
};

/* convert block number to bit mask, owner array index */
#define	GP1_ALLOC_BIT(blk)	((u_long) 0x80000000 >> (blk))
#define	GP1_OWNER_INDEX(blk)	(GP1_BLOCKS - 1 - (blk))

/* register offsets and bit definitions */

#define	GP1_BOARD_IDENT_REG	0
#define	GP1_ID_MASK		0xfe
#define	GP1_ID_VALUE		0xea

#define	GP1_CONTROL_REG		1
#define	GP1_CR_CLRIF		0x8000
#define	GP1_CR_IENBLE		0x0300
#define	GP1_CR_INT_TOGGLE	0x0300
#define	GP1_CR_INT_DISABLE	0x0200
#define	GP1_CR_INT_ENABLE	0x0100
#define	GP1_CR_INT_NOCHANGE	0x0000
#define	GP1_CR_RESET		0x0040
#define	GP1_CR_VP_CONTROL	0x0038
#define	GP1_CR_VP_STRT0		0x0020
#define	GP1_CR_VP_HLT		0x0010
#define	GP1_CR_VP_CONT		0x0008
#define	GP1_CR_PP_CONTROL	0x0007
#define	GP1_CR_PP_STRT0		0x0004
#define	GP1_CR_PP_HLT		0x0002
#define	GP1_CR_PP_CONT		0x0001

#define	GP1_STATUS_REG		1
#define	GP1_SR_IFLG		0x8000
#define	GP1_SR_IEN		0x4000
#define	GP1_SR_RESET		0x0400
#define	GP1_SR_VP_STATE		0x0200
#define	GP1_SR_PP_STATE		0x0100
#define	GP1_SR_VP_STATUS	0x00F0
#define	GP1_SR_PP_STATUS	0x000F

#define	GP1_UCODE_ADDR_REG	2
#define	GP1_UCODE_DATA_REG	3


/* GP2 definitions */
#define	GP2_ID_MASK		0xff
#define	GP2_ID_VALUE		0xec

/* GP2 control register bits (active low) */
#define	GP2_CR_XP_RST		0x40000000
#define	GP2_CR_XP_HLT		0x20000000
#define	GP2_CR_RP_RST		0x10000000
#define	GP2_CR_RP_HLT		0x08000000
#define	GP2_CR_PP_RST		0x02000000
#define	GP2_CR_PP_HLT		0x01000000

/* GP2 status register message shifts */
#define	GP2_SR_XPMSG_SHIFT	16
#define	GP2_SR_RPMSG_SHIFT	8
#define	GP2_SR_PPMSG_SHIFT	0


/* shared memory offsets */
#define	GP1CB_HOST_COUNT	0	/* host command count */
#define	GP1CB_GP_COUNT		1	/* GP command count */
#define	GP1CB_ALLOC_SEM		2	/* buffer allocation semaphore */
#define	GP1CB_HOST_ALLOC_H	3	/* host block allocation bit vector */
#define	GP1CB_HOST_ALLOC_L	4	/*	high/low */
#define	GP1CB_GP_ALLOC_H	5	/* GP block allocation bit vector */
#define	GP1CB_GP_ALLOC_L	6	/*	high/low */
#define	GP1CB_POST_SEM		7	/* command posting semaphore */
#define	GP1CB_POST_BUF		8	/* command posting buffer, 256 words */
#define	GP1CB_BLOCK_OWNERS	264	/* minor device number for owner */
					/* of each block -- reverse order! */

#define	GP1CB_VER_FLAG		2	/* version info is valid if */
					/*	high byte is 1 */
#define	GP1CB_VER_RELEASE	508	/* major/minor release numbers */
#define	GP1CB_VER_SERIAL	509	/* serial number/version flag */

#endif	/*  _sundev_gp1reg_h  */
