/* @(#)gt_rp.h 1.2 91/03/17 SMI */
#ifndef gt_rp_DEFINED
#define gt_rp_DEFINED

/*
 * @(#)gt_rp.h 1.4 90/08/14 Copyright 1989, Sun Microsystems, Inc.
 */


/*
 * RP Control Space Hardware registers - most of the registers here 
 * are used only by diagnostic and Initialisation. 
 */
/* Host Frame Buffer Address Space */
#define HKPBM_HFBAS_ADDR	0x00090000
#define HKPBM_HFBAS_MASK	0x00000007
#define HKPBM_HFBAS_ST_IMAGE	0x00000000	/* Stencil image */
#define HKPBM_HFBAS_ST_DEPTH	0x00000001	/* Stencil depth */
#define HKPBM_HFBAS_ST_WINDOW	0x00000002	/* Stencil window */
#define HKPBM_HFBAS_CURSOR_BYTE	0x00000003	/* Cursor, image-byte */
#define HKPBM_HFBAS_IMAGE	0x00000004	/* Pixel image */
#define HKPBM_HFBAS_DEPTH	0x00000005	/* Pixel depth */
#define HKPBM_HFBAS_WINDOW	0x00000006	/* Pixel window */
#define HKPBM_HFBAS_IMAGE_DEPTH	0x00000007	/* Pixel image/depth */

/* Host control and status register for PBM */
#define HKPBM_HCS_ADDR		0x00090004
#define HKPBM_HCS_WRITE_MASK	0x00000011	/* Writable bits in register */
#define HKPBM_HCS_HFBSS		0x00000001	/* Host FB state set */
#define HKPBM_HCS_VR		0x00000002	/* Vertical retrace */
#define HKPBM_HCS_RPST		0x00000008	/* RP stalled */
#define HKPBM_HCS_HSRP		0x00000010	/* Host stalling RP */
#define HKPBM_HCS_PBB		0x00000020	/* Pixel formatter busy */
#define HKPBM_HCS_PFB		0x00000040	/* Pixel bus busy */

/* Host FB Address Space reg  - RW*/
typedef struct {
	unsigned	: 29;	/* not used*/
	unsigned	as:3;	/* Access mode and plane group */
} H_rp_as_reg_bit;

typedef union {
    H_rp_as_reg_bit	reg;
    unsigned int	wd;
} H_rp_as_reg;
    

/* Host PBM control space reg - RW. 
 * Note on RP - it can be stalled by FE or Host or RP Semaphore or 
 * by wait for vertical retrace 
 */
typedef struct {
	unsigned	: 25;		/* not used*/
	unsigned	pf_busy:1;		/* Pixel bus busy? RO */
	unsigned	pb_busy:1;		/* Pixel formatter busy? RO */
	unsigned	host_stall:1;		/* Host stalling RP? RW*/
	unsigned	H_rp_stall:1;		/* RP stalled? RO */
	unsigned	:1;			/* not used */
	unsigned	v_retrace:1;	/* Vertical retrace? RO */
	unsigned	fb_stateset:1;	/* Host FB stateset. RW */
} H_rp_csr_reg_bit;

typedef union {
    H_rp_csr_reg_bit	reg;
    unsigned int	wd;
} H_rp_csr_reg;

/* RP registers defined as a structure - from a pixrect view point */

typedef struct {
    H_rp_as_reg	    as_reg;			/* 0 */
    H_rp_csr_reg    csr_reg;			/* 1 */
}H_rp;

/* This is the register to control the front end hold */

#define HKFE_HOLD_REQ		0x1	/* Request for FE hold */
#define HKFE_HOLD_ACK		0x2	/* FE acks hold */

typedef struct {
	unsigned	: 30;		/* not used*/
	unsigned	fe_hold:2;	/* Front end hold/ack bits */
} H_fe_hold_reg_bit;

typedef union {
    H_fe_hold_reg_bit	reg;
    unsigned int	wd;
} H_fe_hold_reg;
    
typedef struct {
    H_fe_hold_reg	    fe_hold_reg;			/* 0 */
}H_fe_hold;

#endif	/* gt_rp_DEFINED */
