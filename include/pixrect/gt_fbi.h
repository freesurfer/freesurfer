/* @(#)gt_fbi.h 1.4 91/05/31 SMI */
#ifndef gt_fbi_DEFINED
#define gt_fbi_DEFINED

/*
 * Copyright 1989, Sun Microsystems, Inc.
 */

/* Hawk FB Input control space register definitions */

/* Copy destination address register */
#define HKPF_CDA_ADDR		0x00100000	/* Copy destination */
#define HKPF_CDA_MASK		0x00FFFFFF
#define HKPF_CDA_PG_MASK	0x00C00000	/* Plane group */
#define HKPF_CDA_PG_IMAGE	0x00000000
#define HKPF_CDA_PG_DEPTH	0x00400000
#define HKPF_CDA_PG_WINDOW	0x00800000
#define HKPF_CDA_PG_IMAGE_DEPTH 0x00C00000
#define HKPF_CDA_Y_MASK		0x003FF800
#define HKPF_CDA_X_MASK		0x000007FF
#define HKPF_CDA_XY_MASK	0x003FFFFF
#define HKPF_CDA_Y_SHIFT	11
#define HKPF_CDA_X_SHIFT	0

/* Fill destination address register */
/* Note: COPY_LEFT must be 0 in copy/fill direction register! */
#define HKPF_FDA_ADDR		0x00100004	/* Fill source */
#define HKPF_FDA_MASK		0x00FFFFFF
#define HKPF_FDA_PG_MASK	0x00C00000	/* Plane group */
#define HKPF_FDA_PG_IMAGE	0x00000000
#define HKPF_FDA_PG_DEPTH	0x00400000
#define HKPF_FDA_PG_WINDOW	0x00800000
 /* Image_depth is illegal */
#define HKPF_FDA_Y_MASK		0x003FF800
#define HKPF_FDA_X_MASK		0x000007FF
#define HKPF_FDA_XY_MASK	0x003FFFFF
#define HKPF_FDA_Y_SHIFT	11
#define HKPF_FDA_X_SHIFT	0
/* View clip registers */
#define HKPP_VCX_ADDR		0x0020008C
#define HKPP_VCX_MASK		0x07FF07FF	/* Mask, then shift */
#define HKPP_VCY_ADDR		0x00200094
#define HKPP_VCY_MASK		0x07FF07FF
#define HKPP_VCXY_MASK		0x000007FF
#define HKPP_VCX_LEFT_SHIFT	16
#define HKPP_VCX_RIGHT_SHIFT	0
#define HKPP_VCY_TOP_SHIFT	16
#define HKPP_VCY_BOTTOM_SHIFT	0

/* Pick aperture registers */
#define HKPP_PAX_ADDR		0x00200098	/* X left and right */
#define HKPP_PAX_MASK		0x07FF07FF
#define HKPP_PAY_ADDR		0x002000A4	/* Y top and bottom */
#define HKPP_PAY_MASK		0x07FF07FF
#define HKPP_PAF_ADDR		0x002000A8	/* Front */
#define HKPP_PAB_ADDR		0x002000B0	/* Back */
#define HKPP_PAXY_MASK		0x000007FF
#define HKPP_PAX_LEFT_SHIFT	16
#define HKPP_PAX_RIGHT_SHIFT	0
#define HKPP_PAY_TOP_SHIFT	16
#define HKPP_PAY_BOTTOM_SHIFT	0
#define HKPP_PAZ_MASK		0x00FFFFFF	/* Front or back */

/* Pick control register */
#define HKPP_PC_ADDR		0x002000C4
#define HKPP_PC_MASK		0x00000007
#define HKPP_PC_PICK_ENA	0x00000004
#define HKPP_PC_PICK_WHL_RENDER 0x00000002
#define HKPP_PC_PICK_3D_APER	0x00000001

/* Frame buffer width register */
#define HKPP_FBW_ADDR		0x002000C8
#define HKPP_FBW_MASK		0x00000003
#define HKPP_FBW_1280		0x00000000	/* 1280 (x 1024) */
#define HKPP_FBW_1920		0x00000001	/* 1920 (x 1035) */
#define HKPP_FBW_960		0x00000002	/* 960 (x 680) */
#define HKPP_FBW_960STEREO	0x00000003	/* 960 (x 680) Stereo */

/* Stereo write enable */
#define HKPP_SWE_ADDR		0x002000D0
#define HKPP_SWE_MASK		0x00000003
#define HKPP_SWE_MONO		0x00000000	/* 0 = mono */
#define HKPP_SWE_STEREO		0x00000001	/* 1 = stereo */
#define HKPP_SWE_RIGHT		0x00000010	/* 2 = right eye */


/* Interleave test register */
#define HKPP_IT_ADDR		0x002000E0
#define HKPP_IT_MASK		0x0000000F
#define HKPP_IT_BANK_ID		0x00000007
#define HKPP_IT_TEST_ENABLE	0x00000008

/* Screen door transparency matrix (Note start address not on even boundary) */
#define HKPP_SDTM_ADDR		0x00200104	/* 20 registers +0x04 apart */
#define HKPP_SDTM_ADDR_MASK	0x0000007C
#define HKPP_SDTM_MASK		0x0000FFFF
#define HKPP_SDTM_COLUMNS	20     /* Width of screen door transp matrix */
#define HKPP_SDTM_BITS		16     /* Number of bits for sdtm word */

/* Screen door transparency enable */
#define HKPP_SDTE_ADDR		0x00200180
#define HKPP_SDTE_MASK		0x00000001
#define HKPP_SDTE_ENA		0x00000001	/* Use screen door pattern */

/* Buffer select register */
#define HKPP_BS_ADDR		0x00200204
#define HKPP_BS_MASK	 	0x0000003F
#define HKPP_BS_I_WRITE_A	0x00000001	/* Write image A */
#define HKPP_BS_I_WRITE_B	0x00000002	/* (Both may be on) */
#define HKPP_BS_I_READ_B	0x00000004	/* 0 = read image A */
#define HKPP_BS_O_WRITE_A 	0x00000001	/* Write overlay A */
#define HKPP_BS_O_WRITE_B	0x00000002	/* Both may be on */
#define HKPP_BS_O_READ_B	0x00000004	/* 0 = read overlay A */

/* Fast clear control register */
#define HKPP_FCC_ADDR		0x00200210
#define HKPP_FCC_MASK		0x00000003
#define HKPP_FCC_NORMAL		0x00000000
#define HKPP_FCC_INITIALIZE	0x00000002
#define HKPP_FCC_FILL		0x00000003

/* Current window ID register */
#define HKPP_CWID_ADDR		0x00200400
#define HKPP_CWID_MASK		0x000003FF

/* Window ID control register */
#define HKPP_WIDC_ADDR		0x00200404
#define HKPP_WIDC_MASK		0x800003FF
#define HKPP_WIDC_REPLACE_WID	0x80000000	/* Enable WID replace mode */
#define HKPP_WIDC_WID_MASK	0x000003FF
/* A "1" means the bit is compared */

/* Window Background register */
#define HKPP_WBG_ADDR		0x00200408
/* No mask.  32 bits of alpha, blue, green and red. */

/* Constant Z source */
#define HKPP_CZS_ADDR		0x0020040C
#define HKPP_CZS_MASK		0x00FFFFFF

/* Z control */
#define HKPP_ZC_ADDR		0x00200410
#define HKPP_ZC_MASK		0x0000000F
#define HKPP_ZC_WID_EXT		0x00000001	/* WID extension clip enable */
#define HKPP_ZC_CONST_Z		0x00000002	/* Constant Z enable */
#define HKPP_ZC_WRITE_ENA	0x00000004	/* Z write enable */
#define HKPP_ZC_COMPARE		0x00000008	/* HSR enable (compare) */

/* Image write mask register */
#define HKPP_IWM_ADDR		0x00200414
/* No mask.  32 bits of alpha, blue, green and red. */

/* Window write mask register */
#define HKPP_WWM_ADDR		0x00200418
#define HKPP_WWM_MASK		0x00001FFF
#define HKPP_WWM_FC		0x00001000	/* Fast clear */
#define HKPP_WWM_CD		0x00000800	/* Cursor data */
#define HKPP_WWM_CE		0x00000400	/* Cursor enable */
#define HKPP_WWM_WID		0x000003FF	/* Window ID planes */

/* Byte mode channel select */
#define HKPP_BMCS_ADDR		0x0020041C
#define HKPP_BMCS_MASK		0x00000003
#define HKPP_BMCS_ALPHA_OV	0x00000000
#define HKPP_BMCS_BLUE		0x00000001
#define HKPP_BMCS_GREEN		0x00000002
#define HKPP_BMCS_RED		0x00000003

/* Byte mode write mask register */
#define HKPP_BMWM_ADDR		0x00200420
#define HKPP_BMWM_MASK		0x000000FF
/* ROP/Blend operation register */
#define HKPP_RBO_ADDR		0x00200424
#define HKPP_RBO_MASK		0x0000001F
#define HKPP_RBO_BLEND		0x00000010	/* Enable blending (not ROP) */
#define HKPP_RBO_ROP_MASK	0x0000000F
/* ROP definitions */
#define HKPP_ROP_ZERO		0x00
#define HKPP_ROP_NSRC_AND_NDST	0x01
#define HKPP_ROP_NSRC_AND_DST	0x02
#define HKPP_ROP_NSRC		0x03
#define HKPP_ROP_SRC_AND_NDST	0x04
#define HKPP_ROP_NDST		0x05
#define HKPP_ROP_SRC_XOR_DST	0x06
#define HKPP_ROP_NSRC_OR_NDST	0x07
#define HKPP_ROP_SRC_AND_DST	0x08
#define HKPP_ROP_NSRC_XOR_DST	0x09
#define HKPP_ROP_DST		0x0A
#define HKPP_ROP_NSRC_OR_DST	0x0B
#define HKPP_ROP_SRC		0x0C   /* Direct write */
#define HKPP_ROP_SRC_OR_NDST	0x0D
#define HKPP_ROP_SRC_OR_DST	0x0E
#define HKPP_ROP_ONE		0x0F

/* Mpg and fcs control register */
#define HKPP_MPG_ADDR           0x00200434
#define HKPP_MPG_MASK           0x000000FF
#define HKPP_MPG_FCS_MASK       0x00000007
#define HKPP_MPG_FCS_DISABLE    0x00000004
#define HKPP_MPG_WRITE_EN_MASK  0x000000F8
#define HKPP_MPG_WRITE_RGB      0x000000E0
#define HKPP_MPG_WRITE_OV       0x00000010
#define HKPP_MPG_WRITE_Z        0x00000008
 
/* Stencil mask */
#define HKPF_SM_ADDR		0x00200C00	/* All 32 bits */

/* Stencil foreground color */
#define HKPF_SFC_ADDR		0x00200C04	/* All 32 bits */

/* Stencil background color */
#define HKPF_SBC_ADDR		0x00200C08	/* All 32 bits */

/* Stencil transparency flag */
#define HKPF_STF_ADDR		0x00200C0C
#define HKPF_STF_MASK		0x00000001
#define HKPF_STF_TRANSPARENT	0x00000001	/* Transparent */
#define HKPF_STF_USE_BG		0x00000000	/* Use background color */

/* Copy/fill direction/size register */
#define HKPF_CFDS_ADDR		0x00200C10
#define HKPF_CFDS_MASK		0x803FFFFF
#define HKPF_CFDS_COPY_LEFT	0x80000000	/* Copy right to left */
#define HKPF_CFDS_COPY_HEIGHT_MASK 0x003FF800	/* Pixel count up and down */
#define HKPF_CFDS_COPY_WIDTH_MASK 0x000007FF	/* Pixel count across */
#define HKPF_CFDS_XY_MASK	0x003FFFFF
#define HKPF_CFDS_COPY_HEIGHT_SHIFT 11
#define HKPF_CFDS_COPY_WIDTH_SHIFT 0

/* Copy source address registers */
#define HKPF_CSA_ADDR		0x00200C14	/* Copy source */
#define HKPF_CSA_MASK		0x00FFFFFF
#define HKPF_CSA_PG_MASK	0x00C00000	/* Plane group */
#define HKPF_CSA_PG_IMAGE	0x00000000
#define HKPF_CSA_PG_DEPTH	0x00400000
#define HKPF_CSA_PG_WINDOW	0x00800000
#define HKPF_CSA_PG_IMAGE_DEPTH	0x00C00000
#define HKPF_CSA_Y_MASK		0x003FF800
#define HKPF_CSA_X_MASK		0x000007FF
#define HKPF_CSA_XY_MASK	0x003FFFFF
#define HKPF_CSA_Y_SHIFT	11
#define HKPF_CSA_X_SHIFT	0


/* Every hardware is defined in two ways - 1. The bitwise definition
   of the register (reg) and 2. the whole register as a word (wd).
*/

/* PF Copy destination address reg - RW *
 * gg values : 00 - image plngrp
 *	       01 - depth plngrp
 *	       10 - window plngrp
 *	       11 - image &depth in parallel
 */
typedef struct {
    unsigned       :8;
    unsigned        gg:2;	       /* copy destination plngrp */
    unsigned        start_y:11;
    unsigned        start_x:11;
} H_fbi_cp_dst_reg_bit;

typedef union {
    H_fbi_cp_dst_reg_bit    reg;
    u_int	    wd;
} H_fbi_cp_dst_reg;


/* PF fill destination address reg - RW *
 * gg values : 00 - image plngrp
 *	       01 - depth plngrp
 *	       10 - window plngrp
 *	       11 - ILLEGAL
 */
typedef struct {
    unsigned        :8;
    unsigned        gg:2;	       /* fill destination plngrp */
    unsigned        start_y:11;
    unsigned        start_x:11;
} H_fbi_fl_dst_reg_bit;

typedef union {
    H_fbi_fl_dst_reg_bit    reg;
    u_int	    wd;
} H_fbi_fl_dst_reg;


/* PP View Clip X bounds reg - RW  - specifies VP  LtoR */
typedef struct {
    unsigned        :5;		       
    unsigned        l_bound:11;
    unsigned        :5;		       
    unsigned        r_bound:11;
} H_fbi_vwclp_x_reg_bit;

typedef union {
    H_fbi_vwclp_x_reg_bit   reg;
    u_int	    wd;
} H_fbi_vwclp_x_reg;


/* PP View Clip Y bounds reg - RW  - specifies VP top to bot */
typedef struct {
    unsigned        :5;
    unsigned        t_bound:11;
    unsigned        :5;
    unsigned        b_bound:11;
} H_fbi_vwclp_y_reg_bit;

typedef union {
    H_fbi_vwclp_y_reg_bit   reg;
    u_int	    wd;
} H_fbi_vwclp_y_reg;


/* PP Pick aperture X bounds reg - RW */
typedef struct {
    unsigned        r_bound:11;
    unsigned        :5;
    unsigned        l_bound:11;
    unsigned        :5;
} H_fbi_pick_x_reg_bit;

typedef union {
    H_fbi_pick_x_reg_bit    reg;
    u_int	    wd;
} H_fbi_pick_x_reg;


/* PP Pick aperture Y bounds reg - RW */
typedef struct {
    unsigned        b_bound:11;
    unsigned        :5;		       
    unsigned        t_bound:11;
    unsigned        :5;		       
} H_fbi_pick_y_reg_bit;

typedef union {
    H_fbi_pick_y_reg_bit    reg;
    u_int	    wd;
} H_fbi_pick_y_reg;


/* PP Pick aperture front bounds reg - RW */
typedef struct {
    unsigned		    f_bound:24;
    unsigned		    :8;		       
} H_fbi_pick_f_reg_bit;

typedef union {
    H_fbi_pick_f_reg_bit    reg;
    u_int	    wd;
} H_fbi_pick_f_reg;


/* PP Pick aperture back bounds reg - RW */
typedef struct {
    unsigned		    b_bound:24;
    unsigned		    :8;		       
}               H_fbi_pick_b_reg_bit;

typedef union {
    H_fbi_pick_b_reg_bit    reg;
    u_int	    wd;
} H_fbi_pick_b_reg;

/* PP Pick control reg - RW */
typedef struct {
    unsigned		    ap:1;
    unsigned		    ren_pic:1;
    unsigned		    pic_on:1;
    unsigned		    :29;		       
} H_fbi_pick_ctrl_reg_bit;

typedef union {
    H_fbi_pick_ctrl_reg_bit reg;
    u_int	    wd;
} H_fbi_pick_ctrl_reg;

/* PP FB width reg - RW */
typedef struct {
    unsigned		    :29;		       
    unsigned		    width:2;
} H_fbi_fb_width_reg_bit;

typedef union {
    H_fbi_fb_width_reg_bit  reg;
    u_int	    wd;
} H_fbi_fb_width_reg;


/* PP Stereo Control reg - RW */
typedef struct {
    unsigned		    :29;		       
    unsigned		    st_on:2;	/* Duplicate write enable:0-mono
				        * 1-stereo */
} H_fbi_stereo_reg_bit;

typedef union {
    H_fbi_stereo_reg_bit    reg;
    u_int	    wd;
} H_fbi_stereo_reg;

/* PP buffer select reg - RW */
typedef struct {
    unsigned		    :26;		       
    unsigned		    ao_read_b:1;       /* alpha/overlay: read buf B */
    unsigned		    ao_write_b:1;      /* write buffer B */
    unsigned		    ao_write_a:1;      /* write buffer A */
    unsigned		    rgb_read_b:1;      /* RGB:read buffer B */
    unsigned		    rgb_write_b:1;     /* write buffer B */
    unsigned		    rgb_write_a:1;      /* write buffer A */
} H_fbi_buf_sel_reg_bit;

typedef union {
    H_fbi_buf_sel_reg_bit   reg;
    u_int	    wd;
} H_fbi_buf_sel_reg;


/*  PP Fast Clear Control reg	RW */
typedef struct {
    unsigned		    fcc_ctrl:3;
    unsigned		    :29;		       
} H_fbi_fcc_ctrl_reg_bit;

typedef union {
    H_fbi_fcc_ctrl_reg_bit  reg;
    u_int	    wd;
} H_fbi_fcc_ctrl_reg;


/* PP current window id reg - RW */
typedef struct {
    unsigned		    :22;		       
    unsigned		    cwid:10;
} H_fbi_cur_wid_reg_bit;

typedef union {
    H_fbi_cur_wid_reg_bit  reg;
    u_int		wd;
} H_fbi_cur_wid_reg;

/* PP current WID control reg - RW */
typedef struct {
    unsigned			wid_repl:1;	 /* enable WID replace mode */
    unsigned			:21;
    unsigned			wid_m:10;	 /* WID mask bits */
} H_fbi_wid_mask_reg_bit;

typedef union {
    H_fbi_wid_mask_reg_bit	reg;
    u_int		wd;
} H_fbi_wid_mask_reg;


/* PP Constant Z source reg */
typedef struct {
    unsigned			z_src:24;
    unsigned			:8;		
} H_fbi_con_z_reg_bit;

typedef union {
    H_fbi_con_z_reg_bit		reg;
    u_int		wd;
} H_fbi_con_z_reg;

/* PP Z Control  reg - see documentation for this*/
typedef struct {
    unsigned			x_bit:1;
    unsigned			c_bit:1;
    unsigned			w_bit:1;
    unsigned			h_bit:1;
    unsigned			:28;		       
}               H_fbi_z_ctrl_reg_bit;

typedef union {
    H_fbi_z_ctrl_reg_bit	reg;
    u_int		wd;
} H_fbi_z_ctrl_reg;


/* PP Image write mask reg - RW */
typedef struct {
    unsigned			alpha_over:8;
    unsigned			blue:8;
    unsigned			green:8;
    unsigned			red:8;
} H_fbi_im_wmask_reg_bit;

typedef union {
    H_fbi_im_wmask_reg_bit	reg;
    u_int		wd;
} H_fbi_im_wmask_reg;

/* PP window write mask reg - RW */
typedef struct {
    unsigned			:19;		       /* unused */
    unsigned			fast_clr:1;	       /* fast clear planes */
    unsigned			crs_data:1;	       /* cursor data */
    unsigned			crs_ena:1;	       /* cursor enable */
    unsigned			wid:10;	       /* WID planes */
} H_fbi_win_wmask_reg_bit;

typedef union {
    H_fbi_win_wmask_reg_bit	reg;
    u_int		wd;
} H_fbi_win_wmask_reg;

/* PP byte mode channel select reg - RW .
 * 0-alpha/overlay 1-blue 2-green 3-red
 */
typedef struct {
    unsigned			:20;
    unsigned			bmode:2;
} H_fbi_byte_mode_reg_bit;

typedef union {
    H_fbi_byte_mode_reg_bit	reg;
    u_int		wd;
} H_fbi_byte_mode_reg;

/* PP byte mode write mask reg - RW .
 * 0-alpha/overlay 1-blue 2-green 3-red
 */
typedef struct {
    unsigned			:24;
    unsigned			bmask:8;
} H_fbi_byte_wmask_reg_bit;

typedef union {
    H_fbi_byte_wmask_reg_bit	reg;
    u_int		wd;
} H_fbi_byte_wmask_reg;

/* PP ROP/Blend reg - RW.
 * Bit 4 must be 0 for state set 0.
 */
typedef struct {
    unsigned			:27;
    unsigned			bl_or_rop:1;
    unsigned			brop:4;
} H_fbi_rop_reg_bit;

typedef union {
    H_fbi_rop_reg_bit		reg;
    u_int		wd;
} H_fbi_rop_reg;

typedef struct {
    unsigned			:24;
    unsigned			mpg_set:8;
} H_fbi_mpg_reg_bit;

typedef union {
    H_fbi_mpg_reg_bit		reg;
    u_int		wd;
} H_fbi_mpg_reg;

/* PP Stencil mask RW */
typedef struct {
    unsigned			mask:32;
} H_fbi_sten_mask_reg_bit;

typedef union {
    H_fbi_sten_mask_reg_bit	reg;
    u_int		wd;
} H_fbi_sten_mask_reg;

/* PF Foreground color - RW */
/* Image planes */
typedef struct {
    unsigned			alpha_over:8;
    unsigned			blue:8;
    unsigned			green:8;
    unsigned			red:8;
} H_fbi_im_sten_fg_reg_bit;

typedef union {
    H_fbi_im_sten_fg_reg_bit	reg;
    u_int		wd;
} H_fbi_im_sten_fg_reg;

/* depth planes */
typedef struct {
    unsigned			:8;
    unsigned			depth:24;
} H_fbi_dep_sten_fg_reg_bit;

typedef union {
    H_fbi_dep_sten_fg_reg_bit	reg;
    u_int		wd;
} H_fbi_dep_sten_fg_reg;

/* window id planes */
typedef struct {
    unsigned			:12;
    unsigned			fc:4;
    unsigned			:4;
    unsigned			curs_data:1;
    unsigned			cur_enab:1;
    unsigned			wid_pl:10;
} H_fbi_wid_sten_fg_reg_bit;

typedef union {
    H_fbi_wid_sten_fg_reg_bit	reg;
    u_int		wd;
} H_fbi_wid_sten_fg_reg;

typedef union {
    H_fbi_im_sten_fg_reg_bit	fg_im;
    H_fbi_dep_sten_fg_reg_bit	fg_dep;
    H_fbi_wid_sten_fg_reg_bit	fg_wid;
    u_int		wd;
} H_fbi_sten_fg_reg;

/* PF Stencil background color */
/* Image planes */
typedef struct {
    unsigned			alpha_over:8;
    unsigned			blue:8;
    unsigned			green:8;
    unsigned			red:8;
} H_fbi_im_sten_bg_reg_bit;

typedef union {
    H_fbi_im_sten_bg_reg_bit	reg;
    u_int		wd;
} H_fbi_im_sten_bg_reg;

/* depth planes */
typedef struct {
    unsigned			:8;
    unsigned			depth:24;
} H_fbi_dep_sten_bg_reg_bit;

typedef union {
    H_fbi_dep_sten_bg_reg_bit	reg;
    u_int		wd;
} H_fbi_dep_sten_bg_reg;

/* window id planes */
typedef struct {
    unsigned			:12;
    unsigned			fc:4;
    unsigned			:4;
    unsigned			curs_data:1;
    unsigned			cur_enab:1;
    unsigned			wid_pl:10;
} H_fbi_wid_sten_bg_reg_bit;

typedef union {
    H_fbi_wid_sten_bg_reg_bit	reg;
    u_int		wd;
} H_fbi_wid_sten_bg_reg;

typedef union {
    H_fbi_im_sten_bg_reg_bit	bg_im;
    H_fbi_dep_sten_bg_reg_bit	bg_dep;
    H_fbi_wid_sten_bg_reg_bit	bg_wid;
    u_int			wd
} H_fbi_sten_bg_reg;



/* PF Stencil Transparency flag - RW */
typedef struct {
    unsigned			:31;		
    unsigned			stransp:1;       /* stencil transparency */
} H_fbi_sten_transp_reg_bit;

typedef union {
    H_fbi_sten_transp_reg_bit	reg;
    u_int		wd;
} H_fbi_sten_transp_reg;

/* PF Copy source address reg - RW *
 * gg values : 00 - image plngrp
 *	       01 - depth plngrp
 *	       10 - window plngrp
 *	       11 - image &depth in parallel
 */
typedef struct {    
    unsigned			:8;			
    unsigned			gg:2;	       /* copy src plngrp */
    unsigned			start_y:11;
    unsigned			start_x:11;
} H_fbi_cp_src_reg_bit;

typedef union {
    H_fbi_cp_src_reg_bit	reg;
    u_int		wd;
} H_fbi_cp_src_reg;


/* PF Copy/fill dir/size reg - RW *
 * D- copy direction
 * 0 - outer loop top to bot, inner loop L to R; start@upper L
 * 1 - outer loop bot to top, inner loop R to L; start@lower R
 * D=0 for fills.
 */
typedef struct {
    unsigned			dir:1;	       /* fill dir */
    unsigned			:9;		       
    unsigned			hgt:11;
    unsigned			wdt:11;
} H_fbi_dir_sz_reg_bit;

typedef union {
    H_fbi_dir_sz_reg_bit	reg;
    u_int		wd;
} H_fbi_dir_sz_reg;

/* FBI registers defined as a structure - from a pixrect view point	*/

typedef struct {
    H_fbi_cp_dst_reg	    copy_dst;
    H_fbi_fl_dst_reg	    fill_dst;  
} H_fbi_go;


typedef struct {
#define HPAD0	0x8c
    u_char		    pad0[HPAD0];
    H_fbi_vwclp_x_reg	    vwclp_x;
#define HPAD1	0x4
    u_char		    pad1[HPAD1];
    H_fbi_vwclp_y_reg	    vwclp_y;
#define HPAD2	0x30
    u_char		    pad2[HPAD2];
    H_fbi_fb_width_reg	    fb_wid;
#define HPAD3	0x138
    u_char		    pad3[HPAD3];
    H_fbi_buf_sel_reg	    buf_sel;
    H_fbi_stereo_reg	    stereo;
#define HPAD4	0x1f4
    u_char		    pad4[HPAD4];
    H_fbi_cur_wid_reg	    cur_wid;
    H_fbi_wid_mask_reg	    wid_mask;
#define HPAD5	0x4
    u_char		    pad5[HPAD5];
    H_fbi_con_z_reg	    con_z;
    H_fbi_z_ctrl_reg	    z_ctrl;
    H_fbi_im_wmask_reg	    i_wmask;
    H_fbi_win_wmask_reg	    w_wmask;
    H_fbi_byte_mode_reg	    b_mode;
    H_fbi_byte_wmask_reg    b_wmask;
    H_fbi_rop_reg	    rop;
#define HPAD6	0x0c
    u_char		    pad6[HPAD6];
    H_fbi_mpg_reg	    mpg_set;
#define HPAD7	0x7c8
    u_char		    pad7[HPAD7];
    H_fbi_sten_mask_reg	    sten_mask;
    H_fbi_sten_fg_reg	    fg_col;
    H_fbi_sten_bg_reg	    bg_col;
    H_fbi_sten_transp_reg   transp_flag;
    H_fbi_dir_sz_reg	    dir_size;
    H_fbi_cp_src_reg	    copy_src;


} H_fbi;


#endif	/* !gt_fbi_DEFINED */
