/* @(#)gtvar.h 1.9 91/08/16 SMI */

/* Copyright (c) 1989 by Sun Microsystems, Inc.  */

/*
 *	GT Header File
 *	Describe "Private Data" of GT Pixrect,  also external
 *	references to pixrect operation, and some macros.
 */
#ifndef	gtvar_DEFINED
#define	gtvar_DEFINED
#include <sys/param.h>
#include <pixrect/memvar.h>
#include <pixrect/pr_impl_util.h>
#include <sbusdev/gtreg.h>
#include <sbusdev/gtmcb.h>
#include <pixrect/gt_rp.h>
#include <pixrect/gt_fbi.h>
#include <pixrect/pr_planegroups.h>
#include <sun/fbio.h>

#define TRUE			1
#define	GT_NFBS			6	/* number of frame buffers in hawk */
#define GT_DEPTH_24		32
#define GT_DEPTH_WID		32
#define GT_DEPTH_8		8
#define GT_DEPTH_1		1
#define GT_OVERLAY_SHIFT24	24
#define GT_IMAGE_LINEBYTES	8192
#define GT_OVERLAY_LINEBYTES	2048
#define GT_BYTES_PER_WORD	4
#define GT_SRC1_PIX_PER_WORD	32
#define GT_SRC1_PIX_PER_BYTE	8
#define GT_STEN_FAC		4      /* In sten mode addr incr= 4*(#pixels)*/
#define GT_STEN_INCR		128    /* addr incr for u_long =4*32	*/
#define GT_COPY_HACK		15
#define GT_DBL_BUF_OFF		0x0
#define GT_DBL_BUF_ON		0x1
#define GT_DBL_FORE		2       /* Just like PW_DBL_FORE in pw_dbl.h */
#define GT_DBL_BACK		3
#define GT_DBL_BOTH		4
#define GT_DBL_NONE		5
#define GT_DISP_STEREO		0x1	/* Display/monitor is stereo */


 
#define FALSE			0
/* description of single GT frame buffer */
struct gt_fb {
    int             group;
    int             depth;
    struct mprp_data mprp;  
    int             state;	       
};

/* GT  Pixrect private data */
struct gt_data {
    struct mprp_data mprp;		    /* memory pixrect simulator */
    int			    flag;	    /* multi purpose flag */
    int			    planes;	    /* current group and mask */
    int			    fd;		    /* file descriptor */
    int			    buf;	    /* selected buffer */
    struct gt_fb	    fb[GT_NFBS];    /* frame buffer info */
    short		    active;	    /* active fb no. */
    short		    wid_part;	    /* Partition of wids bet I/O */
    struct fb_wid_dbl_info  wid8;	    /* for 8bit windows */
    struct fb_wid_dbl_info  wid24;	    /* for 24bit windows */
    int			    clut_id;	    /* CLUT attached to pr */
    int			    mpg_fcs;	    /* FCS and enabling writes */
    int			    windowfd;	    /* Inherit from egret - for SV !*/
    H_rp		    *gt_rp;	    /* RP control space */
    H_fbi_go		    *gt_fbi_go;	    /* FBI "go" regs" */
    H_fbi		    *gt_fbi;	    /* FB input section */
    H_fe_hold		    *gt_fe_h;	    /* FE hold register section */
    short		    gt_fe_frozen;   /* Flag indicating FE is dead */
    short		    gt_fe_timeout;  /* Get from environment var -
					     * Default 15 secs - */
    short		    gt_stereo;	    /* default mono = 0 */
};

int		gt_putcolormap();
int		gt_putattributes();
int		gt_rop();
int		gt_ioctl();

#ifndef KERNEL
extern struct pixrectops gtp_ops;

Pixrect		*gt_make();
int		gt_destroy();
Pixrect		*gt_region();
int		gt_getcolormap();
int		gt_getattributes();
int		gt_vector();
int		gt_get();
int		gt_put();
int		gt_stencil();
int		gt_batchrop();
#endif	!KERNEL


/* Hawk Macros */
int gt_dummy_extern;

#define gt_d(pr)	((struct gt_data *) ((pr)->pr_data))

#define gt_waitidle(rp, timeout)					\
	  {								\
	   int spin;							\
           int cnt = 0;							\
	   while ((rp)->csr_reg.wd & HKPBM_HCS_PFB) {			\
		/* spin a bit, give other bus accesses some space */    \
                for (spin=0; spin < cnt; spin++) {                      \
                    /*                                                  \
                     * Some busy work using an external variable to     \
                     * ensure optimizers will not take out loop.        \
                     */                                                 \
                    gt_dummy_extern = spin;                             \
                }                                                       \
                                                                        \
                if (cnt++ > 2000000) {                                  \
                    /*                                                  \
                     * Only get here if major problem with GT system.   \
                     * Slow the loop, reduce CPU usage.                 \
                     */							\
		    if (cnt > 2000 + timeout) {				\
                        /* very dead, abort */                          \
			while ((rp)->csr_reg.wd & HKPBM_HCS_PFB);	\
		    }							\
		}							\
	   }								\
	} 

#define gt_set_as_reg(rp, as) \
	(rp)->as_reg.wd = (as);

#define gt_get_as_reg(rp, as) \
	(as) = (rp)->as_reg.wd;

#define gt_set_buf_reg(fbi, buf) \
	(fbi)->buf_sel.wd = (buf);

#define gt_get_buf_reg(fbi, buf) \
	(buf) = (fbi)->buf_sel.wd;

#define gt_set_sten_mask(fbi, mask) \
	(fbi)->sten_mask.wd = (mask);

#define gt_set_vwclip(fbi, xmin, xmax) \
	(fbi)->vwclp_x.wd = (HKPP_VCX_MASK & ((xmin) << HKPP_VCX_LEFT_SHIFT) | (xmax));

#define gt_set_z_ctrl(fbi) \
	(fbi)->z_ctrl.wd = 0;

#define gt_set_fcs_copy(fbi) \
	(fbi)->w_wmask.wd |= HKPP_WWM_FC;
	
#define gt_set_sten_mode_col(rp, fbi, mo, col) \
	(fbi)->fg_col.wd = (col);		\
	(rp)->as_reg.wd = (mo);

#define gt_set_common_reg(dp, col, mo, orp, b_sel, w_mask) \
      {						\
    	struct gt_data *mgtd = gt_d((dp));		\
	H_rp *mrp = (H_rp *)mgtd->gt_rp;		\
	H_fbi *mfbi = (H_fbi *)mgtd->gt_fbi;		\
							\
	gt_waitidle(mrp, mgtd->gt_fe_timeout);		\
	mrp->as_reg.wd = (mo);				\
	mfbi->fg_col.wd = (col);			\
	mfbi->bg_col.wd = (0);				\
	mfbi->buf_sel.wd = (HKPP_BS_MASK &(b_sel));	\
	mfbi->rop.wd = (orp);				\
	mfbi->mpg_set.wd = (HKPP_MPG_MASK & mgtd->mpg_fcs);	\
							\
	switch (mo) {					\
	    case HKPBM_HFBAS_WINDOW :			\
		mfbi->w_wmask.wd = (w_mask);		\
 	    break;					\
	    case HKPBM_HFBAS_IMAGE:			\
            case HKPBM_HFBAS_ST_IMAGE:			\
		mfbi->i_wmask.wd = (w_mask);		\
	    break;					\
	    case HKPBM_HFBAS_CURSOR_BYTE:		\
		mfbi->b_wmask.wd = (w_mask);		\
		mfbi->i_wmask.wd = (w_mask)<<GT_OVERLAY_SHIFT24;\
	    break;					\
	    case HKPBM_HFBAS_DEPTH  :			\
		mfbi->i_wmask.wd = (w_mask);		\
		mfbi->z_ctrl.wd = HKPP_ZC_WRITE_ENA;	\
	    break;					\
	 }						\
    }

#define gt_acc_fe_hold(fh)					\
	(fh)->fe_hold_reg.wd |= HKFE_HOLD_REQ;

#ifndef KERNEL
#define gt_wait_fe_hold(fh, frozen_flag, timeout)		    \
	  {							    \
	    int spin;						    \
	    int cnt = 0;					    \
	    if ((fh)->fe_hold_reg.wd & HKFE_HOLD_ACK) {		    \
		frozen_flag = FALSE;				    \
	    } else {						    \
		while (!frozen_flag &&				    \
		    (!((fh)->fe_hold_reg.wd & HKFE_HOLD_ACK))) {    \
		    /* spin a bit, give bus accesses some space */  \
		    for (spin = 0; spin < cnt; spin ++) {	    \
			/*                                          \
			 * Some busy work using an external	    \
			 * variable to ensure optimizers will not   \
			 * take out loop.			    \
			*/					    \
			gt_dummy_extern = spin;			    \
		    }						    \
		    if (cnt++ > 2000) {				    \
			/*                                          \
			 * Get here if FE is having some problems.  \
			 * Slow the loop, reduce CPU usage.         \
			*/					    \
			sleep (1);				    \
			if (cnt > 2000 + timeout) {		    \
			   /* front end kaput */		    \
			   frozen_flag = TRUE;			    \
			   printf("GT Hardware timeout (front end)  \
				   from gt_wait_fe_hold\n");	    \
			}					    \
		    }						    \
		}						    \
	    } 							    \
	}
#else
#define gt_wait_fe_hold(fh, frozen_flag, timeout)		    \
	  {							    \
	    int spin;						    \
	    int cnt = 0;					    \
	    if ((fh)->fe_hold_reg.wd & HKFE_HOLD_ACK) {		    \
		frozen_flag = FALSE;				    \
	    } else {						    \
		while (!frozen_flag &&				    \
		    (!((fh)->fe_hold_reg.wd & HKFE_HOLD_ACK))) {    \
		    /* spin a bit, give bus accesses some space */  \
		    for (spin = 0; spin < cnt; spin ++) {	    \
			/*                                          \
			 * Some busy work using an external	    \
			 * variable to ensure optimizers will not   \
			 * take out loop.			    \
			*/					    \
			gt_dummy_extern = spin;			    \
		    }						    \
		    if (cnt++ > 2500) {				    \
			   frozen_flag = TRUE;			    \
		    }						    \
		}						    \
	    } 							    \
	}
#endif KERNEL

#define gt_rel_fe_hold(fh)					\
	(fh)->fe_hold_reg.wd &= ~HKFE_HOLD_REQ;

#define	gt_set_fill_regs(rp, fbi, fbi_g, fe_hold, wid, hgt, x0, y0, gg, frozen_flag, max_timeout)							 \
      {								 \
									 \
	(rp)->csr_reg.wd = HKPBM_HCS_HSRP;				 \
	gt_acc_fe_hold(fe_hold);					 \
	gt_dummy_extern = fbi->i_wmask.wd;						 \
	(fbi)->dir_size.wd = ((HKPF_CFDS_COPY_WIDTH_MASK & (wid)) |	 \
	(HKPF_CFDS_COPY_HEIGHT_MASK &					 \
		    ((hgt) << HKPF_CFDS_COPY_HEIGHT_SHIFT)));		 \
	gt_wait_fe_hold(fe_hold, frozen_flag, max_timeout);		 \
	(fbi_g)->fill_dst.wd =						 \
		    ((HKPF_FDA_Y_MASK &	((y0) << HKPF_FDA_Y_SHIFT))	 \
		    | (HKPF_FDA_X_MASK & (x0))| (gg));			 \
	(rp)->csr_reg.wd = 0;						 \
    }

#define gt_set_copy_regs(rp, fbi, fbi_g, fe_hold, xs, ys, gg, wid, hgt,	dir, xd, yd, frozen_flag, max_timeout)						 \
       {								 \
									 \
	(rp)->csr_reg.wd = HKPBM_HCS_HSRP;				 \
	gt_acc_fe_hold(fe_hold);					 \
	gt_dummy_extern = fbi->i_wmask.wd;				 \
	(fbi)->copy_src.wd  = ((HKPF_CSA_X_MASK & (xs)) |		 \
		(HKPF_CSA_Y_MASK & ((ys) << HKPF_CSA_Y_SHIFT)) | (gg));	 \
	(fbi)->dir_size.wd = ((HKPF_CFDS_COPY_WIDTH_MASK & (wid)) |	 \
		(HKPF_CFDS_COPY_HEIGHT_MASK &				 \
		((hgt) << HKPF_CFDS_COPY_HEIGHT_SHIFT)) |		 \
		((dir) ? HKPF_CFDS_COPY_LEFT:0 ));			 \
	gt_wait_fe_hold(fe_hold, frozen_flag, max_timeout);		 \
	(fbi_g)->copy_dst.wd = ((HKPF_CDA_X_MASK & (xd)) |		 \
		(HKPF_CDA_Y_MASK & ((yd) << HKPF_CDA_Y_SHIFT)) | (gg));	 \
	(rp)->csr_reg.wd = 0;						 \
    } 

#define gt_sten_rop_1(gt_a, mem_a, type, off, by_a, s_lb, val)		\
      {								\
	*(gt_a) = (*(mem_a) << (off));					\
	(gt_a)  = (type *) PTR_ADD((gt_a), GT_IMAGE_LINEBYTES);		\
	(mem_a) = (type *)(((int)(by_a) += (s_lb)) & (val));		\
    }

#define gt_sten_rop_2(gt_a, mem_a)					\
      {								\
	*(gt_a) = *(mem_a)++;						\
	(gt_a) = (u_long *)PTR_ADD((gt_a), GT_STEN_INCR);		\
    }
#define gt_set_stereo(fbi, ster)					\
	(fbi)->stereo.wd = (ster);

#define gt_pr_to_mem(gtpr, mempr, op)				    \
    {									\
	struct gt_data *prd = gt_d(gtpr);				\
									\
	mempr = *gtpr;							\
	mempr.pr_ops = &mem_ops;					\
	if (PIX_ATTRGROUP(prd->planes) == PIXPG_8BIT_COLOR) {		\
	    gt_set_common_reg(gtpr, PIXOP_COLOR(op),			\
		    HKPBM_HFBAS_CURSOR_BYTE, HKPP_ROP_SRC,		\
		    prd->buf,(PIX_ALL_PLANES & (prd->planes)));		\
	} else if (PIX_ATTRGROUP(prd->planes) == PIXPG_24BIT_COLOR) {	\
	    gt_set_common_reg(gtpr, PIXOP_COLOR(op),			\
	           HKPBM_HFBAS_IMAGE, HKPP_ROP_SRC,			\
		   prd->buf,(PIX_ALL_PLANES &(prd->planes)));		\
	}								\
    }	    
#endif	gtvar_DEFINED
