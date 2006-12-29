/* @(#)cg8var.h 1.19 of 7/1/91 SMI */

/*
 * Copyright 1988 by Sun Microsystems, Inc.
 */

#ifndef cg8var_DEFINED
#define    cg8var_DEFINED

#include <sys/ioccom.h>            /* Define _IOR and _IOW */
#include <pixrect/pixrect.h>    /* Definition for struct rect */
#include <pixrect/cg4var.h>
#include <sbusdev/memfb.h>        /* Needed by sbusdev/cg8reg.h. */
#ifdef sun4
#include <sundev/cg8reg.h>
#endif
#include <sunwindow/cms.h>

/* FBIOSATTR device specific array indices, copied from cg4var.h */
#define    FB_ATTR_CG8_SETOWNER_CMD    0    /* 1 indicates PID is valid */
#define    FB_ATTR_CG8_SETOWNER_PID    1    /* new owner of device */


#define    CG8_NFBS    8

#define CG8_PRIMARY         0x01    /* Mark the PRIMARY Pixrect    */
#define CG8_OVERLAY_CMAP    0x02    /* Overlay CMAP to be changed    */
#define CG8_24BIT_CMAP      0x04    /* 24 Bit CMAP to be changed     */
#define CG8_KERNEL_UPDATE   0x08    /* kernel vs. user ioctl    */
/* 0x10 & 0x20 are dbl buf in cg9 */
#define CG8_SLEEPING        0x40    /* Denote if wake_up is necessary */
#define CG8_COLOR_OVERLAY   0x80    /* view overlay & enable as 2 bits */
#define CG8_UPDATE_PENDING  0x100
#define CG8_PIP_PRESENT     0x200   /* PIP is present. */
#define CG8_STOP_PIP        0x400   /* Stop PIP when accessing this pln. grp. */
#define CG8_EIGHT_BIT_PRESENT 0x800 /* There is an 8-bit frame buffer. */

struct cg8_data
{
  struct mprp_data   mprp;          /* memory pixrect simulator */
  int                flags;         /* misc. flags */
  int                planes;        /* current group and mask */
  int                fd;            /* file descriptor */
  short              active;        /* active fb no. */
  int                num_fbs;       /* number of frame buffers in "fb. */
  struct cg4fb       fb[CG8_NFBS];  /* frame buffer info */
  int                windowfd;      /* if 8-bit indexed pw */
  struct colormapseg cms;           /* if 8-bit indexed pr */
  int                real_windowfd; /* if 8-bit indexed pw */
};

#define    cg8_d(pr)    ((struct cg8_data *)((pr)->pr_data))

#define CG8_PR_TO_MEM(src, mem)                     \
    if (src && src->pr_ops != &mem_ops)                 \
    {                                   \
 mem.pr_ops = &mem_ops;                      \
 mem.pr_size = src->pr_size;                 \
 mem.pr_depth = src->pr_depth;                   \
 mem.pr_data = (char *) &cg8_d(src)->mprp;           \
 src = &mem;                         \
    }

#define CG8_PR_TO_MEM32(pr, mem32_pr, mem32_pr_data)            \
    if (pr && pr->pr_ops != &mem_ops)                   \
    {                                   \
 mem32_pr.pr_ops = &mem32_ops;                   \
 mem32_pr.pr_size = pr->pr_size;                 \
 mem32_pr.pr_depth = pr->pr_depth;               \
 mem32_pr_data.mprp = cg8_d(pr)->mprp;               \
 mem32_pr_data.plane_group = cg8_d(pr)->fb[cg8_d(pr)->active].group;\
 mem32_pr_data.fd = cg8_d(pr)->fd;               \
 mem32_pr_data.windowfd = cg8_d(pr)->windowfd;           \
 mem32_pr_data.cms = cg8_d(pr)->cms;             \
 mem32_pr.pr_data = (char *) &mem32_pr_data;         \
 pr = &mem32_pr;                         \
    }


extern struct pixrectops cg8_ops;

int        cg8_putcolormap();
int        cg8_putattributes();
int        cg8_ioctl();

#ifndef KERNEL

Pixrect    *cg8_make();
int        cg8_destroy();
Pixrect    *cg8_region();
int        cg8_getcolormap();
int        cg8_getattributes();
int        cg8_vector();
int        cg8_get();
int        cg8_put();
int        cg8_rop();

#endif    !KERNEL

#ifndef pipio_DEFINED

#define pipio_DEFINED

/*  IOCTL definitions for the SBus True Card PIP
 *
 */

/*    First unused ioctl number is 43.*/

/*    EEPROM Write Byte operation control information:
 */
typedef struct
{
  int        address;        /* Address in EEPROM to write to. */
  int        value;            /* Value to write (in low-order 8 bits.) */
}
EEPROM_Write_Byte;

/*    Frame buffer and memory mapping information ioctl:
 */
typedef struct    Fb_Description
{ /* Frame buffer description: */
  short    group;                /* ... Sun "plane group" of frame buffer. */
  short    width;                /* ... width of frame buffer. */
  short    height;               /* ... height of frame buffer. */
  short    depth;                /* ... depth of frame buffer. */
  u_int    linebytes;            /* ... # of bytes per scan line for fb. */
  u_int    mmap_size;            /* ... size of mapping for frame buffer. */
  u_int    mmap_offset;          /* ... offset for memory map of fb. */
}
Fb_Description;

#define PIP_NFBS 10  /* # of frame buffer descriptions in Pipio_Fb_Info. */
#define FB_NPGS  12  /* # of plane groups possible. */

typedef struct Pipio_Fb_Info
{              /* Frame buffer info record: */
  int            frame_buffer_count;        /* ... # of fbs supported. */
  u_int          total_mmap_size;           /* ... memory map size of all fbs */
  Fb_Description fb_descriptions[PIP_NFBS]; /* ... individual fb descriptions */
}
Pipio_Fb_Info;

/*    Frame buffer emulation ioctl:
 */
typedef struct Pipio_Emulation
{      /* Emulation control layout: */
  u_char  plane_groups[FB_NPGS];    /* ... plane groups to enable. */
  u_short timing;                   /* ... timing/size regimen. */
}
Pipio_Emulation;
#define NATIVE_TIMING 0        /* Provide standard (default) timing. */
#define NTSC_TIMING   1        /* Provide NTSC timing. */
#define PAL_TIMING    2        /* Provide PAL timing. */

/*    I/O controls used by Sunview Pixrect library routines.
 */
#define PIPIO_G_FB_INFO      _IOR(X, 1, Pipio_Fb_Info) /* Get info about fbs. */
#define PIPIO_G_EMULATION_MODE _IOR(X, 3, Pipio_Emulation)
/* Return current emulation mode. */
#define PIPIO_S_EMULATION_MODE _IOW(X, 4, Pipio_Emulation)
/* Set the device being emulated. */
#define PIPIO_G_PIP_ON_OFF     _IOR(X, 5, int)
/* Get the value of the pip on bit. */
#define PIPIO_S_PIP_ON_OFF  _IOW(X, 7, int) /* Set or clear pip on bit. */
#define PIPIO_G_PIP_ON_OFF_RESUME  _IOR(X, 9, int)
/* Resume (pop) pip operations, return new status. */
#define PIPIO_G_PIP_ON_OFF_SUSPEND _IOR(X, 10, int)
/* Get pip status, & suspend pip ops.*/

#define PIPIO_G_CURSOR_COLOR_FREEZE _IOR(X, 40, int)
/* Get setting of cursor color frozen switch. */
#define PIPIO_S_CURSOR_COLOR_FREEZE _IOW(X, 41, int)
/* Set cursor color frozen switch. */
#define PIPIO_S_MAP_SLOT  _IOW(X, 42, int)    /* Map SBus slot at offset 0x900000. */
#define PIPIO_G_TEST      _IOR(X, 43, int)    /* For testing purposes. */
#define PIPIO_S_TEST      _IOW(X, 44, int)    /* For testing purposes. */

#endif pipio_DEFINED

#endif cg8var_DEFINED
