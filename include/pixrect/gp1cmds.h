/* @(#)gp1cmds.h 1.28 91/05/31 SMI  */

/*
 * Copyright 1986-1990 by Sun Microsystems, Inc.
 */

#ifndef gp1cmds_DEFINED
#define gp1cmds_DEFINED

/* GP low level command set */

#define GP1_EOCL   (0 << 8)
#define GP1_USE_CONTEXT   (1 << 8)
#define GP1_PR_VEC   (2 << 8)
#define GP1_PR_ROP_NF   (3 << 8)
#define GP1_PR_ROP_FF   (4 << 8)

#define GP1_PR_PGON_SOL   (5 << 8)
#define GP1_SET_ZBUF   (6 << 8)
#define GP1_SET_HIDDEN_SURF  (7 << 8)
#define GP1_SET_MAT_NUM   (8 << 8)
#define GP1_MUL_POINT_FLT_2D  (9 << 8)

#define GP1_MUL_POINT_FLT_3D  (10 << 8)
#define GP1_XF_PGON_FLT_3D  (11 << 8)
#define GP1_XF_PGON_FLT_2D  (12 << 8)

#define GP1_SET_CLIP_LIST  (15 << 8)
#define GP1_SET_FB_NUM   (16 << 8)
#define GP1_SET_VWP_3D   (17 << 8)
#define GP1_SET_VWP_2D   (18 << 8)
#define GP1_SET_ROP   (19 << 8)

#define GP1_SET_CLIP_PLANES  (20 << 8)
#define GP1_MUL_POINT_INT_2D  (21 << 8)
#define GP1_MUL_POINT_INT_3D  (22 << 8)
#define GP1_SET_FB_PLANES  (23 << 8)
#define GP1_SET_MAT_3D   (24 << 8)

#define GP1_SET_MAT_OFFSET  (26 << 8)
#define GP1_SET_COLOR   (28 << 8)
#define GP1_SET_MAT_2D   (29 << 8)

#define GP1_XF_PGON_INT_3D  (30 << 8)
#define GP1_COPY_MAT   (31 << 8)
#define GP1_MUL_MAT_2D   (32 << 8)
#define GP1_MUL_MAT_3D   (33 << 8)
#define GP1_GET_MAT_2D   (34 << 8)

#define GP1_GET_MAT_3D   (35 << 8)
#define GP1_PROC_LINE_FLT_3D  (36 << 8)
#define GP1_PROC_PGON_FLT_3D  (37 << 8)
#define GP1_PROC_PGON_FLT_2D  (38 << 8)

#define GP2_NOP                         (39 << 8)

#define GP1_PR_LINE   (40 << 8)
#define GP1_PR_POLYLINE   (41 << 8)
#define GP1_SET_LINE_TEX  (42 << 8)
#define GP1_SET_LINE_WIDTH  (43 << 8)
#define GP1_CGI_LINE   (44 << 8)

#define GP1_XF_LINE_FLT_2D  (45 << 8)
#define GP1_XF_LINE_FLT_3D  (46 << 8)
#define GP1_XF_LINE_INT_3D  (47 << 8)
#define GP1_PR_PGON_TEX   (48 << 8)

#define GP1_PR_ROP_TEX   (50 << 8)
#define GP1_SET_PGON_TEX_BLK  (51 << 8)
#define GP1_SET_PGON_TEX  (52 << 8)
#define GP1_SET_PGON_TEX_ORG_SCR (53 << 8)
#define GP1_SET_PGON_TEX_ORG_XF_2D (54 << 8)

#define GP1_SET_PGON_TEX_ORG_XF_3D (55 << 8)
#define GP1_XF_LINE_INT_2D  (57 << 8)
#define GP1_XF_PGON_INT_2D  (58 << 8)
#define GP1_PROC_PGON_INT_2D  (59 << 8)

#define GP1_PROC_LINE_FLT_2D  (60 << 8)
#define GP1_PROC_LINE_INT_2D  (61 << 8)
#define GP1_PROC_LINE_INT_3D  (62 << 8)
#define GP1_PROC_PGON_INT_3D  (63 << 8)
#define GP1_SET_PICK_ID   (64 << 8)
#define GP1_SET_PICK_WINDOW  (65 << 8)
#define GP1_GET_PICK   (66 << 8)
#define GP1_SET_PICK   (67 << 8)
#define GP1_CLEAR_PICK   (68 << 8)

#define GP1_XF_POINT_INT_2D  (69 << 8)
#define GP1_XF_POINT_INT_3D  (70 << 8)
#define GP1_XF_POINT_FLT_2D  (71 << 8)
#define GP1_XF_POINT_FLT_3D  (72 << 8)

#define GP1_CLEAR_CONTEXT  (73 << 8)
#define GP1_SET_DEPTH_CUE  (74 << 8)
#define GP1_SET_DEPTH_CUE_COLORS        (75 << 8)
#define GP2_SET_DEPTH_CUE_PARAMETERS    (76 << 8)

/* number of commands defined so far */
#define GP1_NCMDS   76

/*
 * Commands that exist only on the GP2 start here.
 * True except for some minor exceptions.
 * (Leave a little gap for gp1 expansion.)
 */

#define GP2_SET_RGB_COLOR               (80 << 8)
#define GP2_PROC_PGON_FLT_3D_RGB        (81 << 8)
#define GP2_PROC_PGON_INT_3D_RGB        (82 << 8)
#define GP2_XF_PGON_FLT_3D_RGB          (83 << 8)
#define GP2_XF_PGON_INT_3D_RGB          (84 << 8)
#define GP2_SET_LIGHT                   (85 << 8)
#define GP2_SET_LIGHT_MATRIX            (86 << 8)
#define GP2_SET_REFLECTANCE             (87 << 8)
#define GP2_SET_EYE   (88 << 8)

#define GP2_SET_TRANSPARENCY  (89 << 8)
#define GP2_SET_DITHER   (90 << 8)
#define GP2_SET_LIGHT_OPTIONS           (91 << 8)

#define GP2_XF_SHADE_LINE_FLT_3D        (99 << 8)

#define GP2_XF_RECT_INT_2D              (100 << 8)
#define GP2_XF_RECT_FLT_2D              (101 << 8)
#define GP2_SET_TEXT_TYPE               (102 << 8)
#define GP2_SET_TEXT_PRECISION          (103 << 8)
#define GP2_SET_TEXT_PATH               (104 << 8)
#define GP2_SET_FONT                    (105 << 8)
#define GP2_SET_TEXT_EXPANSION          (106 << 8)
#define GP2_SET_TEXT_SPACING            (107 << 8)
#define GP2_SET_TEXT_MAP                (108 << 8)
#define GP2_LOAD_FONT                   (109 << 8)
#define GP2_SET_TEXT_ATTRIBUTES         (110 << 8)
#define GP2_TEXT_INT_2D                 (111 << 8)
#define GP2_TEXT_FLT_2D                 (112 << 8)
#define GP2_TEXT_INT_3D                 (113 << 8)
#define GP2_TEXT_FLT_3D                 (114 << 8)
#define GP2_SET_TEXT_MAT_NUM            (115 << 8)
#define GP2_SET_CMAP_OFFSET             (116 << 8)
#define GP2_XF_TRISTRIP_FLT_3D_RGB      (117 << 8)

/* number of gp2 commands defined so far. */
#define GP2_NCMDS   117

/*
 * Commands that exist only on the GP3 (cgtwelve)
 * start here.
 *
 */

#define GP3_SET_DB_PLANES_RGB  (98 << 8)
#define GP3_SET_FB_PLANES_RGB  (100 << 8)
#define GP3_SET_ANTI_ALIAS  (102 << 8) /* 0x66 */

#define GP3_PR_ROP24_FF                 (118 << 8)
#define GP3_PR_ROP24_NF                 (119 << 8)
#define GP3_PR_ROP24_TEX                (120 << 8)
#define GP3_PR_ROP24_BATCH              (121 << 8)
#define GP3_PR_VEC24                    (122 << 8)

#define GP3_TRISTAR_FLT_3D_RGB  (125 << 8)

/* number of gp3 commands defined so far. */
#define GP3_NCMDS                       125


/* Constants for option field of commands */
#define GP1_FREEBLKS  0x80 /* for EOCL */

#define GP1_NOHIDDENSURF 0 /* for SET_HIDDENSURF */
#define GP1_ZBHIDDENSURF 1 /* depth buffer polygons */
#define GP1_ZBLINES  2 /* depth buffer lines */
#define GP1_ZBMARKERS  4 /* depth buffer markers */
#define GP1_ZBALL  7 /* depth buffer all primitives */

#define GP1_SHADE_CONSTANT 0 /* for XF_POLYGON */
#define GP1_SHADE_GOURAUD 1
#define GP1_SHADE_TEX  2

#define GP1_CLIP_PLANE_LEFT 0x20 /* for SET_CLIP_PLANES */
#define GP1_CLIP_PLANE_RIGHT 0x10
#define GP1_CLIP_PLANE_BOTTOM 0x8
#define GP1_CLIP_PLANE_TOP 0x4
#define GP1_CLIP_PLANE_HITHER 0x2
#define GP1_CLIP_PLANE_YON 0x1
#define GP1_CLIP_PLANES_2D 0x3C
#define GP1_CLIP_PLANES_3D 0x3F

#define GP1_PICK_OFF  0 /* for SET_PICK and CLEAR_PICK */
#define GP1_PICK_DRAW  1
#define GP1_PICK_NODRAW  2

#define GP1_MARKER_LINE  0
#define GP1_MARKER_POLY  1
#define GP1_MARKER_CIRC  2
#define GP1_MARKER_FILL_CIRC 3

#define GP1_MARKER_PHIGS_DOT 4
#define GP1_MARKER_PHIGS_1 4
#define GP1_MARKER_PHIGS_PLUS 5
#define GP1_MARKER_PHIGS_2 5
#define GP1_MARKER_PHIGS_STAR 6
#define GP1_MARKER_PHIGS_3 6
#define GP1_MARKER_PHIGS_CIRCLE 7
#define GP1_MARKER_PHIGS_4 7
#define GP1_MARKER_PHIGS_CROSS 8
#define GP1_MARKER_PHIGS_5 8
#define GP1_MARKER_PHIGS_SQUARE 9
#define GP1_MARKER_PHIGS_6 9
#define GP1_MARKER_PHIGS_BOW_NE 10
#define GP1_MARKER_PHIGS_7 10
#define GP1_MARKER_PHIGS_BOW_NW 11
#define GP1_MARKER_PHIGS_8 11
#define GP1_MARKER_PHIGS(i) (GP1_MARKER_FILL_CIRC + i)


#define GP1_DEPTH_CUE_OFF 0 /* for GP1_SET_DEPTH_CUE */
#define GP1_DEPTH_CUE_ON 1


#define GP2_SET_CMAP_DEFAULT    0
#define GP2_SET_CMAP_ZERO       1

#define GP2_INDEX_COLOR  0   /* For GP2_XF_PGON_XXX_3D_RGB */
#define GP2_RGB_COLOR_TRIPLE    1   /*     GP1_SET_DEPTH_CUE_COLORS */
#define GP2_RGB_COLOR_PACK      2   /* and GP2_SET_RGB_COLOR */
#define GP2_VERTEX_NORMALS 4
#define GP2_POLYGON_NORMAL 8

/* For GP2_SET_RGB_COLOR */
#define GP3_SET_BACK_COLOR 128 /* 0x80 */

/* For GP2_XF_TRISTRIP_FLT_3D_RGB */
/* #define GP2_VERTEX_NORMALS   4 same as above */
#define GP2_FACET_NORMALS 8
#define GP2_STRIP_NORMAL 16
#define GP3_FACET_COLORS 32

/* For GP1_SET_ZBUF */
#define GP2_SET_Z 0     /* a 1 will work the same */
#define GP2_SET_I 2
#define GP2_SET_ZI 3
#define GP3_SET_I_DITH 6     /* will dither if dithering is on */
#define GP3_SET_ZI_DITH 7



/* Light source types or states */

#define GP2_LIGHT_AMBIENT       0 /* 0x00 */
#define GP2_LIGHT_DIRECTIONAL   1 /* 0x01 */
#define GP2_LIGHT_POSITIONAL 2 /* 0x02 */
#define GP2_LIGHT_SPOT      4 /* 0x04 */

#define GP2_FRONT_PROPERTIES 0 /* 0x00 */
#define GP2_BACK_PROPERTIES 1 /* 0x01 */

#define GP2_NO_FACE_REJ  0 /* 0x00 */
#define GP2_BACK_FACE_REJ       1 /* 0x01 */
#define GP2_FRONT_FACE_REJ      2 /* 0x02 */

#define GP2_FACE_REJECTION      0 /* 0x00 */
#define GP2_BACK_FACE_SHADE     1 /* 0x01 */
#define GP2_COPY_PNORM_TO_VNORM 2 /* 0x02 */
#define GP2_GENERATE_PNORMAL    4 /* 0x04 */
#define GP2_LIGHT_ON  8 /* 0x08 */
#define GP3_LIGHT_FRONT  9 /* 0x09 */
#define GP3_LIGHT_BACK  10 /* 0x0a */
#define GP2_LIGHT_OFF  16 /* 0x10 */
#define GP2_IGNORE_FLAG         32 /* 0x20 */

#define GP2_IGNORE_VERTEX_NORMAL 33 /* 0x21 */
#define GP2_IGNORE_VERTEX_COLOR  34 /* 0x22 */
#define GP2_IGNORE_VERTEX_DATA   35 /* 0x23 */
#define GP2_IGNORE_LIGHTING      20 /* 0x14 */
#define GP2_IGNORE_FRONT_FACE    0 /* 0x00 */
#define GP2_IGNORE_BACK_FACE     1 /* 0x01 */

#define GP2_EYE_DIRECTIONAL      0 /* 0x00 */
#define GP2_EYE_POSITIONAL  1 /* 0x01 */

/* For use in conjunction with GP2_SET_LIGHT_OPTIONS
   with the GP2_GENERATE_NORMAL option where you can have
   (ON or OFF) OR'd with (LEFTHAND or RIGHTHAND) to specify
   by which rule the application would like normals generated.
   This option is only available through the GP3.              */

#define GP3_USE_LEFTHAND_RULE 1 /* 0x01 */
#define GP3_USE_RIGHTHAND_RULE 2  /* 0x02 */

/* For use with all XF_PGON commands
 *
 * "or"ing this bit into a polygon command will cause quads to be split into
 * two triangles for increased performance.  Note that this will cause concave
 * quads to be drawn incorrectly.
 */
#define GP2_DIVIDE_QUADS 0x40

/* Constants for text commands */
#define GP2_TEXT_CLEAR_FONT    0
#define GP2_TEXT_FONT0         0
#define GP2_TEXT_FONT1         1
#define GP2_TEXT_FONT2         2
#define GP2_TEXT_FONT3         3
#define GP2_TEXT_FONT4         4
#define GP2_TEXT_FONT5         5
#define GP2_TEXT_FONT6         6
#define GP2_TEXT_FONT7         7
#define GP2_TEXT_RIGHT         0
#define GP2_TEXT_LEFT          1
#define GP2_TEXT_UP            2
#define GP2_TEXT_DOWN          3
#define GP2_TEXT_STRING        0
#define GP2_TEXT_CHAR          1
#define GP2_TEXT_STROKE        2
#define GP2_TEXT_VECTOR        0
#define GP2_TEXT_FILLED        1
#define GP2_TEXT_SPLINE        2

/* Constants for size limits on commands */

/* maximum # of 16 bit texture words for PR_PGON_TEX, PR_ROP_TEX* */
#define GP1_MAXPRTEXSHORTS 2048

/* limits for PR_LINE, PR_POLYLINE */
#define GP1_MAX_BRUSH_WIDTH 181  /* sqrt(32767) */
#define GP1_MAX_PAT_SEGS 16

/* Constant to indicate version information is available
   in 3.2FCS and later releases of the pixrect library */
#define GP1_VERSION_QUERY 1


/* Old command names still available for compatibility */

#define GP1_USEFRAME (1 << 8)
#define GP1_PRVEC (2 << 8)
#define GP1_PRROPNF (3 << 8)
#define GP1_PRROPFF (4 << 8)
#define GP1_PRPOLYSOL (5 << 8)
#define GP1_SETZBUF (6 << 8)
#define GP1_SETHIDDENSURF (7 << 8)
#define GP1_SELECTMATRIX (8 << 8)
#define GP1_XFORMPT_2D (9 << 8)
#define GP1_XFORMPT_3D (10 << 8)
#define GP1_XFPOLYGON_3D (11 << 8)
#define GP1_XFPOLYGON_2D (12 << 8)
#define GP1_CORENDCVEC_3D (13 << 8)
#define GP1_CGIVEC (14 << 8)
#define GP1_SETCLPLST (15 << 8)
#define GP1_SETFBINDX (16 << 8)
#define GP1_SETVWP_3D (17 << 8)
#define GP1_SETVWP_2D (18 << 8)
#define GP1_SETROP (19 << 8)
#define GP1_SETCLIPPLANES (20 << 8)
#define GP1_SETPIXPLANES (23 << 8)
#define GP1_SET_MATRIX_3D (24 << 8)
#define GP1_XFVEC_3D (25 << 8)
#define GP1_XFVEC_2D (27 << 8)
#define GP1_SETCOLOR (28 << 8)
#define GP1_SET_MATRIX_2D (29 << 8)
#define GP1_CORENDCPOLY_3D (30 << 8)
#define GP1_MATMUL_2D (32 << 8)
#define GP1_MATMUL_3D (33 << 8)
#define GP1_GETMATRIX_2D (34 << 8)
#define GP1_GETMATRIX_3D (35 << 8)
#define GP1_COREWLDVECNDC_3D (36 << 8)
#define GP1_COREWLDPOLYNDC_3D (37 << 8)
#define GP1_SET_EF_TEX (39 << 8)

/*
 * Macros for 32 bit accesses to GP shared memory.
 *
 * "p" is a pointer and "a" is a float or int argument.
 */

#define GP1_PUT_S(p, a) (*p++ = a)
#define GP1_GET_S(p, a) (a = *p++)

#ifdef mc68000

#define GP1_PUT_F(p, a) (* (float *) (p) = (a), \
    (p) += sizeof (float) / sizeof *(p))
#define GP1_PUT_I(p, a) (* (int *) (p) = (a), \
    (p) += sizeof (int) / sizeof *(p))
#define GP1_GET_F(p, a) ((a) = * (float *) (p), \
    (p) += sizeof (float) / sizeof *(p))
#define GP1_GET_I(p, a) ((a) = * (int *) (p), \
    (p) += sizeof (int) / sizeof *(p))

#else  mc68000

#define GP1_PUT_F(p, a) (((short *) (p))[0] = ((short *)&(a))[0], \
    ((short *) (p))[1] = ((short *) &(a))[1] , \
    (p) += sizeof (float) / sizeof *(p))

#define GP1_PUT_I(p, a) (((short *) (p))[0] = ((short *)&(a))[0], \
    ((short *) (p))[1] = ((short *) &(a))[1], \
    (p) += sizeof (int) / sizeof *(p))

#define GP1_GET_F(p, a) (((short *) &(a))[0] = ((short *)(p))[0], \
    ((short *) &(a))[1] = ((short *) (p))[1], \
    (p) += sizeof (float) / sizeof *(p))

#define GP1_GET_I(p, a) (((short *) &(a))[0] = ((short *)(p))[0], \
    ((short *) &(a))[1] = ((short *) (p))[1], \
    (p) += sizeof (int) / sizeof *(p))

#endif mc68000


/* The names of these macros have been changed to avoid
 * conflict with the pixwin calls
 */
#define gp1_pw_width(p) ((p)->pw_clipdata->pwcd_prmulti->pr_size.x)
#define gp1_pw_height(p) ((p)->pw_clipdata->pwcd_prmulti->pr_size.y)
#define gp1_pw_offset_X(p) (gp1_d((p)->pw_clipdata->pwcd_prmulti)->cgpr_offset.x)
#define gp1_pw_offset_Y(p) (gp1_d((p)->pw_clipdata->pwcd_prmulti)->cgpr_offset.y)

#endif gp1cmds_DEFINED
