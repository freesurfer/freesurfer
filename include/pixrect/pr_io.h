/*      @(#)pr_io.h 1.7 88/02/08 SMI      */

/*
 * Copyright (c) 1984 by Sun Microsystems, Inc.
 */

#ifndef pr_io_DEFINED
#define pr_io_DEFINED

/*
 * Data definitions and constant definitions for clients of utilities
 *   to do i/o between pixrects and rasterfiles.
 *
 */

typedef struct {
	int		 type;
	int		 length;
	unsigned char 	*map[3];
} colormap_t;
/*
 * The type should be one of the supported ras_maptype's in rasterfile.h.
 * Note: the interpretation of the length and map array depends on the type.
 */


/* error messages */
#define PR_IO_ERR_BADDEPTH	"Cannot replicate color pixrect if depth != 8"
#define PR_IO_ERR_CDEPTH	\
			"Color raster depth does not match display depth"
#define PR_IO_ERR_DISPLAY	"Cannot open display device"
#define PR_IO_ERR_INFILE	"Cannot open input file"
#define PR_IO_ERR_NOBWDISPLAY	\
			"Cannot open display device for monochrome raster"
#define PR_IO_ERR_NOCDISPLAY	"Cannot open display device for color raster"
#define PR_IO_ERR_NONULL	"Cannot process null argument"
#define PR_IO_ERR_OUTFILE	"Cannot open output file"
#define PR_IO_ERR_PIXRECT	"Internal pixrect operation failed"
#define PR_IO_ERR_RASREAD	"Error reading rasterfile input data"
#define PR_IO_ERR_RASREADH	"Error reading rasterfile header"
#define PR_IO_ERR_RASWRITE	"Error writing rasterfile output data"
#define PR_IO_ERR_RASWRITEH	"Error writing rasterfile header"

/* declarations */
extern int pr_dump();
extern int pr_dump_header();
extern int pr_dump_image();
extern Pixrect * pr_dump_init();
extern Pixrect * pr_load();
extern int pr_load_colormap();
extern int pr_load_header();
extern Pixrect * pr_load_image();
extern Pixrect * pr_load_std_image();

#endif pr_io_DEFINED
