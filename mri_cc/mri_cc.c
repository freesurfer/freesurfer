///////////////////////////////////////////
// mri_cc.c
// 
// written by Peng Yu
// date: 01/27/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2006/09/15 14:43:49 $
// Revision       : $Revision: 1.4 $
////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>

#include "mri.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "mrimorph.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "volume_io/geom_structs.h"
#include "transform.h"
#include "talairachex.h"

//static char vcid[] = "$Id: mri_cc.c,v 1.4 2006/09/15 14:43:49 fischl Exp $";


int             main(int argc, char *argv[]) ; 
static int      get_option(int argc, char *argv[]) ; 
static char     *wmvolume = "mri/wm" ;
char            *Progname ;             
MRI             *mri_wm, *mri_cc_tal ; 
int             dxi;

static Real cc_tal_x = 0.0 ;
static Real cc_tal_y = 0.0 ;
static Real cc_tal_z = 27.0 ;
static LTA *lta = 0;
static int find_cc_slice(MRI *mri, Real *pccx, Real *pccy, Real *pccz, const LTA *lta) ;
static int find_corpus_callosum(MRI *mri, Real *ccx, Real *ccy, Real *ccz, const LTA *lta) ;
static int labels[] =   
{ THICKEN_FILL, NBHD_FILL, VENTRICLE_FILL, DIAGONAL_FILL, DEGENERATE_FILL };
#define NLABELS  sizeof(labels) / (sizeof(labels[0]))
#define MAX_SLICES        15  /* 41*/
#define HALF_SLICES       ((MAX_SLICES-1)/2)
#define CUT_WIDTH         1
#define HALF_CUT          ((CUT_WIDTH-1)/2)
#define SEARCH_STEP       3
#define MAX_OFFSET        50
#define CC_VAL            127

/* aspect ratios are dy/dx */
#define MIN_CC_AREA       350  /* smallest I've seen is 389 */
#define MAX_CC_AREA      1400  /* biggest I've seen is 1154 */
#define MIN_CC_ASPECT     0.1
#define MAX_CC_ASPECT     0.75

double findMinSize(MRI *mri)
{
  double xsize, ysize, zsize, minsize;
  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;
  // there are 3! = 6 ways of ordering
  //             xy  yz  zx
  // x > y > z    z min
  // x > z > y    y min  
  // z > x > y    y min
  //////////////////////////
  // y > x > z    z min
  // y > z > x    x min
  // z > y > x    x min
  if (xsize > ysize)
    minsize = (ysize > zsize) ? zsize : ysize;
  else
    minsize = (zsize > xsize) ? xsize : zsize;

  return minsize;
}



int 
main(int argc, char *argv[]) 
{ 
	char        ifname[STRLEN], ofname[STRLEN],  data_dir[STRLEN], *cp ; 
	int         nargs, msec; 
	int         y, z, xi, yi_low=256, yi_high=0, zi_low=256, zi_high=0, temp;
	int         volume[5], i, j, k;           
	struct timeb  then ;
	MRI         *mri_tal, *mri_talheader;
	Real        xv, yv, zv;
	FILE        *fp;
	MATRIX      *inverse_transform_matrix;
	LT          *lt = 0;
	VOL_GEOM    vgtmp;

	Progname = argv[0] ; 
	DiagInit(NULL, NULL, NULL) ; 
	ErrorInit(NULL, NULL, NULL) ; 

	for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) 
	{ 
		nargs = get_option(argc, argv) ; 
		argc -= nargs ; 
		argv += nargs ; 
	} 
	
	if (argc < 2) 
		ErrorExit(ERROR_BADPARM, 
							"usage: %s <input volume>", Progname); 
	
	TimerStart(&then) ; 

	cp = getenv("SUBJECTS_DIR");
  if (cp==NULL)
  {
	printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
	exit(1);
  }
	strcpy(data_dir, cp) ;
	
	sprintf(ifname,"%s/cc_volume_%d.txt",data_dir,dxi) ;  
	if((fp = fopen(ifname, "a")) == NULL)
  {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "cc volume measurement: file %s does not exist!", ifname));
  }
	print("writing results to %s\n",ifname);

	sprintf(ifname,"%s/%s/%s",data_dir,argv[1],wmvolume) ;  
	fprintf(stderr,"reading white matter volume from %s\n", ifname);	
	mri_wm = MRIread(ifname) ; 
	
	sprintf(ifname,"%s/%s/mri/transforms/talairach.xfm",data_dir,argv[1]) ;  
	lta = LTAreadEx(ifname);
	if (lta==0)
		ErrorExit(ERROR_BADPARM,"ERROR: cound not load lta from %s.\n", ifname);      
	fprintf(stderr, "INFO: Using %s and its offset for Talairach volume ...\n", ifname);
	
  for (i = 0 ; i < NLABELS ; i++)
  {
    MRIreplaceValues(mri_wm, mri_wm, labels[i], 0) ;
  }
	sprintf(ofname,"%s/%s/mri/wmpeng.mgz",cp,argv[1]) ; 
	fprintf(stderr, "writing wm volume to %s...\n", ofname) ; 
	MRIwrite(mri_wm, ofname) ;

	mri_talheader = MRIallocHeader(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_wm, mri_talheader); // not allocate memory, though
 
	ModifyTalairachCRAS(mri_talheader, lta);
 
	mri_tal = MRIalloc(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_talheader, mri_tal);
  // now fill the talairach volume values
  MRItoTalairachEx(mri_wm, mri_tal, lta); 
	
  // binalize the talairach volume (mri_tal)
  MRIbinarize(mri_tal, mri_tal, DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;

	sprintf(ofname,"%s/%s/mri/wm_tal.mgz",cp,argv[1]) ; 
	fprintf(stderr, "writing talairach volume to %s...\n", ofname) ; 
	MRIwrite(mri_tal, ofname) ; 

	mri_cc_tal = MRIcopy(mri_tal, NULL) ;	
	MRIcopyHeader(mri_talheader, mri_cc_tal);
	MRIvalueFill(mri_cc_tal, 0) ;

	find_corpus_callosum(mri_tal,&cc_tal_x,&cc_tal_y,&cc_tal_z, lta);

	xi = nint(cc_tal_x);

	/* Then find the bounding box of CC in the midsaggital plane */

	for (y = 0 ; y < mri_cc_tal->height ; y++) 
	{ 
		for (z = 0 ; z < mri_cc_tal->depth ; z++) 
		{ 	 
			if ( MRIvox(mri_cc_tal, xi, y, z) )
			{
				if (y < yi_low)
					yi_low = y ;
				if (z < zi_low)
					zi_low = z ;
				if (y > yi_high )
					yi_high = y ;
				if (z > zi_high)
					zi_high = z ;
			} 	
		} 
	} 

	volume[0] = 0; volume[1] = 0; volume[2] = 0;
	volume[3] = 0; volume[4] = 0;

	MRIcopy(mri_cc_tal, mri_tal);
	MRIvalueFill(mri_tal, 0) ;	

	for (i = xi-dxi ; i <= xi+dxi ; i++) 
	{		
		for (j = 0 ; j < mri_cc_tal->height ; j++) 
		{ 
			for (k = 0 ; k < mri_cc_tal->depth ; k++) 
			{
				MRIvox(mri_tal, i, j, k) = MRIvox(mri_cc_tal, i, j, k) ;	
			} 	
		}
	} 

#if 0
	MRIfromTalairachEx(mri_tal, mri_wm, lta);
#else
	inverse_transform_matrix = MatrixInverse(lta->xforms[0].m_L,
																					 NULL);
	if(inverse_transform_matrix == NULL)
	{
		fprintf(stderr, "ERROR: inverting transform\n");
		MatrixPrint(stdout,lta->xforms[0].m_L);
		exit(1);
	}
	
	MatrixFree(&(lta->xforms[0].m_L));
	lta->xforms[0].m_L = inverse_transform_matrix;
	// reverse src and dst target info.
	// since it affects the c_ras values of the result
	// in LTAtransform()
	// question is what to do when transform src info is invalid.
	lt = &lta->xforms[0];
	copyVolGeom(&lt->dst, &vgtmp);
	copyVolGeom(&lt->src, &lt->dst);
	copyVolGeom(&vgtmp, &lt->src);
	mri_wm = MRIclone(mri_tal, NULL) ;
	mri_wm->c_r = 0;
	mri_wm->c_a = 0;
	mri_wm->c_s = 0;  
	mri_wm = LTAtransformInterp(mri_tal, mri_wm, lta, SAMPLE_NEAREST);
#endif
	sprintf(ofname,"%s/%s/mri/cc.mgz",cp,argv[1]) ; 
	fprintf(stderr, "writing output to %s...\n", ofname) ; 
	MRIwrite(mri_wm, ofname) ;
 
	//MRItoTalairachEx(mri_wm,mri_tal,lta);
	sprintf(ofname,"%s/%s/mri/cc_tal.mgz",cp,argv[1]) ; 
	fprintf(stderr, "writing output to %s...\n", ofname) ; 
	MRIwrite(mri_tal, ofname) ;

		/* Then count the number of CC in the bounding box in this plane */
	for (i = 0 ; i <= mri_wm->width ; i++) 
	{		
		for (j = 0 ; j < mri_wm->height ; j++) 
		{ 
			for (k = 0 ; k < mri_wm->depth ; k++) 
			{ 	 
				if (MRIvox(mri_wm, i, j, k))
				{
					MRIvoxelToTalairachVoxelEx(mri_wm, i, j, k, &xv, &yv, &zv, lta);					
					if ( yv>=yi_low && yv<=yi_high && zv>=zi_low && zv<=zi_high )
					{
						temp = (int) (zv-zi_low)/((zi_high-zi_low+1)/5);
						volume[temp]++ ;
					}
					else printf ("error in talairach transform\n"); 	
				} 
			}
		}
	}

	fprintf(fp, "%s %d %d %d %d %d %d %d %d %d \n", argv[1], volume[4], volume[3],volume[2], volume[1],volume[0], yi_low, yi_high, zi_low, zi_high);
	fprintf(stderr, "%s %d %d %d %d %d %d %d %d %d \n", argv[1], volume[4], volume[3],volume[2], volume[1],volume[0], yi_low, yi_high, zi_low, zi_high);

	MRIfree(&mri_tal) ;
	MRIfree(&mri_wm); 
	MRIfree(&mri_talheader);
	MRIfree(&mri_cc_tal) ; 
	msec = TimerStop(&then) ; 
	fprintf(stderr, "corpus callosum matter segmentation took %2.1f minutes\n", (float)msec/(1000.0f*60.0f)); 
	fclose(fp);	
	exit(0) ; 
	return(0) ; 
} 

#define CC_SPREAD       10
#define MIN_THICKNESS   3
#define MAX_THICKNESS   20

static int
find_corpus_callosum(MRI *mri_tal, Real *pccx, Real *pccy, Real *pccz, const LTA *lta)
{
  int         xv, yv, zv, max_y, max_thick=0, thickness=0, y1, xcc, ycc, x, y,x0, extension=50 ;
  Real        xr, yr, zr ;
  MRI_REGION  region ;
  int cc_spread, min_thickness, max_thickness, slice_size;
	double voxsize=findMinSize(mri_tal);
  
  cc_spread = ceil(CC_SPREAD/voxsize);
  min_thickness = ceil(MIN_THICKNESS/voxsize); // estimate bigger
  max_thickness = ceil(MAX_THICKNESS/voxsize);
  slice_size = mri_tal->width;

  MRIboundingBox(mri_tal, 1, &region) ;
  // bounding box center position in voxel coords
  x0 = region.x+region.dx/2 ;

  // this function is called with mri being talairached volume
  // get the talairach coords (0,0,0) in the voxel space
  if (mri_tal->linear_transform || lta)
  {
    MRIworldToVoxel(mri_tal, 0.0, 0.0, 0.0, &xr, &yr, &zr);   /* everything is now in tal coords */
    xv = nint(xr) ; yv = nint(yr) ; zv = nint(zr) ;
  }
  else
  {
    xv = x0; yv = region.y+region.dy/2; zv = region.z+region.dz/2; 
  }

	fprintf(stderr, "original seed found at x=%d, y=%d z=%d \n", x0, yv, zv );
  /* find the column with the lowest starting y value of any sign. thick. */
  xcc = ycc = max_y = 0 ; 
  for (x = x0-cc_spread ; x <= x0+cc_spread ; x++)
  {
    /* search for first non-zero pixel */
    // in the talairach origin coronal slice from the top
		while (thickness==0 && yv-extension >= 0 && yv+extension <= 256)
		{
			for (y = yv-extension ; y < yv+extension ; y++)
			{
				if (MRIvox(mri_tal, x, y, zv) >= WM_MIN_VAL)  
					break ;
			}
			// find y which is greater than WM_MIN_VAL
			/* check to make sure it as reasonably thick */
			if (y < yv+extension ) // within the region and bigger than so far
			{
				for (y1 = y, thickness = 0 ; y1 < slice_size ; y1++, thickness++)
					if (!MRIvox(mri_tal, x, y1, zv)) // if becomes zero, then break out
						break ;
				if ( thickness > min_thickness && thickness < max_thickness )
				{
					if ( y > max_y || (y == max_y && thickness > max_thick) )
					{
						// found the zero voxel at y1 -> thinckness
						xcc = x ; ycc = y+thickness/2 ;  /* in middle of cc */
						max_y = y ;            // mark starting y position
						max_thick = thickness;
						if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
							fprintf(stderr, "potential cc found at (%d, %d), thickness = %d\n",
											xcc, ycc, thickness) ;
					}
				}
			}
			extension += 10; 
		}
		thickness = 0; extension = 50;
  }

  if (!max_y)
    return(ERROR_BADPARM) ;

  /* now convert the in-plane coords to Talairach coods */
  MRIvoxelToWorld(mri_tal, xcc, ycc, zv, pccx, pccy, pccz) ;
	fprintf(stderr, "%d, %d, %d\n", xcc, ycc, zv);

  find_cc_slice(mri_tal, pccx, pccy, pccz, lta) ;

  return(NO_ERROR) ;
}


static int
find_cc_slice(MRI *mri_tal, Real *pccx, Real *pccy, Real *pccz, const LTA *lta)
{
  // here we can handle only up to .5 mm voxel size
  int         area[MAX_SLICES*2], min_area, min_slice, slice, offset,xv,yv,zv,
              xo, yo ;
  MRI         *mri_slice, *mri_filled ;
  Real        aspect, x_tal, y_tal, z_tal, x, y, z, xvv, yvv, zvv;
  MRI_REGION  region ;
  char        fname[STRLEN] ;
  int half_slices;
  double voxsize = findMinSize(mri_tal);
  int slice_size = mri_tal->width;
  int max_slices = ceil(MAX_SLICES/voxsize);
  int max_cc_area = ceil(MAX_CC_AREA/(voxsize*voxsize));
  int min_cc_area = floor(MIN_CC_AREA/(voxsize*voxsize));

  half_slices = floor(HALF_SLICES/voxsize);
  if ( half_slices <= 0)
    half_slices = 1;

  x_tal = *pccx ; y_tal = *pccy ; z_tal = *pccz ;
  offset = 0 ;
  xo = yo = (slice_size-1)/2 ;  /* center point of the slice */
  for (slice = 0 ; slice < max_slices ; slice++)
  {
    offset = slice - half_slices ;
    x = x_tal + offset ; y = y_tal ; z = z_tal ;
    MRIworldToVoxel(mri_tal, x, y,  z, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    mri_slice = MRIextractPlane(mri_tal, NULL, MRI_SAGITTAL, xv);
    mri_filled =  MRIfillFG(mri_slice, NULL, zv, yv,0,WM_MIN_VAL,CC_VAL,&area[slice]);
    MRIboundingBox(mri_filled, 1, &region) ;
    aspect = (Real)region.dy / (Real)region.dx ;

    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= slice_size -1 ||
        region.y+region.dy >= slice_size-1)
      area[slice] = 0 ;
    
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "slice[%d] @ (%d, %d, %d): area = %d\n", 
              slice, xv, yv, zv, area[slice]) ;

    if ((Gdiag & DIAG_WRITE) && !(slice % 1) && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "cc_slice%d.mgz", slice);
      MRIwrite(mri_slice, fname) ;
      sprintf(fname, "cc_filled%d.mgz", slice);
      MRIwrite(mri_filled, fname) ;
    }
		MRIfillPlane(mri_filled, mri_cc_tal, MRI_SAGITTAL, xv, CC_VAL);

    MRIfree(&mri_filled) ; MRIfree(&mri_slice) ;
  }
  
  min_area = 10000 ; min_slice = -1 ;
  for (slice = 1 ; slice < max_slices-1 ; slice++)
  {
    if (area[slice] < min_area && 
        (area[slice] >= min_cc_area && area[slice] <= max_cc_area))
    {
      min_area = area[slice] ;
      min_slice = slice ;
    }
  }
  
  /* couldn't find a good slice - don't update estimate */
  if (min_slice < 0)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "%s: could not find valid seed for the cc",
                 Progname));
  
  offset = min_slice - half_slices ;
  *pccx = x = x_tal + offset ; *pccy = y = y_tal ; *pccz = z = z_tal ;
  
  // just for debugging
  MRIworldToVoxel(mri_tal, x, y,  z, &xvv, &yvv, &zvv) ;
	*pccx = xvv ; *pccy = yvv ; *pccz = zvv ;

	fprintf(stderr, "updating initial cc seed to Tal vol (%.2f, %.2f, %.2f) TAL (%.2f, %.2f, %.2f)\n",
	    xvv, yvv, zvv, x, y, z);

  return(NO_ERROR) ;
}




/*----------------------------------------------------------------------
	
	Parameters: 
	
	Description: 
	----------------------------------------------------------------------*/


static int 
get_option(int argc, char *argv[]) 
{ 
	int  nargs = 0 ; 
	char *option ; 
	
	option = argv[1] + 1 ;            /* past '-' */ 
	switch (toupper(*option)) 
	{ 
	case '?': 
	case 'U': 
		fprintf(stderr, 
						"usage: %s <input volumes> <output volume>\n", 
						Progname) ; 
		exit(1) ; 
		break ;
	case 'T':
		dxi = atoi(argv[2]);
		fprintf(stderr,"change thickness to %d mm\n", 2*dxi+1);
		nargs = 1;
		break; 
	default: 
		fprintf(stderr, "unknown option %s\n", argv[1]) ; 
		exit(1) ; 
		break ; 
	} 
	
	return(nargs) ; 
}






