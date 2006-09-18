///////////////////////////////////////////
// mri_cc.c
// 
// written by Peng Yu
// date: 01/27/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2006/09/18 14:22:35 $
// Revision       : $Revision: 1.7 $
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
#include "matrix.h"
#include "mriTransform.h"

//static char vcid[] = "$Id: mri_cc.c,v 1.7 2006/09/18 14:22:35 fischl Exp $";


int             main(int argc, char *argv[]) ; 
static int      get_option(int argc, char *argv[]) ; 
static char     *wmvolume = "mri/wm" ;
char            *Progname ;             
MRI             *mri_wm, *mri_cc_tal ; 
int             dxi=0;
int             x_edge=0, y_edge=0;


static Real cc_tal_x = 0.0 ;
static Real cc_tal_y = 0.0 ;
static Real cc_tal_z = 27.0 ;
static LTA *lta = 0;
static int find_cc_slice(MRI *mri, Real *pccx, Real *pccy, Real *pccz, const LTA *lta, MRI *mri_tal_cc) ;
static int find_corpus_callosum(MRI *mri, Real *ccx, Real *ccy, Real *ccz, const LTA *lta, MRI *mri_tal_cc) ;
static MRI *remove_fornix(MRI *mri_filled, int xv, int yv, int zv);
static int edge_detection(MRI *mri_temp, int edge_count,int signal);
static int labels[] =   
{ THICKEN_FILL, NBHD_FILL, VENTRICLE_FILL, DIAGONAL_FILL, DEGENERATE_FILL };
#define NLABELS  sizeof(labels) / (sizeof(labels[0]))
#define MAX_SLICES        21  /* 41*/
#define HALF_SLICES       ((MAX_SLICES-1)/2)
#define CUT_WIDTH         1
#define HALF_CUT          ((CUT_WIDTH-1)/2)
#define SEARCH_STEP       3
#define MAX_OFFSET        50
#define CC_VAL            100

/* aspect ratios are dy/dx */
#define MIN_CC_AREA       350  /* smallest I've seen is 389 */
#define MAX_CC_AREA      1400  /* biggest I've seen is 1154 */
#define MIN_CC_ASPECT     0.1
#define MAX_CC_ASPECT     0.75

static char sdir[STRLEN] = "" ;

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
	Real        xc,yc,zc;
	MATRIX      *mrot, *mtrans;
	int         i, j, k;           
	struct timeb  then ;
	MRI         *mri_tal, *mri_talheader, *mri_header, *mri_cc;
	Real        xv, yv, zv;
	FILE        *fp;
	LTA         *lta2 = 0;	
	float       volume[5];

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

  if (!strlen(sdir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  strcpy(data_dir, sdir) ;
	
	sprintf(ifname,"%s/cc_volume_%d.txt",data_dir,dxi) ;  
	if((fp = fopen(ifname, "a")) == NULL)
  {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "cc volume measurement: file %s does not exist!", ifname));
  }
	print("writing results to %s\n",ifname);

	sprintf(ifname,"%s/%s/%s",data_dir,argv[1],wmvolume) ;  
	print("reading white matter volume from %s\n", ifname);	
	mri_wm = MRIread(ifname) ; 
	
	sprintf(ifname,"%s/%s/mri/transforms/talairach.xfm",data_dir,argv[1]) ;  
	lta = LTAreadEx(ifname);
	if (lta==0)
		ErrorExit(ERROR_BADPARM,"ERROR: cound not load lta from %s.\n", ifname);      
	fprintf(stdout, "INFO: Using %s and its offset for Talairach volume ...\n", ifname);
	
  for (i = 0 ; i < NLABELS ; i++)
  {
    MRIreplaceValues(mri_wm, mri_wm, labels[i], 0) ;
  }

	mri_talheader = MRIallocHeader(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_wm, mri_talheader); // not allocate memory, though

	ModifyTalairachCRAS(mri_talheader, lta);
 
	mri_tal = MRIalloc(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_talheader, mri_tal);
  // now fill the talairach volume values
  MRItoTalairachEx(mri_wm, mri_tal, lta); 
	
  // binalize the talairach volume (mri_tal)
  MRIbinarize(mri_tal, mri_tal, DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(ofname,"%s/%s/mri/wm_tal.mgz",data_dir,argv[1]) ; 
    fprintf(stdout, "writing talairach transformed white matter volume to %s...\n", ofname) ; 
    MRIwrite(mri_tal, ofname) ; 
  }

	//find the transform matrix
	mtrans = MatrixAlloc(4, 4, MATRIX_REAL) ;
	mrot = MatrixAlloc(4, 4, MATRIX_REAL) ;

	//try method 2 to get the rotation matrix
	sprintf(ifname,"%s/%s/mri/transforms/talairach.xfm",data_dir,argv[1]) ;  
	lta2 = LTAreadEx(ifname);
	mtrans=lta2->xforms[0].m_L;
	Trns_ExtractRotationMatrix (mtrans,mrot);
	*MATRIX_RELT(mrot, 1, 4) = mtrans->rptr[1][4];
	*MATRIX_RELT(mrot, 2, 4) = mtrans->rptr[2][4];
	*MATRIX_RELT(mrot, 3, 4) = mtrans->rptr[3][4];
	lta2->xforms[0].m_L=mrot;

	//rotation wm volume to be upright, using cc volume temporarily
	mri_header = MRIallocHeader(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_wm, mri_header); 
	ModifyTalairachCRAS(mri_header, lta2);
	mri_cc = MRIcopy(mri_wm, NULL) ;	
	MRIcopyHeader(mri_header, mri_cc);
	MRItoTalairachEx(mri_wm, mri_cc, lta2); 
  // binalize the rotated wm  volume (mri_cc)
  MRIbinarize(mri_cc, mri_cc, DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;
  
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(ofname,"%s/%s/mri/wm.mgz",data_dir,argv[1]) ; 
    fprintf(stdout, "writing rotated white matter volume to %s...\n", ofname) ; 
    MRIwrite(mri_cc, ofname) ;
  }

	//now start cc segmentation in talairach space
 	mri_cc_tal = MRIcopy(mri_tal, NULL) ;	
	MRIcopyHeader(mri_talheader, mri_cc_tal);
	MRIvalueFill(mri_cc_tal, 0) ;

	//most of the work is done in find_corpus_callosum function
	find_corpus_callosum(mri_tal,&cc_tal_x,&cc_tal_y,&cc_tal_z, lta, mri_cc_tal);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(ofname,"%s/%s/mri/cc_tal.mgz",data_dir,argv[1]) ; 
    fprintf(stdout, "writing output to %s...\n", ofname) ; 
    MRIwrite(mri_cc_tal, ofname) ;
  }

	//starting the volume measurement
	volume[0] = 0.0; volume[1] = 0.0; volume[2] = 0.0;
	volume[3] = 0.0; volume[4] = 0.0;

	//transform cc volume from talairach space to normal space
	MRIfromTalairachEx(mri_cc_tal, mri_wm, lta);
 // binalize the rotated cc volume (mri_wm)
  MRIbinarize(mri_wm, mri_wm, CC_VAL/2-1, 0, 100) ;
	sprintf(ofname,"%s/%s/mri/cc_org.mgz",data_dir,argv[1]) ; 
	fprintf(stdout, "writing corpus callosum in original space to %s...\n", ofname) ; 
	MRIwrite(mri_wm, ofname) ;

	//trying to find the position of mid-sagital plane
	MRIcopyHeader(mri_talheader, mri_cc_tal);
	MRItalairachVoxelToVoxelEx(mri_cc_tal, cc_tal_x, cc_tal_y, cc_tal_z, &xv, &yv, &zv, lta) ;

	//rotate the cc volume by the rotation matrix calculated above
	MRIcopyHeader(mri_header, mri_cc);
	MRItoTalairachEx(mri_wm, mri_cc, lta2);
	// binalize the rotated cc volume (mri_cc)
  MRIbinarize(mri_cc, mri_cc, CC_VAL/2-1, 0, 100) ;

	MRIvoxelToTalairachVoxelEx(mri_cc, xv, yv, zv, &xc, &yc, &zc, lta2) ;

	//find the mid-sagital plane there
	xi=nint(xc);
	fprintf(stdout,"cc center is found at %d %d %d\n",xi, nint(yc),nint(zc));

	//find the bounding box
	for (y = 0 ; y < mri_cc->height ; y++) 
	{ 
		for (z = 0 ; z < mri_cc->depth ; z++) 
		{ 	 
			if ( MRIvox(mri_cc, xi, y, z) )
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

	//count the number if equally segmented five parts	
	for (i = xi-dxi ; i <= xi+dxi ; i++) 
	{		
		for ( j = 0 ; j <= 255 ; j++) 
		{ 
			for ( k = 0 ; k <= 255 ; k++) 
			{
				if ( MRIvox(mri_cc, i, j, k)>0) 
				{
					if ( k>=zi_low-10 )
					{
						temp = floor((k-zi_low)/((zi_high-zi_low+1)/5));
						if (temp < 0) temp = 0;
						if (temp >= 5) temp = 4;
						volume[temp] +=1 ;
						MRIvox(mri_cc, i, j, k)=(temp+1)*20+10;
					}
				} 	
			}
		} 
	}


	fprintf(fp, "%s %d %d %d %d %d %d %d %d %d \n", argv[1], nint(volume[4]), nint(volume[3]),nint(volume[2]), nint(volume[1]),nint(volume[0]), yi_low, yi_high, zi_low, zi_high);
	fprintf(stdout, "%s %d %d %d %d %d %d %d %d %d \n", argv[1], nint(volume[4]), nint(volume[3]),nint(volume[2]), nint(volume[1]),nint(volume[0]), yi_low, yi_high, zi_low, zi_high);

	sprintf(ofname,"%s/%s/mri/cc.mgz",data_dir,argv[1]) ; 
	fprintf(stdout, "writing corpus callosum output to %s...\n", ofname) ; 
	MRIwrite(mri_cc, ofname) ;

	MRIfree(&mri_tal) ;
	MRIfree(&mri_wm); 
	MRIfree(&mri_cc); 
	MRIfree(&mri_talheader);
	MRIfree(&mri_header);
	MRIfree(&mri_cc_tal) ; 
	MatrixFree(&mtrans);
	MatrixFree(&mrot);
	msec = TimerStop(&then) ; 
	fprintf(stdout, "corpus callosum matter segmentation took %2.1f minutes\n", (float)msec/(1000.0f*60.0f)); 
	fclose(fp);	
	exit(0) ; 
	return(0) ; 
} 

#define CC_SPREAD       10
#define MIN_THICKNESS   3
#define MAX_THICKNESS   20

static int
find_corpus_callosum(MRI *mri_tal, Real *pccx, Real *pccy, Real *pccz, const LTA *lta, MRI *mri_cc_tal)
{
  int         xv, yv, zv, max_y, max_thick=0, thickness=0, y1, xcc, ycc, x, y,x0, extension=50 ;
	int         flag=0, counts=0;
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

	fprintf(stdout, "original seed found at x=%d, y=%d z=%d \n", x0, yv, zv );
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
							fprintf(stdout, "potential cc found at (%d, %d), thickness = %d\n",
											xcc, ycc, thickness) ;
					}
					else if ( y==max_y && thickness== max_thick && (x==xcc+1 || flag ==1 ))
					{
						flag = 1;
						counts ++;
					}
					else if (flag == 1)
					{
						xcc = xcc+ nint(counts/2);
						flag = 0;
						counts = 0;
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
	fprintf(stdout, "%d, %d, %d\n", xcc, ycc, zv);

  find_cc_slice(mri_tal, pccx, pccy, pccz, lta, mri_cc_tal) ;

  return(NO_ERROR) ;
}


static int
find_cc_slice(MRI *mri_tal, Real *pccx, Real *pccy, Real *pccz, const LTA *lta, MRI *mri_cc_tal)
{
  // here we can handle only up to .5 mm voxel size
  int         area[MAX_SLICES*2], flag[MAX_SLICES*2], min_area, min_slice, slice, offset,xv,yv,zv,
              xo, yo ,i, total_area=0, left=0, right=0;
  MRI         *mri_slice, *mri_filled ;
  Real        aspect, x_tal, y_tal, z_tal, x, y, z, xvv, yvv, zvv;
  MRI_REGION  region ;
  char        fname[STRLEN] ;
  int         half_slices, ii, jj;
  double      voxsize = findMinSize(mri_tal);
  int         slice_size = mri_tal->width;
  int         max_slices = ceil(MAX_SLICES/voxsize);
  int         max_cc_area = ceil(MAX_CC_AREA/(voxsize*voxsize));
  int         min_cc_area = floor(MIN_CC_AREA/(voxsize*voxsize));

  half_slices = floor(HALF_SLICES/voxsize);
  if ( half_slices <= 0)
    half_slices = 1;

  x_tal = *pccx ; y_tal = *pccy ; z_tal = *pccz ;
  offset = 0 ;
  xo = yo = (slice_size-1)/2 ;  /* center point of the slice */
  for (slice = 0 ; slice < max_slices ; slice++)
  {
    offset = slice - half_slices ;

		i=0; area[slice]=0;
		while (area[slice]<100 && i<=5)
		{
			x = x_tal + offset ; y = y_tal ; z = z_tal ;
			MRIworldToVoxel(mri_tal, x, y,  z, &x, &y, &z) ;
			xv = nint(x) ; yv = nint(y)-i ; zv = nint(z) ;
			mri_slice = MRIextractPlane(mri_tal, NULL, MRI_SAGITTAL, xv);
			mri_filled =  MRIfillFG(mri_slice, NULL, zv, yv,0,WM_MIN_VAL,CC_VAL,&area[slice]);
			MRIboundingBox(mri_filled, 1, &region) ;
			aspect = (Real)region.dy / (Real)region.dx ;
			if(i++) fprintf(stdout,"moved %d in slice %d \n", i-1, slice);
		}
    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= slice_size -1 ||
        region.y+region.dy >= slice_size-1)
      area[slice] = 0 ;

		if ( !(area[slice]>1100&&(nint(y)-region.y>11)) || region.dy>=3.5*(y-region.y) )
		{	
			mri_filled = remove_fornix(mri_filled,xv,yv,zv);    
			
			area[slice] = 0;
			flag[slice] = 0;
			
			for (ii = 0 ; ii < mri_filled->width ; ii++) 
			{ 
				for (jj = 0 ; jj < mri_filled->height ; jj++) 
				{ 	 
					if ( MRIvox(mri_filled, ii, jj, 0)>0 ) area[slice]++;
				} 
			} 
		}
		else 
		{
			flag[slice] = 1;
			//if (offset>=-5&&offset<0) left++;
			//else if (offset<=5&&offset>0) right++;
		}
		
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "slice[%d] @ (%d, %d, %d): area = %d\n", 
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
  
#if 0
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
#else
  min_area = 10000*5 ; min_slice = -1 ;
  for (slice = 6 ; slice <= 14 ; slice++)
  {
		for (i=-2, total_area =0; i <=2; i++) 
			total_area += area[slice+i];

    if (total_area < min_area && 
        (total_area >= min_cc_area*5 && total_area <= max_cc_area*5))
    {
      min_area = total_area ;
      min_slice = slice ;
    }
  }
#endif
  
  /* couldn't find a good slice - don't update estimate */
  if (min_slice < 0)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "%s: could not find valid seed for the cc",
                 Progname));
  
  offset = floor((min_slice - half_slices)/2) ;
  fprintf(stdout, "find offset as %d using area\n", offset);

	//another way to move the central slice
  for (slice = 0 ; slice < max_slices ; slice++)
  {
	if ( (slice-(half_slices+offset)>=-6) && (slice<half_slices+offset) && flag[slice]==1 ) left++;
	else if ( (slice-(half_slices+offset)<=6) && (slice>half_slices+offset) && flag[slice]==1 ) right++;
  }
  offset = offset+left-right;
  fprintf(stdout, "find offset as %d using shifting\n", left-right);
  if (abs(offset)>=5) offset = 5*(offset/abs(offset));
  *pccx = x = x_tal+floor(offset) ; *pccy = y = y_tal ; *pccz = z = z_tal ;
  
  // just for debugging
  MRIworldToVoxel(mri_tal, x, y,  z, &xvv, &yvv, &zvv) ;
  *pccx = xvv ; *pccy = yvv ; *pccz = zvv ;

  fprintf(stdout, "updating initial cc seed to Tal vol (%.2f, %.2f, %.2f) TAL (%.2f, %.2f, %.2f)\n",
	    xvv, yvv, zvv, x, y, z);

  return(NO_ERROR) ;
}


static MRI *
remove_fornix(MRI *mri_filled, int xv, int yv, int zv)
{
	int    x, y, xi_low=255, xi_high=0, yi_low=255, yi_high=0, edge_count = 0, length=0;
	int    temp=0, temp2=0, old_temp =0;
	int    x1=0, y1=0, x2=0, y2=0, height=mri_filled->height, width=mri_filled->width, flag =0;
	int    i, section1[150];
	MRI *mri_temp1, *mri_temp2;
	char        fname[STRLEN] ;

	mri_temp1=MRIcopy(mri_filled,NULL);
	mri_temp2=MRIcopy(mri_filled,NULL);
	MRIvalueFill(mri_temp1, 0) ;	
	MRIvalueFill(mri_temp2, 0) ;	

	for (x = 1 ; x < width-1 ; x++) 
	{ 
		for (y = 1 ; y < height-1 ; y++) 
		{ 	 
			if ( MRIvox(mri_filled, x-1, y, 0) || MRIvox(mri_filled, x, y-1, 0) || MRIvox(mri_filled, x, y+1, 0) || MRIvox(mri_filled, x+1, y, 0) || MRIvox(mri_filled, x, y, 0) )
				MRIvox(mri_temp1,x,y,0)=100; 	
		} 
	}

	for (x = 1 ; x < width-1 ; x++) 
	{ 
		for (y = 1 ; y < height-1 ; y++) 
		{ 	 
			if ( MRIvox(mri_temp1, x-1, y, 0) && MRIvox(mri_temp1, x, y-1, 0) && MRIvox(mri_temp1, x, y+1, 0) && MRIvox(mri_temp1, x+1, y, 0) && MRIvox(mri_temp1, x, y, 0) )
				MRIvox(mri_temp2,x,y,0)=100; 	
		} 
	}


	for (x = 0 ; x < width ; x++) 
	{ 
		for (y = 0 ; y < height ; y++) 
		{ 	 
			if ( MRIvox(mri_temp2, x, y, 0) )
			{
				if (x < xi_low)
					xi_low = x ;
				if (y < yi_low)
					yi_low = y ;
				if (x > xi_high )
					xi_high = x ;
				if (y > yi_high)
					yi_high = y ;
			} 	
		} 
	}

	if (yi_high>yi_low+50) 
	{
		yi_high=yi_low+40;
		for (x = 0 ; x < width ; x++) 
		{
			for (y = yi_high ; y < height ; y++) 
			{ 	 
				MRIvox(mri_temp2, x, y, 0) = 0;
				MRIvox(mri_temp1, x, y, 0) = 0;
			}
		}
	}

	sprintf(fname, "/space/neo/2/recon/buckner/001015_vc5442/mri/cc_dilation.mgz");
	//MRIwrite(mri_temp1, fname) ;
	mri_filled=MRIcopy(mri_temp2,NULL);
	sprintf(fname, "/space/neo/2/recon/buckner/001015_vc5442/mri/cc_filled.mgz");
	//MRIwrite(mri_temp2, fname) ;

	/*find the first edge of the spike */

	x_edge =0;
	y_edge =0;

	/*the first  edge of the spike */
	flag=0;
  for (x_edge=xi_high-nint((xi_high-xi_low)/4); x_edge>=xi_low+nint((xi_high-xi_low)/7); x_edge--)
	{
		edge_count=0; y_edge=0;
		while ( edge_count<3 && y_edge<height-2 )	
		{
			length = edge_detection(mri_temp2,edge_count,0);
			if (length>1) edge_count++ ; 
			if (edge_count==1&&length>1) 
			{
				if (length>=20&&x_edge>xi_low+15&&x_edge<zv+20 ) 
				{
					temp=old_temp;
					if (x1<x_edge)
					{			
						y1=old_temp;
						x1=x_edge;
						for (y=temp; y<256; y++)
						MRIvox(mri_filled, x_edge, y, 0) =0;
					}
					flag = 2;
				}
				else	temp=y_edge;
				
				if ( length<=1 || y_edge<yv-3 )
				{
					//for (y=y_edge+1; y>0; y--)
					//	MRIvox(mri_filled, x_edge, y, 0) =0;
					edge_count-=1;
				}
				else if (length>2&&x_edge>xi_low+15&&x_edge<zv+17)
				{
					for (y=y_edge+1; y<256; y++)
						MRIvox(mri_filled, x_edge, y, 0) =0;
				}
			}
			else if (length>1&&edge_count>1 && y_edge<yi_high+1 && y_edge>yv+15 )
				edge_count -=1;
		}
		if (edge_count>=2&&flag==0) flag=1;
		else if (edge_count<=1&&flag==1&&x_edge>xi_low+13) 
		{
			flag = 0;
			y1=old_temp;
			x1=x_edge;
			//			if (x1<zv+13) break;
		}
		old_temp = temp;
	}
	//fprintf(stdout, "first point found at %d %d \n", x1, y1);
	x_edge =0;
	y_edge =0;

	/*the second edge of the spike */

	flag=0;
	//for (y_edge=yi_high-nint((yi_high-yi_low)/3); y_edge>=yi_low+4; y_edge--)
	for (y_edge=yv+20; y_edge>=yv; y_edge--)
	{
		edge_count=0; x_edge=0;
		i=yv+20-y_edge; 
		section1[i]=0;
		while (x_edge<width-1)	
		{
			length = edge_detection(mri_temp2,edge_count,1);
			if (length >=2) 				edge_count++ ;
			if (edge_count==1) 
			{
				temp=x_edge;
				if (!section1[i]) section1[i]=length;
			}
			if (edge_count==2) temp2=x_edge-length;
		}

		if (edge_count>=3&&flag==0) 				flag=1;
		else if (edge_count<=2&&flag==1) 
		{
			flag = 0;
			x2=old_temp;
			y2=y_edge;
			if (y2<=yi_low+20) break;
		}
		else if( x2==0&&i>=0&&(section1[i]>=19||(4*section1[i-1]<3*section1[i]))  )
		{
			x2=old_temp;
			y2=y_edge;
			if (y2<=yi_low+22) break;
		}
		if ( (edge_count>=4) && temp2<x1) old_temp =  temp2-1;
		else		old_temp = temp;
	}
	//fprintf(stdout, "second point found at %d %d \n", x2, y2);

	if ( x2>0 && x1>xi_low+8 && x1>x2 && x1<zv+5)
	{
		if (x2<x1)
		{
			temp=x1; x1=x2; x2=temp;
			temp=y1; y1=y2; y2=temp;
		}
		for (x=x1; x<=x2; x++)
		{
			for (y=nint(y1+(y2-y1)*(x-x1)/(x2-x1))+1;y<height;y++) MRIvox(mri_filled, x, y, 0)=0;
		}
	}

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "/space/neo/2/recon/buckner/010601_vc6977/mri/cc_cut.mgz");
    MRIwrite(mri_filled, fname) ;
  }
	MRIfree(&mri_temp1);
	MRIfree(&mri_temp2);	
	return(mri_filled) ;
}


static int
edge_detection(MRI *mri_temp, int edge_count, int signal)
{
	int length = 0, gap = 0;

	if (signal==1)
	{
		while (gap<=2)
		{
			gap=0;
			while ( x_edge < 256)
			{	
				if (MRIvox(mri_temp, x_edge, y_edge, 0))  
					break ;
				x_edge++;
			}
			
			while ( x_edge < 256 )
			{
				if (!MRIvox(mri_temp, x_edge, y_edge, 0))  
					break ;
				else length++ ;
				x_edge++;
			}	

			while ( x_edge < 256)
			{	
				if (MRIvox(mri_temp, x_edge, y_edge, 0))  
					break ;
				else gap++;
				x_edge++;
			}
			
			if (gap<=2&&x_edge<256) length += gap;
			else 
			{
				x_edge -= gap;
				break;
			}
		}
	}
	else
	{
		while (y_edge < 256)
		{	
			if (MRIvox(mri_temp, x_edge, y_edge, 0))  
				break ;
			y_edge++;
		}
		
		while ( y_edge < 256 )
		{
			if (!MRIvox(mri_temp, x_edge, y_edge, 0))  
				break ;
			else length++ ;
			y_edge++;
		}		
	}
	return(length);	
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
		fprintf(stdout, 
						"usage: %s <input volumes> <output volume>\n", 
						Progname) ; 
		exit(1) ; 
		break ;
	case 'T':
		dxi = atoi(argv[2]);
		fprintf(stdout,"change thickness to %d mm\n", 2*dxi+1);
		nargs = 1;
		break; 
	default: 
		fprintf(stdout, "unknown option %s\n", argv[1]) ; 
		exit(1) ; 
		break ; 
	} 
	
	return(nargs) ; 
}






