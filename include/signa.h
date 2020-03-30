/**
 * @brief signa mri file format 
 *
 * Defines the sizes, image file offsets, and the header word
 * offsets into the study, series, image and dss headers.
 */
/*
 * Original Author: G.E. Medical Systems, Gregory L. Meyers
 *
 */


#ifndef SIGNA_H
#define SIGNA_H

#include "mri.h"

/*       "DESC:1
 *  ***********************************************************************
 *
 *         INCLUDE FILE    IDBM_HEADER_DEF
 *
 *         AUTHOR          G.E. Medical Systems
 *                         NMR Software Engineering
 *                         Gregory L. Meyers
 *
 *         PURPOSE
 *
 *                 Defines the sizes, image file offsets, and the header word
 *         offsets into the study, series, image and dss headers.
 *
 *         INCLUDE FILES   None
 *
 *         DETAILS         None
 *
 *         GENERATION      DATE            AUTHOR NAME(S)
 *         ----------      ------------    --------------
 *
 *         01.01.00        July 28, 1983   Gregory L. Meyers
 *         01.06.01        August 17, 1984 Gregory L. Meyers
 *         01.06.02        January 22, 1985 Julie A. Donnell
 *
 * This file ('idbm_hdr_def.h')  Converted to C Language header file
 * format.  W.M. Leue  11-8-85
 *
 *  ***********************************************************************
 *        "ETX
 */

/* ---------------------------------------*/
/* Define the offsets in the study header.*/
/* ---------------------------------------*/


#define STHDR_ID        0       /* Study Header Identifier*/
#define STHDR_REV       7       /* Study Header Revision Number*/
#define STHDR_BLKS      11      /* Number of Study Header Blocks*/
#define STHDR_CRTRP     12      /* Study Header Creator (Process*/
#define STHDR_CRTRT     28      /* Study Header Creator (Task*/
#define STHDR_RAWNM     29      /* Raw Data Study Number*/
#define STHDR_STNUM     32      /* Study Number*/
#define STHDR_RAWID     35      /* Raw Data System ID*/
#define STHDR_SGENID    37      /* System Generation ID*/
#define STHDR_DATE      39      /* Date of Study (ascii*/
#define STHDR_IDATE     44      /* Date of Study (integer*/
#define STHDR_TIME      47      /* Time of Study (ascii*/
#define STHDR_ITIME     51      /* Time of Study (integer*/
#define STHDR_PNM       54      /* Patient Name*/
#define STHDR_PID       70      /* Patient ID*/
#define STHDR_PIDTMP    76      /* Patient ID padding for future exp.*/
#define STHDR_AGE       78      /* Age of patient*/
#define STHDR_SEX       80      /* Sex of patient*/
#define STHDR_WGHT      81      /* Weight of the patient in grams*/
#define STHDR_RFR       83      /* Refered by*/
#define STHDR_DGN       99      /* Diognostician*/
#define STHDR_OP        115     /* Operator*/
#define STHDR_DESC      131     /* Description*/
#define STHDR_HIST      161     /* History*/
#define STHDR_STIME     221     /* Creation time in seconds.*/
#define STHDR_HOSP      223     /* Hospital name*/
#define STHDR_RSRV1     239     /* GE NMR Reserved Area*/
#define STHDR_RSRV2     256     /* GE NMR Reserved Area*/
#define STHDR_CHECK     511     /* Study Header Checksum*/



/* ----------------------------------------*/
/* Define the offsets in the series header.*/
/* ----------------------------------------*/



#define SEHDR_ID        0  /* Series Header Identifier*/
#define SEHDR_REV       7  /* Series Header Revision Number*/
#define SEHDR_BLKS      11  /* Number of Series Header Blocks*/
#define SEHDR_CRTRP     12  /* Series Header Creator (Proc*/
#define SEHDR_CRTRT     28  /* Series Header Creator (Task*/
#define SEHDR_RAWNM     29  /* Original Series Number*/
#define SEHDR_SERNUM    31  /* Series Number*/
#define SEHDR_RAWID     33  /* Raw Data System ID*/
#define SEHDR_SGENID    35  /* System Generation ID*/
#define SEHDR_DATE      37  /* Date of series (ascii*/
#define SEHDR_IDATE     42  /* Date of series (integer*/
#define SEHDR_TIME      45  /* Time of Series (ascii*/
#define SEHDR_ITIME     49  /* Time of Series (integer*/
#define SEHDR_DESC      52  /* Series Description*/
#define SEHDR_TYPE      112  /* Series Type*/
#define SEHDR_CTYPE     113  /* Coil Type*/
#define SEHDR_CNAME     114  /* Coil Name*/
#define SEHDR_CNTRDESC  122  /* Contrast Description*/
#define SEHDR_PTYPE     138  /* Plane Type*/
#define SEHDR_PNAME     139  /* Plane Name*/
#define SEHDR_IMODE     147  /* Image Mode*/
#define SEHDR_FSTREN    148  /* Magnetic Field Strength*/
#define SEHDR_PSEQ      149  /* Pulse Sequence*/
#define SEHDR_PSSTYPE   150  /* Pulse sequence subtype*/
#define SEHDR_FOV       151  /* Field of view*/
#define SEHDR_CENTER    153  /* Center*/
#define SEHDR_C_R       153
#define SEHDR_C_A       155
#define SEHDR_C_S       157

#define SEHDR_ORIEN     159  /* Orientation*/
#define SEHDR_POS       160  /* Position*/
#define SEHDR_ANREF     161  /* Longitudinal Anotomical Reference*/
#define SEHDR_VANREF    177  /* Vertical Anotomical Reference*/
#define SEHDR_VERLAN    193  /* Vertical Landmark*/
#define SEHDR_HORLAN    195  /* Horizontal Landmark*/
#define SEHDR_TBLLOC    197  /* Physical Table Location*/
#define SEHDR_SMATRIX   199  /* Scan Matrix*/
#define SEHDR_IMATRIX   201  /* Image Matrix*/
#define SEHDR_IALLOC    202  /* No. of Images Allocated*/
#define SEHDR_GTYP      203  /* Gating Type*/
#define SEHDR_RSRV1     204  /* GE NMR Reserved*/
#define SEHDR_RSRV2     256  /* GE NMR Reserved*/
#define SEHDR_CHECK     511  /* Checksum for Series Header*/



/* --------------------------------*/
/* Define the image header offsets.*/
/* --------------------------------*/


#define IHDR_ID         0  /* Image Header Identifier*/
#define IHDR_REV        7  /* Image Header Revision Number*/
#define IHDR_BLKS       11  /* Number of Image Header Blocks*/
#define IHDR_CRTRP      12  /* Image Header Creator (Proc*/
#define IHDR_CRTRT      28  /* Image Header Creator (Task*/
#define IHDR_DATE       29  /* Image Creation Date (ascii*/
#define IHDR_IDATE      34  /* Image Creation Date (integer*/
#define IHDR_TIME       37  /* Image Creation Time (ascii*/
#define IHDR_ITIME      41  /* Image Creation Time (integer*/
#define IHDR_IMNUM      44  /* Image Number*/
#define IHDR_SERNM      46  /* Series Number of Image*/
#define IHDR_RAWID      48  /* Raw Data System ID*/
#define IHDR_SGENID     50  /* System Generation ID*/
#define IHDR_STRTX      52  /* Start Location X, Right min*/
#define IHDR_ENDX       54  /* End Location X, Right max*/
#define IHDR_STRTY      56  /* Start Location Y, Anterior min*/
#define IHDR_ENDY       58  /* End Location Y, Anterior max*/
#define IHDR_STRTZ      60  /* Start Location Z, Superior min*/
#define IHDR_ENDZ       62  /* End Location Z, Superior max*/
#define IHDR_OBLIQUE    64  /* Reserved for future use.*/
#define IHDR_LOCATN     73  /* Image Location*/
#define IHDR_TBLPOS     75  /* Table Position*/
#define IHDR_THICK      77  /* Thickness*/
#define IHDR_SPACE      79  /* Spacing*/
#define IHDR_ROUND      81  /* Round*/
#define IHDR_TR         82  /* Repititon/Recovery Time*/
#define IHDR_TS         84  /* Scan Time*/
#define IHDR_TE         86  /* Echo Delay*/
#define IHDR_TI         88  /* Inversion Time*/
#define IHDR_TY         90  /* Reserved for future use.*/
#define IHDR_NECHO      98  /* Number of echos.*/
#define IHDR_ECHON      99  /* Echo number.*/
#define IHDR_SLQUANT    100  /* Number of slices in scan group.*/
#define IHDR_NAVE       101  /* Number of averages.*/
#define IHDR_RSRCH      102  /* Research mode used ?*/
#define IHDR_PNAME      103  /* Name of PSD file.*/
#define IHDR_PSDDT      119/* Creation Date of PSD file.*/
#define IHDR_GPRE       125  /* Graphically Prescribed ?*/
#define IHDR_PSERIES    126  /* Prescribed Series Numbers*/
#define IHDR_PIMAGES    131  /* Prescribed Image Numbers*/
#define IHDR_SHAPE      136  /* Image Shape*/
#define IHDR_X          137  /* X pixel dimension*/
#define IHDR_Y          138  /* Y pixel dimension*/
#define IHDR_PIXSIZ     139  /* Pixel Size*/
#define IHDR_CMPRS      141  /* Image Compressed ?*/
#define IHDR_BITPIX     142  /* Bits per Pixel*/
#define IHDR_WINDOW     143  /* Default Window*/
#define IHDR_LEVEL      144  /* Default Level*/
#define IHDR_IFBLKS     145  /* Number of Blocks in File*/
#define IHDR_NEX        146  /* Number of excitations - real*/
#define IHDR_PSAR       148  /* Value of peak SAR - real*/
#define IHDR_ASAR       150  /* Value of average SAR - real*/
#define IHDR_MONITOR    152  /* SAR monitored ?*/
#define IHDR_CONTIG     153  /* Contiguous slices ?*/
#define IHDR_HRT_RT     154  /* Cardiac Heart Rate*/
#define IHDR_DEL_TRG    155  /* Total Delay Time After Trigger*/
#define IHDR_ARR        157  /* Arrhythmia Rejection Ratio*/
#define IHDR_RTIME      158  /* Cardiac Rep Time*/
#define IHDR_IMGS_PCY   159  /* Images per Cardiac Cycle*/
#define IHDR_ARRS_SCN   160  /* Number of ARR's during the Scan*/
#define IHDR_XMTATTN    162  /* Transmit attenuator setting*/
#define IHDR_RCVATTN    163  /* Recieve attenuator setting*/
#define IHDR_FLDSTR     164  /* Magnetic Field Strength*/
#define IHDR_RSRV1      166  /* GE NMR Reserved*/
#define IHDR_RSRV2      256  /* GE NMR Reserved*/
#define IHDR_CHECK      511  /* Image Header Checksum*/


/* ------------------------------*/
/* Define the DSS header offsets.*/
/* ------------------------------*/


#define DSSHDR_ID       0  /* DSS Header Identifier*/
#define DSSHDR_REV      7  /* DSS Header Revision Number*/
#define DSSHDR_BLKS     11  /* Number of DSS Header Blocks*/
#define DSSHDR_CRTRP    12  /* DSS Header Creator (process*/
#define DSSHDR_CRTRT    28  /* DSS Header Creator (task*/
#define DSSHDR_DATE     29  /* Date of Creation (ascii*/
#define DSSHDR_IDATE    34  /* Date of Creation (integer*/
#define DSSHDR_TIME     37  /* Time of Creation (ascii*/
#define DSSHDR_ITIME    41  /* Time of Creation (integer*/
#define DSSHDR_TLM      44  /* Time the DSSHDR was last modified*/
#define DSSHDR_ASTAT    46  /* Study Archive Status*/
#define DSSHDR_APEND    47  /* Study archive pending count.*/
#define DSSHDR_SPEND    48  /* Substructure archive pending count.*/
#define DSSHDR_LOCK     49  /* Resource locking queue desc*/
#define DSSHDR_STLST    57  /* IDBM Formated Study List Entry*/
#define DSSHDR_STIME    121  /* Creation time in seconds.*/
#define DSSHDR_MSERIES  123  /* Maximum series number created to date */
#define DSSHDR_SERIES   124  /* Bit map for series in study.*/
#define DSSHDR_DEFSER   126  /* Default series and iii for the study.*/
#define DSSHDR_EOF      129  /* Points to the next block to allocate.*/
#define DSSHDR_SDIR     130  /* Series number and sgenid referenced by max */
#define DSSHDR_IMAP     223  /* Image block number map for each series*/
#define DSSHDR_RSRV1    471  /* GE NMR Reserved Area*/


/* ---------------------------------------------------------------*/
/* Define the offsets for each series entry in the DSS file series*/
/* directory blocks.                                              */
/* ---------------------------------------------------------------*/


#define DSSHDR_SERIES_WPENTRY   33  /* Words per series entry. */
#define DSSHDR_SERIES_MIMAGE    0  /* Maximum image number created to date.*/
#define DSSHDR_SERIES_IMAGE     1  /* Number of images in the series.*/
#define DSSHDR_SERIES_RSRV      2  /* Number of images reserved for space.*/
#define DSSHDR_SERIES_PTYPE     3  /* Series plane type, AXIAL, SAG ...*/
#define DSSHDR_SERIES_PSEQ      4  /* Pulse sequence type. */
#define DSSHDR_SERIES_ANREF     5  /* Anatomical reference, 2 chars.*/
#define DSSHDR_SERIES_VANREF    6  /* Vertical Anatomical reference.*/
#define DSSHDR_SERIES_FOV       7  /* Field of view. in I*2 mm * 10 */
#define DSSHDR_SERIES_DESC      8  /* First 22 chars of series description.*/
#define DSSHDR_SERIES_ASTAT     19  /* Archive status.*/
#define DSSHDR_SERIES_APEND     20  /* Archive pending status.*/
#define DSSHDR_SERIES_IMATRIX   21  /* Image matrix size.*/
#define DSSHDR_SERIES_MINLOC    22  /* Minimum location, I*2 mm * 10.*/
#define DSSHDR_SERIES_MAXLOC    23  /* Maximim location, I*2 mm * 10.*/
#define DSSHDR_SERIES_LOCCHAR   24  /* Characters whcih represent +,- pos.*/
#define DSSHDR_SERIES_TLM       25  /* Time the series was last modified.*/
#define DSSHDR_SERIES_CHECK     27  /* Checksum for the series header.*/
#define DSSHDR_SERIES_DEFIM     28  /* Default image for the series.*/
#define DSSHDR_SERIES_TIMAGE    29  /* Total image count (partial & complete*/


/* ------------------------------------------------------------------------*/
/* Define the offsets for each image entry in the DSS file image directory */
/* blocks.                                                                 */
/* ------------------------------------------------------------------------*/

#define DSSHDR_IMAGE_WPENTRY    16
#define DSSHDR_IMAGE_LOCATN     0  /* Image location in I*2 mm * 10 */
#define DSSHDR_IMAGE_THICK      1  /* Image thickness in I*2 mm * 10*/
#define DSSHDR_IMAGE_NEXT       2  /* Next image by location.*/
#define DSSHDR_IMAGE_PRIOR      3  /* Prior image by location.*/
#define DSSHDR_IMAGE_EXIST      4  /* None, header or data exists. */
#define DSSHDR_IMAGE_T1         5
#define DSSHDR_IMAGE_T2         7
#define DSSHDR_IMAGE_TRGD       9  /* Trigger Delay */

/* ------------------------------------------------------ */
/* Define Block Offsets for the Various Header Components */
/* ------------------------------------------------------ */

#define SYSCON_START    0*256
#define SITCUS_START    4*256
#define STHDR_START     6*256
#define SEHDR_START     8*256
#define IHDR_START      10*256
#define RDBM_START      12*256
#define PSD_START       16*256
#define PIXMAP_START    26*256
#define IDATA_START     28*256

/* ------------------------------------------------------- */
/* Define Two Useful Offsets in the SYS_CONFIGURATION File */
/* ------------------------------------------------------- */

#define SCON_SYSID      6               /* System ID String */
#define SCON_HNAME      16              /* Hospital Name */

/* -------------------*/
/* End of Header File */
/* -------------------*/

MRI *read_signa(char *fname, char *h, unsigned char **im,float scale);
#define HLENGTH  IDATA_START*2
#define IMA_HLENGTH  6144
#define IMA2_HLENGTH  4096

typedef struct HINFO_
{
  int imnr0,imnr1;
  int x,y,ptype;
  float tr,ti,te,strtx,endx,strty,endy,strtz,endz;
  float locatn,fov,center,psiz,thick;
  float c_r, c_a, c_s, orien,pos ;
  int  plane_type ;
  int  num_echoes ;
}
HINFO;

#define SIGNA_AXIAL      0
#define SIGNA_SAGITTAL   1
#define SIGNA_CORONAL    2
#define SIGNA_OBLIQUE    3
#define SIGNA_SCREENSAVE 4

char *ghc( char *destin, char *header, int offset, int byte_length );
int ghi( char *header, int offset );
float ghf( char *header, int offset );

int is_signa(char *fname) ;
MRI *signaRead(char *fname, int read_volume_flag) ;

#endif
