/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef FIELD_CODE_INCLUDED
#define FIELD_CODE_INCLUDED

#define NUMBER_OF_VECTORIAL_FIELDS 14

/* field code: definition of the field type*/
#define INFLATED_CURV_CORR_FRAME       0
#define SULC_CORR_FRAME                1
#define CURVATURE_CORR_FRAME           2
#define GRAYMID_CORR_FRAME             3
#define T1MID_CORR_FRAME               4
#define T2MID_CORR_FRAME               5
#define PDMID_CORR_FRAME               6
#define AMYGDALA_CORR_FRAME            7
#define HIPPOCAMPUS_CORR_FRAME         8
#define PALLIDUM_CORR_FRAME            9
#define PUTAMEN_CORR_FRAME            10
#define CAUDATE_CORR_FRAME            11
#define LAT_VENTRICLE_CORR_FRAME      12
#define INF_LAT_VENTRICLE_CORR_FRAME  13
#define OVERLAY_FRAME                 14
#define DISTANCE_TRANSFORM_FRAME      15

/* surface names */
#define INFLATED_CURVATURE_NAME         NULL      /* directly computed */
#define SULC_NAME                      "sulc"
#define CURVATURE_NAME                 NULL       /* directly computed */
/* GRAYMID_NAME  should already have been defined in mrisurf.h  */
#ifndef GRAYMID_NAME
#define GRAYMID_NAME                   "graymid"
#endif
#define  T1MID_NAME                   "T1mid"
#define  T2MID_NAME                   "T2mid"
#define  PDMID_NAME                   "PDmid"
#define  AMYGDALA_DIST_NAME           "amygdala_dist"
#define  HIPPOCAMPUS_DIST_NAME        "hippocampus_dist"
#define  PALLIDUM_DIST_NAME           "pallidum_dist"
#define  PUTAMEN_DIST_NAME            "putamen_dist"
#define  CAUDATE_DIST_NAME            "caudate_dist"
#define  LAT_VENTRICLE_DIST_NAME      "latventricle_dist"
#define  INF_LAT_VENTRICLE_DIST_NAME  "inflatventricle_dist"
#define   OVERLAY_NAME                 "overlay"
#define   DISTANCE_TRANSFORM_NAME      "distance_transform"
#ifndef MAX_OF_TWO
#define MAX_OF_TWO(a,b) ((a)>(b) ? (a):(b))
#endif

const char *ReturnFieldName(int which_field);
int IsDistanceField(int which_field);

typedef struct
{
  int     field;              /* see field code above */
  int     frame;              /* corresponding frame in mrisp */
  int     type;               /* field type (default,distance field,...) */
  float   l_corr;             /* correlation coefficient */
  float   l_pcorr;            /* polar correlation coefficient */
  float   sse;                /* sse associated with this field */
  char    *name ;             // if user specified
  int     navgs ;
  int     which_norm ;
}
FIELD_LABEL;

int InitFieldLabel(FIELD_LABEL *label);

int SetFieldLabel(FIELD_LABEL *label,
                  int field,
                  int frame,
                  float l_corr,
                  float l_pcorr,
                  int navgs,
                  int which_norm);

int SetFieldName(FIELD_LABEL *label, char *name) ;

#endif
