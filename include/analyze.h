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


#ifndef ANALYZE_H
#define ANALYZE_H

/*******************************************************************/
/*struct dsr--the ANALYZE .hdr struct                                   */
/********************************************************************/

#define DT_NONE          0
#define DT_BINARY        1
#define DT_UNSIGNED_CHAR 2
#define DT_SIGNED_SHORT  4
#define DT_SIGNED_INT    8
#define DT_FLOAT        16
#define DT_DOUBLE       64

typedef struct
{
  int sizeof_hdr;                 /*required--byte size of header file*/
  char data_type[10];
  char db_name[18];
  int extents;                    /*required--16384*/
  short int session_error;
  char regular;                   /*required--'r'=regular*/
  char hkey_un0;
}
header_key;

typedef struct
{
  short int dim[8];               /*required*/
  char vox_units[4];
  char cal_units[8];
  short int unused1;
  short int datatype;             /*required: 0=unk,1=1 bit/pix,2=8bits,4=16 bits*/
  /*8=32 bits (signed int),16=32 bits (floating pt)*/
  /*32=64 bits (2 floats),64=64 bits (double)     */
  short int bitpix;               /*bits/pixel*/
  short int dim_un0;
  float pixdim[8];                /*real world values of dimensions mm ms*/
  float vox_offset;
  float roi_scale;
  float funused1;
  float funused2;
  float cal_max;
  float cal_min;
  int compressed;
  int verified;
  int glmax,glmin;                /*required*/
}
image_dimension;

typedef struct
{
  char descrip[80];               /*Will be displayed when loading*/
  char aux_file[24];
  char orient;
  unsigned char originator[10];
  char generated[10];
  char scannum[10];
  char patient_id[10];
  char exp_date[10];
  char exp_time[10];
  char hist_un0[3];
  int views;
  int vols_added;
  int start_field;
  int field_skip;
  int omax,omin;
  int smax,smin;
}
data_history;

typedef struct
{
  header_key hk;
  image_dimension dime;
  data_history hist;
}
dsr;


#endif
