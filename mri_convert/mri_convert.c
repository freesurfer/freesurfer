#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include "mri.h"
#include "error.h"
#include "mri_identify.h"
#include "utils.h"
#include "transform.h"
#include "mrimorph.h"
#include "DICOMRead.h"
#include "unwarpGradientNonlinearity.h"
#include "version.h"

/* ----- determines tolerance of non-orthogonal basis vectors ----- */
#define CLOSE_ENOUGH  (5e-3)

void get_ints(int argc, char *argv[], int *pos, int *vals, int nvals);
void get_floats(int argc, char *argv[], int *pos, float *vals, int nvals);
void get_string(int argc, char *argv[], int *pos, char *val);
void usage_message(FILE *stream);
void usage(FILE *stream);
float findMinSize(MRI *mri);

int debug=0;

extern int errno;

char *Progname;

int main(int argc, char *argv[])
{
  int nargs = 0;
  MRI *mri_unwarped;
  MRI *mri, *mri2, *template, *mri_in_like;
  int i;
  int reorder_vals[3];
  float invert_val;
  int in_info_flag, out_info_flag;
  int template_info_flag;
  int conform_flag;
  int conform_min;  // conform to the smallest dimension
  int parse_only_flag;
  int reorder_flag;
  int in_stats_flag, out_stats_flag;
  int read_only_flag, no_write_flag;
  char in_name[STRLEN], out_name[STRLEN];
  int in_volume_type, out_volume_type;
  char resample_type[STRLEN];
  int resample_type_val;
  int in_i_size_flag, in_j_size_flag, in_k_size_flag;
  int out_i_size_flag, out_j_size_flag, out_k_size_flag;
  float in_i_size, in_j_size, in_k_size;
  float out_i_size, out_j_size, out_k_size;
  int sizes_good_flag;
  float in_i_directions[3], in_j_directions[3], in_k_directions[3];
  float out_i_directions[3], out_j_directions[3], out_k_directions[3];
  int in_i_direction_flag, in_j_direction_flag, in_k_direction_flag;
  int out_i_direction_flag, out_j_direction_flag, out_k_direction_flag;
  int in_tr_flag = 0;
  float in_tr = 0;
  float magnitude;
  float i_dot_j, i_dot_k, j_dot_k;
  float in_center[3], out_center[3];
  int in_center_flag, out_center_flag;
  int out_data_type;
  char out_data_type_string[STRLEN];
  int out_n_i, out_n_j, out_n_k;
  int out_n_i_flag, out_n_j_flag, out_n_k_flag;
  float fov_x, fov_y, fov_z;
  int force_in_type_flag, force_out_type_flag;
  int forced_in_type, forced_out_type;
  char in_type_string[STRLEN], out_type_string[STRLEN];
  char subject_name[STRLEN];
  char reslice_like_name[STRLEN];
  int reslice_like_flag;
  int frame_flag;
  int frame;
  char in_name_only[STRLEN];
  char transform_fname[STRLEN];
  int transform_flag, invert_transform_flag;
  LTA *lta_transform;
  M3D *m3d_transform;
  MRI *mri_transformed = NULL;
  int transform_type;
  MATRIX *inverse_transform_matrix;
  int smooth_parcellation_flag, smooth_parcellation_count;
  int in_like_flag;
  char in_like_name[STRLEN];
  int in_n_i, in_n_j, in_n_k;
  int in_n_i_flag, in_n_j_flag, in_n_k_flag;
  int fill_parcellation_flag;
  int read_parcellation_volume_flag;
  int zero_outlines_flag;
  int read_otl_flags;
  int color_file_flag;
  char color_file_name[STRLEN];
  int no_scale_flag;
  int temp_type;
  int roi_flag;
  FILE *fptmp;
  int j,translate_labels_flag;
  int force_ras_good = FALSE;
  char gdf_image_stem[STRLEN];
  int in_matrix_flag, out_matrix_flag;
  float minSize;

  for(i=0;i<argc;i++) printf("%s ",argv[i]);
  printf("\n");
  fflush(stdout);

  for(i=0;i<argc;i++){
    if(strcmp(argv[i],"--debug")==0){
      fptmp = fopen("debug.gdb","w");
      fprintf(fptmp,"# source this file in gdb to debug\n");
      fprintf(fptmp,"file %s \n",argv[0]);
      fprintf(fptmp,"run ");
      for(j=1;j<argc;j++){
  if(strcmp(argv[j],"--debug")!=0)
    fprintf(fptmp,"%s ",argv[j]);
      }
      fprintf(fptmp,"\n");
      fclose(fptmp);
      break;
    }
  }

  /* ----- keep the compiler quiet ----- */
  mri2 = NULL;
  forced_in_type = forced_out_type = MRI_VOLUME_TYPE_UNKNOWN;
  invert_transform_flag = FALSE;

  /* ----- get the program name ----- */
  Progname = strrchr(argv[0], '/');
  Progname = (Progname == NULL ? argv[0] : Progname + 1);

  /* ----- pass the command line to mriio ----- */
  mriio_command_line(argc, argv);

  /* ----- catch no arguments here ----- */
  if(argc == 1)
  {
    usage(stdout);
    exit(1);
  }

  /* ----- initialize values ----- */
  in_name[0] = '\0';
  out_name[0] = '\0';
  invert_val = -1.0;
  in_info_flag = FALSE;
  out_info_flag = FALSE;
  conform_flag = TRUE;
  parse_only_flag = FALSE;
  reorder_flag = FALSE;
  in_stats_flag = FALSE;
  out_stats_flag = FALSE;
  read_only_flag = FALSE;
  no_write_flag = FALSE;
  resample_type_val = RESAMPLE_INTERPOLATE;
  in_i_size_flag = in_j_size_flag = in_k_size_flag = FALSE;
  out_i_size_flag = out_j_size_flag = out_k_size_flag = FALSE;
  in_i_direction_flag = in_j_direction_flag = in_k_direction_flag = FALSE;
  out_i_direction_flag = out_j_direction_flag = out_k_direction_flag = FALSE;
  in_center_flag = FALSE;
  out_center_flag = FALSE;
  out_data_type = -1;
  out_n_i_flag = out_n_j_flag = out_n_k_flag = FALSE;
  template_info_flag = FALSE;
  out_volume_type = MRI_CORONAL_SLICE_DIRECTORY;
  force_in_type_flag = force_out_type_flag = FALSE;
  subject_name[0] = '\0';
  reslice_like_flag = FALSE;
  frame_flag = FALSE;
  transform_flag = FALSE;
  smooth_parcellation_flag = FALSE;
  in_like_flag = FALSE;
  in_n_i_flag = in_n_j_flag = in_n_k_flag = FALSE;
  fill_parcellation_flag = FALSE;
  zero_outlines_flag = FALSE;
  color_file_flag = FALSE;
  no_scale_flag = FALSE;
  roi_flag = FALSE;
  translate_labels_flag = TRUE;
  gdf_image_stem[0] = '\0';
  in_matrix_flag = FALSE;
  out_matrix_flag = FALSE;
  conform_min = FALSE;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_convert.c,v 1.52 2003/07/03 22:00:30 tosa Exp $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  for(i = 1;i < argc;i++)
  {
    if(strcmp(argv[i], "-version2") == 0)
      exit(97);
    if(strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reorder") == 0)
    {
      get_ints(argc, argv, &i, reorder_vals, 3);
      reorder_flag = TRUE;
    }
    else if(strcmp(argv[i], "--debug") == 0) debug = 1;
    else if(strcmp(argv[i], "--invert_contrast") == 0)
      get_floats(argc, argv, &i, &invert_val, 1);
    else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input_volume") == 0)
      get_string(argc, argv, &i, in_name);
    else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output_volume") == 0)
      get_string(argc, argv, &i, out_name);
    else if(strcmp(argv[i], "-nc") == 0 || strcmp(argv[i], "--no_conform") == 0)
      conform_flag = FALSE;
    else if (strcmp(argv[i], "-cm") == 0 || strcmp(argv[i], "--conform_min") == 0)
      conform_min = TRUE;
    else if(strcmp(argv[i], "-po") == 0 || strcmp(argv[i], "--parse_only") == 0)
      parse_only_flag = TRUE;
    else if(strcmp(argv[i], "-ii") == 0 || strcmp(argv[i], "--in_info") == 0)
      in_info_flag = TRUE;
    else if(strcmp(argv[i], "-oi") == 0 || strcmp(argv[i], "--out_info") == 0)
      out_info_flag = TRUE;
    else if(strcmp(argv[i], "-ti") == 0 || strcmp(argv[i], "--template_info") == 0)
      template_info_flag = TRUE;
    else if(strcmp(argv[i], "-is") == 0 || strcmp(argv[i], "--in_stats") == 0)
      in_stats_flag = TRUE;
    else if(strcmp(argv[i], "-os") == 0 || strcmp(argv[i], "--out_stats") == 0)
      out_stats_flag = TRUE;
    else if(strcmp(argv[i], "-ro") == 0 || strcmp(argv[i], "--read_only") == 0)
      read_only_flag = TRUE;
    else if(strcmp(argv[i], "-nw") == 0 || strcmp(argv[i], "--no_write") == 0)
      no_write_flag = TRUE;
    else if(strcmp(argv[i], "-im") == 0 || strcmp(argv[i], "--in_matrix") == 0)
      in_matrix_flag = TRUE;
    else if(strcmp(argv[i], "-om") == 0 || strcmp(argv[i], "--out_matrix") == 0)
      out_matrix_flag = TRUE;
    else if(strcmp(argv[i], "--force_ras_good") == 0) force_ras_good = TRUE;
    else if(strcmp(argv[i], "-at") == 0 || strcmp(argv[i], "--apply_transform") == 0 || strcmp(argv[i], "-T") == 0)
    {
      get_string(argc, argv, &i, transform_fname);
      transform_flag = TRUE;
      invert_transform_flag = FALSE;
    }
    else if(strcmp(argv[i], "-ait") == 0 || strcmp(argv[i], "--apply_inverse_transform") == 0)
    {
      get_string(argc, argv, &i, transform_fname);
      transform_flag = TRUE;
      invert_transform_flag = TRUE;
    }
    else if(strcmp(argv[i], "-iis") == 0 || 
      strcmp(argv[i], "--in_i_size") == 0)
    {
      get_floats(argc, argv, &i, &in_i_size, 1);
      in_i_size_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ijs") == 0 || strcmp(argv[i], "--in_j_size") == 0)
    {
      get_floats(argc, argv, &i, &in_j_size, 1);
      in_j_size_flag = TRUE;
    }
    else if(strcmp(argv[i], "-iks") == 0 || strcmp(argv[i], "--in_k_size") == 0)
    {
      get_floats(argc, argv, &i, &in_k_size, 1);
      in_k_size_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ois") == 0 || strcmp(argv[i], "--out_i_size") == 0)
    {
      get_floats(argc, argv, &i, &out_i_size, 1);
      out_i_size_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ojs") == 0 || strcmp(argv[i], "--out_j_size") == 0)
    {
      get_floats(argc, argv, &i, &out_j_size, 1);
      out_j_size_flag = TRUE;
    }
    else if(strcmp(argv[i], "-oks") == 0 || strcmp(argv[i], "--out_k_size") == 0)
    {
      get_floats(argc, argv, &i, &out_k_size, 1);
      out_k_size_flag = TRUE;
    }
    else if(strcmp(argv[i], "-iid") == 0 || strcmp(argv[i], "--in_i_direction") == 0)
    {
      get_floats(argc, argv, &i, in_i_directions, 3);
      magnitude = sqrt(in_i_directions[0]*in_i_directions[0] + in_i_directions[1]*in_i_directions[1] + in_i_directions[2]*in_i_directions[2]);
      if(magnitude == 0.0)
      {
        fprintf(stderr, "\n%s: directions must have non-zero magnitude; in_i_direction = (%g, %g, %g)\n", Progname, in_i_directions[0], in_i_directions[1], in_i_directions[2]);
        usage_message(stdout);
        exit(1);
      }
      if(magnitude != 1.0)
      {
        printf("normalizing in_i_direction: (%g, %g, %g) -> ", in_i_directions[0], in_i_directions[1], in_i_directions[2]);
        in_i_directions[0] = in_i_directions[0] / magnitude;
        in_i_directions[1] = in_i_directions[1] / magnitude;
        in_i_directions[2] = in_i_directions[2] / magnitude;
        printf("(%g, %g, %g)\n", in_i_directions[0], in_i_directions[1], in_i_directions[2]);
      }
      in_i_direction_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ijd") == 0 || strcmp(argv[i], "--in_j_direction") == 0)
    {
      get_floats(argc, argv, &i, in_j_directions, 3);
      magnitude = sqrt(in_j_directions[0]*in_j_directions[0] + in_j_directions[1]*in_j_directions[1] + in_j_directions[2]*in_j_directions[2]);
      if(magnitude == 0.0)
      {
        fprintf(stderr, "\n%s: directions must have non-zero magnitude; in_j_direction = (%g, %g, %g)\n", Progname, in_j_directions[0], in_j_directions[1], in_j_directions[2]);
        usage_message(stdout);
        exit(1);
      }
      if(magnitude != 1.0)
      {
        printf("normalizing in_j_direction: (%g, %g, %g) -> ", in_j_directions[0], in_j_directions[1], in_j_directions[2]);
        in_j_directions[0] = in_j_directions[0] / magnitude;
        in_j_directions[1] = in_j_directions[1] / magnitude;
        in_j_directions[2] = in_j_directions[2] / magnitude;
        printf("(%g, %g, %g)\n", in_j_directions[0], in_j_directions[1], in_j_directions[2]);
      }
      in_j_direction_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ikd") == 0 || strcmp(argv[i], "--in_k_direction") == 0)
    {
      get_floats(argc, argv, &i, in_k_directions, 3);
      magnitude = sqrt(in_k_directions[0]*in_k_directions[0] + in_k_directions[1]*in_k_directions[1] + in_k_directions[2]*in_k_directions[2]);
      if(magnitude == 0.0)
      {
        fprintf(stderr, "\n%s: directions must have non-zero magnitude; in_k_direction = (%g, %g, %g)\n", Progname, in_k_directions[0], in_k_directions[1], in_k_directions[2]);
        usage_message(stdout);
        exit(1);
      }
      if(magnitude != 1.0)
      {
        printf("normalizing in_k_direction: (%g, %g, %g) -> ", in_k_directions[0], in_k_directions[1], in_k_directions[2]);
        in_k_directions[0] = in_k_directions[0] / magnitude;
        in_k_directions[1] = in_k_directions[1] / magnitude;
        in_k_directions[2] = in_k_directions[2] / magnitude;
        printf("(%g, %g, %g)\n", in_k_directions[0], in_k_directions[1], in_k_directions[2]);
      }
      in_k_direction_flag = TRUE;
    }
    else if(strcmp(argv[i], "-oid") == 0 || strcmp(argv[i], "--out_i_direction") == 0)
    {
      get_floats(argc, argv, &i, out_i_directions, 3);
      magnitude = sqrt(out_i_directions[0]*out_i_directions[0] + out_i_directions[1]*out_i_directions[1] + out_i_directions[2]*out_i_directions[2]);
      if(magnitude == 0.0)
      {
        fprintf(stderr, "\n%s: directions must have non-zero magnitude; out_i_direction = (%g, %g, %g)\n", Progname, out_i_directions[0], out_i_directions[1], out_i_directions[2]);
        usage_message(stdout);
        exit(1);
      }
      if(magnitude != 1.0)
      {
        printf("normalizing out_i_direction: (%g, %g, %g) -> ", out_i_directions[0], out_i_directions[1], out_i_directions[2]);
        out_i_directions[0] = out_i_directions[0] / magnitude;
        out_i_directions[1] = out_i_directions[1] / magnitude;
        out_i_directions[2] = out_i_directions[2] / magnitude;
        printf("(%g, %g, %g)\n", out_i_directions[0], out_i_directions[1], out_i_directions[2]);
      }
      out_i_direction_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ojd") == 0 || strcmp(argv[i], "--out_j_direction") == 0)
    {
      get_floats(argc, argv, &i, out_j_directions, 3);
      magnitude = sqrt(out_j_directions[0]*out_j_directions[0] + out_j_directions[1]*out_j_directions[1] + out_j_directions[2]*out_j_directions[2]);
      if(magnitude == 0.0)
      {
        fprintf(stderr, "\n%s: directions must have non-zero magnitude; out_j_direction = (%g, %g, %g)\n", Progname, out_j_directions[0], out_j_directions[1], out_j_directions[2]);
        usage_message(stdout);
        exit(1);
      }
      if(magnitude != 1.0)
      {
        printf("normalizing out_j_direction: (%g, %g, %g) -> ", out_j_directions[0], out_j_directions[1], out_j_directions[2]);
        out_j_directions[0] = out_j_directions[0] / magnitude;
        out_j_directions[1] = out_j_directions[1] / magnitude;
        out_j_directions[2] = out_j_directions[2] / magnitude;
        printf("(%g, %g, %g)\n", out_j_directions[0], out_j_directions[1], out_j_directions[2]);
      }
      out_j_direction_flag = TRUE;
    }
    else if(strcmp(argv[i], "-okd") == 0 || strcmp(argv[i], "--out_k_direction") == 0)
    {
      get_floats(argc, argv, &i, out_k_directions, 3);
      magnitude = sqrt(out_k_directions[0]*out_k_directions[0] + out_k_directions[1]*out_k_directions[1] + out_k_directions[2]*out_k_directions[2]);
      if(magnitude == 0.0)
      {
        fprintf(stderr, "\n%s: directions must have non-zero magnitude; out_k_direction = (%g, %g, %g)\n", Progname, out_k_directions[0], out_k_directions[1], out_k_directions[2]);
        usage_message(stdout);
        exit(1);
      }
      if(magnitude != 1.0)
      {
        printf("normalizing out_k_direction: (%g, %g, %g) -> ", out_k_directions[0], out_k_directions[1], out_k_directions[2]);
        out_k_directions[0] = out_k_directions[0] / magnitude;
        out_k_directions[1] = out_k_directions[1] / magnitude;
        out_k_directions[2] = out_k_directions[2] / magnitude;
        printf("(%g, %g, %g)\n", out_k_directions[0], out_k_directions[1], out_k_directions[2]);
      }
      out_k_direction_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ic") == 0 || strcmp(argv[i], "--in_center") == 0)
    {
      get_floats(argc, argv, &i, in_center, 3);
      in_center_flag = TRUE;
    }
    else if(strcmp(argv[i], "-oc") == 0 || strcmp(argv[i], "--out_center") == 0)
    {
      get_floats(argc, argv, &i, out_center, 3);
      out_center_flag = TRUE;
    }
    else if(strcmp(argv[i], "-oni") == 0 || strcmp(argv[i], "-oic") == 0 || strcmp(argv[i], "--out_i_count") == 0)
    {
      get_ints(argc, argv, &i, &out_n_i, 1);
      out_n_i_flag = TRUE;
    }
    else if(strcmp(argv[i], "-onj") == 0 || strcmp(argv[i], "-ojc") == 0 || strcmp(argv[i], "--out_j_count") == 0)
    {
      get_ints(argc, argv, &i, &out_n_j, 1);
      out_n_j_flag = TRUE;
    }
    else if(strcmp(argv[i], "-onk") == 0 || strcmp(argv[i], "-okc") == 0 || strcmp(argv[i], "--out_k_count") == 0)
    {
      get_ints(argc, argv, &i, &out_n_k, 1);
      out_n_k_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ini") == 0 || strcmp(argv[i], "-iic") == 0 || strcmp(argv[i], "--in_i_count") == 0)
    {
      get_ints(argc, argv, &i, &in_n_i, 1);
      in_n_i_flag = TRUE;
    }
    else if(strcmp(argv[i], "-inj") == 0 || strcmp(argv[i], "-ijc") == 0 || strcmp(argv[i], "--in_j_count") == 0)
    {
      get_ints(argc, argv, &i, &in_n_j, 1);
      in_n_j_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ink") == 0 || strcmp(argv[i], "-ikc") == 0 || strcmp(argv[i], "--in_k_count") == 0)
    {
      get_ints(argc, argv, &i, &in_n_k, 1);
      in_n_k_flag = TRUE;
    }
    else if( strcmp(argv[i], "-tr") == 0 )
    {
      get_floats(argc, argv, &i, &in_tr, 1);
      in_tr_flag = TRUE;
    }

    else if(strcmp(argv[i], "-odt") == 0 || strcmp(argv[i], "--out_data_type") == 0)
    {
      get_string(argc, argv, &i, out_data_type_string);
      if(strcmp(StrLower(out_data_type_string), "uchar") == 0)
        out_data_type = MRI_UCHAR;
      else if(strcmp(StrLower(out_data_type_string), "short") == 0)
        out_data_type = MRI_SHORT;
      else if(strcmp(StrLower(out_data_type_string), "int") == 0)
        out_data_type = MRI_INT;
      else if(strcmp(StrLower(out_data_type_string), "float") == 0)
        out_data_type = MRI_FLOAT;
      else
      {
        fprintf(stderr, "\n%s: unknown data type \"%s\"\n", Progname, argv[i]);
        usage_message(stdout);
        exit(1);
      }
    }
    else if(strcmp(argv[i], "-rt") == 0 || strcmp(argv[i], "--resample_type") == 0)
    {
      get_string(argc, argv, &i, resample_type);
      if(strcmp(StrLower(resample_type), "interpolate") == 0)
        resample_type_val = RESAMPLE_INTERPOLATE;
      else if(strcmp(StrLower(resample_type), "nearest") == 0)
        resample_type_val = RESAMPLE_NEAREST;
      else if(strcmp(StrLower(resample_type), "weighted") == 0)
        resample_type_val = RESAMPLE_WEIGHTED;
      else if(strcmp(StrLower(resample_type), "sinc") == 0)
        resample_type_val = RESAMPLE_SINC;
      else
      {
        fprintf(stderr, "\n%s: unknown resample type \"%s\"\n", Progname, argv[i]);
        usage_message(stdout);
        exit(1);
      }
    }
    else if(strcmp(argv[i], "-it") == 0 || strcmp(argv[i], "--in_type") == 0)
    {
      get_string(argc, argv, &i, in_type_string);
      forced_in_type = string_to_type(in_type_string);
      force_in_type_flag = TRUE;
    }
    else if(strcmp(argv[i], "-ot") == 0 || strcmp(argv[i], "--out_type") == 0)
    {
      get_string(argc, argv, &i, out_type_string);
      forced_out_type = string_to_type(out_type_string);/* see mri_identify.c */
      force_out_type_flag = TRUE;
    }
    else if(strcmp(argv[i], "-sn") == 0 || strcmp(argv[i], "--subject_name") == 0)
    {
      get_string(argc, argv, &i, subject_name);
    }
    else if(strcmp(argv[i], "-gis") == 0 || strcmp(argv[i], "--gdf_image_stem") == 0)
    {
      get_string(argc, argv, &i, gdf_image_stem);
      if(strlen(gdf_image_stem) == 0)
      {
        fprintf(stderr, "\n%s: zero length GDF image stem given\n", Progname);
        usage_message(stdout);
        exit(1);
      }
    }
    else if(strcmp(argv[i], "-rl") == 0 || strcmp(argv[i], "--reslice_like") == 0)
    {
      get_string(argc, argv, &i, reslice_like_name);
      reslice_like_flag = TRUE;
    }
    else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--frame") == 0)
    {
      get_ints(argc, argv, &i, &frame, 1);
      frame_flag = TRUE;
    }
    else if(strcmp(argv[i], "-il") == 0 || strcmp(argv[i], "--in_like") == 0)
    {
      get_string(argc, argv, &i, in_like_name);
      in_like_flag = TRUE;
    }
    else if(strcmp(argv[i], "-roi") == 0 || strcmp(argv[i], "--roi") == 0)
    {
      roi_flag = TRUE;
    }
    else if(strcmp(argv[i], "-fp") == 0 || strcmp(argv[i], "--fill_parcellation") == 0)
    {
      fill_parcellation_flag = TRUE;
    }
    else if(strcmp(argv[i], "-zo") == 0 || strcmp(argv[i], "--zero_outlines") == 0)
    {
      zero_outlines_flag = TRUE;
    }
    else if(strcmp(argv[i], "-sp") == 0 || strcmp(argv[i], "--smooth_parcellation") == 0)
    {
      get_ints(argc, argv, &i, &smooth_parcellation_count, 1);
      if(smooth_parcellation_count < 14 || smooth_parcellation_count > 26)
      {
        fprintf(stderr, "\n%s: clean parcellation count must be between 14 and 26, inclusive\n", Progname);
        usage_message(stdout);
        exit(1);
      }
      smooth_parcellation_flag = TRUE;
    }
    else if(strcmp(argv[i], "-cf") == 0 || strcmp(argv[i], "--color_file") == 0)
    {
      get_string(argc, argv, &i, color_file_name);
      color_file_flag = TRUE;
    }
    else if(strcmp(argv[i], "-nt") == 0 || strcmp(argv[i], "--no_translate") == 0)
    {
      translate_labels_flag = FALSE;
    }
    else if(strcmp(argv[i], "-ns") == 0 || strcmp(argv[i], "--no_scale") == 0)
    {
      get_ints(argc, argv, &i, &no_scale_flag, 1);
      no_scale_flag = (no_scale_flag == 0 ? FALSE : TRUE);
    }
    else if(strcmp(argv[i], "--unwarp_gradient_nonlinearity") == 0)
      {
  /* !@# start */
  unwarp_flag = 1;

  /* Determine gradient type: sonata or allegra */
  get_string(argc, argv, &i, unwarp_gradientType);
  if( (strcmp(unwarp_gradientType, "sonata")  != 0) &&
      (strcmp(unwarp_gradientType, "allegra") != 0) &&
      (strcmp(unwarp_gradientType, "GE") != 0) )
    {
      fprintf(stderr, "\n%s: must specify gradient type ('sonata' or 'allegra' or 'GE')\n", Progname);
      usage_message(stdout);
      exit(1);
    }
  
  /* Determine whether or not to do a partial unwarp */
  get_string(argc, argv, &i, unwarp_partialUnwarp);
  if( (strcmp(unwarp_partialUnwarp, "fullUnwarp") != 0) &&
      (strcmp(unwarp_partialUnwarp, "through-plane")  != 0) )
    {
      fprintf(stderr, "\n%s: must specify unwarping type ('fullUnwarp' or 'through-plane')\n", Progname);
      usage_message(stdout);
      exit(1);
    }

  /* Determine whether or not to do jacobian correction */
  get_string(argc, argv, &i, unwarp_jacobianCorrection);
  if( (strcmp(unwarp_jacobianCorrection, "JacobianCorrection")  != 0) &&
      (strcmp(unwarp_jacobianCorrection, "noJacobianCorrection") != 0) )
    {
      fprintf(stderr, "\n%s: must specify intensity correction type ('JacobianCorrection' or 'noJacobianCorrection')\n", Progname);
      usage_message(stdout);
      exit(1);
    }

  /* Determine interpolation type: linear or sinc */
  get_string(argc, argv, &i, unwarp_interpType);
  if( (strcmp(unwarp_interpType, "linear") != 0) &&
      (strcmp(unwarp_interpType, "sinc")   != 0) )
    {
      fprintf(stderr, "\n%s: must specify interpolation type ('linear' or 'sinc')\n", Progname);
      usage_message(stdout);
      exit(1);
    }

  /* Get the HW for sinc interpolation (if linear interpolation,
           this integer value is not used) */
  get_ints(argc, argv, &i, &unwarp_sincInterpHW, 1);

  /* OPTIONS THAT THERE ARE NO PLANS TO SUPPORT */
  /* Jacobian correction with through-plane only correction */ 
  if( (strcmp(unwarp_jacobianCorrection, "JacobianCorrection") == 0) &&
      (strcmp(unwarp_partialUnwarp, "through-plane") == 0) )      
    {
      fprintf(stderr, "\n%s: Jacobian correction not valid for 'through-plane' only unwarping)\n", Progname);
      exit(1);
    }

  /* OPTIONS NOT CURRENTLY SUPPORTED (BUT W/ PLANS TO SUPPORT) */
  /* 1) GE unwarping not supported until we have offset data */
  if( strcmp(unwarp_gradientType, "GE") == 0 )
    {
      fprintf(stderr, "\n%s: unwarping data from GE scanners not supported at present.\n", Progname);
      exit(1);
    }    
  
  /* 2) for GE: through-plane correction requires rewarping the
           in-plane unwarped image, which requires map inversion */
  if( strcmp(unwarp_partialUnwarp, "through-plane") == 0 )      
    {
      fprintf(stderr, "\n%s: through-plane only unwarping not supported at present.\n", Progname);
      exit(1);
    }
  /* !@# end */
  
      }
    else if(strcmp(argv[i], "-u") == 0 || strcmp(argv[i], "--usage") == 0)
    {
      usage(stdout);
      exit(0);
    }
    /*-------------------------------------------------------------*/
    else if(strcmp(argv[i], "--status") == 0 || 
      strcmp(argv[i], "--statusfile") == 0)
    {
      /* File name to write percent complete for Siemens DICOM */
      if( (argc-1) - i < 1 ){
  fprintf(stderr,"ERROR: option --statusfile requires one argument\n");
  exit(1);
      }
      i++;
      SDCMStatusFile = (char *) calloc(strlen(argv[i])+1,sizeof(char));
      memcpy(SDCMStatusFile,argv[i],strlen(argv[i]));
      fptmp = fopen(SDCMStatusFile,"w");
      if(fptmp == NULL){
  fprintf(stderr,"ERROR: could not open %s for writing\n",
    SDCMStatusFile);
  exit(1);
      }
      fprintf(fptmp,"0\n");
      fclose(fptmp);
    }
    /*-------------------------------------------------------------*/
    else if(strcmp(argv[i], "--sdcmlist") == 0)
    {
      /* File name that contains a list of Siemens DICOM files
   that are in the same run as the one listed on the
   command-line. If not present, the directory will be scanned,
   but this can take a while.
      */
      if( (argc-1) - i < 1 ){
  fprintf(stderr,"ERROR: option --sdcmlist requires one argument\n");
  exit(1);
      }
      i++;
      SDCMListFile = (char *) calloc(strlen(argv[i])+1,sizeof(char));
      memcpy(SDCMListFile,argv[i],strlen(argv[i]));
      fptmp = fopen(SDCMListFile,"r");
      if(fptmp == NULL){
  fprintf(stderr,"ERROR: could not open %s for reading\n",
    SDCMListFile);
  exit(1);
      }
      fclose(fptmp);
    }
    /*-------------------------------------------------------------*/
    else if( (strcmp(argv[i], "--nspmzeropad") == 0) ||
       (strcmp(argv[i], "--out_nspmzeropad") == 0))
    {
      /* Choose the amount of zero padding for spm output files */
      if( (argc-1) - i < 1 ){
  fprintf(stderr,"ERROR: option --out_nspmzeropad requires one argument\n");
  exit(1);
      }
      i++;
      sscanf(argv[i],"%d",&N_Zero_Pad_Output);
    }
    /*-------------------------------------------------------------*/
    else if( (strcmp(argv[i], "--in_nspmzeropad") == 0))
    {
      /* Choose the amount of zero padding for spm input files */
      if( (argc-1) - i < 1 ){
  fprintf(stderr,"ERROR: option --in_nspmzeropad requires one argument\n");
  exit(1);
      }
      i++;
      sscanf(argv[i],"%d",&N_Zero_Pad_Input);
    }
    /*-------------------------------------------------------------*/
    else
    {
      if(argv[i][0] == '-')
      {
        fprintf(stderr, "\n%s: unknown flag \"%s\"\n", Progname, argv[i]);
        usage_message(stdout);
        exit(1);
      }
      else
      {
        if(in_name[0] == '\0')
          strcpy(in_name, argv[i]);
        else if(out_name[0] == '\0')
          strcpy(out_name, argv[i]);
        else
        {
          if(i + 1 == argc)
            fprintf(stderr, "\n%s: extra argument (\"%s\")\n", Progname, argv[i]);
          else
            fprintf(stderr, "\n%s: extra arguments (\"%s\" and following)\n", Progname, argv[i]);
          usage_message(stdout);
          exit(1);
        }
      }

    }
  }
  /**** Finished parsing command line ****/
  /* option inconsistency checks */
  if(force_ras_good && (in_i_direction_flag || in_j_direction_flag ||
      in_k_direction_flag)){
    fprintf(stderr, "ERROR: cannot use --force_ras_good and --in_?_direction_flag\n");
    exit(1);
  }
  if (conform_flag == FALSE && conform_min == TRUE)
  {
    fprintf(stderr, "You cannot use both -nc (--no_conform) and -cm (--conform_min) at the same time.\n");
    exit(1);
  }
  /* ----- catch zero or negative voxel dimensions ----- */

  sizes_good_flag = TRUE;

  if((in_i_size_flag && in_i_size <= 0.0) || 
     (in_j_size_flag && in_j_size <= 0.0) || 
     (in_k_size_flag && in_k_size <= 0.0) || 
     (out_i_size_flag && out_i_size <= 0.0) || 
     (out_j_size_flag && out_j_size <= 0.0) || 
     (out_k_size_flag && out_k_size <= 0.0))
  {
    fprintf(stderr, "\n%s: voxel sizes must be greater than zero\n", Progname);
    sizes_good_flag = FALSE;
  }

  if(in_i_size_flag && in_i_size <= 0.0)
    fprintf(stderr, "in i size = %g\n", in_i_size);
  if(in_j_size_flag && in_j_size <= 0.0)
    fprintf(stderr, "in j size = %g\n", in_j_size);
  if(in_k_size_flag && in_k_size <= 0.0)
    fprintf(stderr, "in k size = %g\n", in_k_size);
  if(out_i_size_flag && out_i_size <= 0.0)
    fprintf(stderr, "out i size = %g\n", out_i_size);
  if(out_j_size_flag && out_j_size <= 0.0)
    fprintf(stderr, "out j size = %g\n", out_j_size);
  if(out_k_size_flag && out_k_size <= 0.0)
    fprintf(stderr, "out k size = %g\n", out_k_size);

  if(!sizes_good_flag)
  {
    usage_message(stdout);
    exit(1);
  }

  /* ----- catch missing input or output volume name ----- */
  if(in_name[0] == '\0')
  {
    fprintf(stderr, "\n%s: missing input volume name\n", Progname);
    usage_message(stdout);
    exit(1);
  }

  if(out_name[0] == '\0' && !(read_only_flag || no_write_flag))
  {
    fprintf(stderr, "\n%s: missing output volume name\n", Progname);
    usage_message(stdout);
    exit(1);
  }

  /* ----- copy file name (only -- strip '@' and '#') ----- */
  MRIgetVolumeName(in_name, in_name_only);

  /* ----- catch unknown volume types ----- */
  if(force_in_type_flag && forced_in_type == MRI_VOLUME_TYPE_UNKNOWN)
  {
    fprintf(stderr, "\n%s: unknown input volume type %s\n", Progname, in_type_string);
    usage_message(stdout);
    exit(1);
  }

  /* ----- warn if read only is desired and an output volume is specified or the output info flag is set ----- */
  if(read_only_flag && out_name[0] != '\0')
    fprintf(stderr, "%s: warning: read only flag is set; nothing will be written to %s\n", Progname, out_name);
  if(read_only_flag && (out_info_flag || out_matrix_flag))
    fprintf(stderr, "%s: warning: read only flag is set; no output information will be printed\n", Progname);

  /* ----- catch the parse-only flag ----- */
  if(parse_only_flag)
  {

    printf("input volume name: %s\n", in_name);
    printf("input name only: %s\n", in_name_only);
    printf("output volume name: %s\n", out_name);
    printf("parse_only_flag = %d\n", parse_only_flag);
    printf("conform_flag = %d\n", conform_flag);
    printf("in_info_flag = %d\n", in_info_flag);
    printf("out_info_flag = %d\n", out_info_flag);
    printf("in_matrix_flag = %d\n", in_matrix_flag);
    printf("out_matrix_flag = %d\n", out_matrix_flag);

    if(force_in_type_flag)
      printf("input type is %d\n", forced_in_type);
    if(force_out_type_flag)
      printf("output type is %d\n", forced_out_type);

    if(subject_name[0] != '\0')
      printf("subject name is %s\n", subject_name);

    if(invert_val >= 0)
      printf("inversion, value is %g\n", invert_val);

    if(reorder_flag)
      printf("reordering, values are %d %d %d\n", reorder_vals[0], reorder_vals[1], reorder_vals[2]);

    printf("translation of otl labels is %s\n", translate_labels_flag ? "on" : "off");

    exit(0);

  }

  /* ----- get the type of the output ----- */
  if(!force_out_type_flag)
  {
    if(!read_only_flag && !no_write_flag)
    {
      out_volume_type = mri_identify(out_name);
      if(out_volume_type == MRI_VOLUME_TYPE_UNKNOWN)
      {
        fprintf(stderr, "%s: can't determine type of output volume\n", Progname);
        exit(1);
      }
    }
  }
  else
    out_volume_type = forced_out_type;

  /* ----- check for a gdf image stem if the output type is gdf ----- */
  if(out_volume_type == GDF_FILE && strlen(gdf_image_stem) == 0)
  {
    fprintf(stderr, "%s: GDF output type, but no GDF image file stem\n", Progname);
    exit(1);
  }

  /* ----- read the in_like volume ----- */
  if(in_like_flag)
  {
    printf("reading info from %s...\n", in_like_name);
    mri_in_like = MRIreadInfo(in_like_name);
    if(mri_in_like == NULL)
      exit(1);
  }

  /* ----- read the volume ----- */
  if(force_in_type_flag)
    in_volume_type = forced_in_type;
  else
    in_volume_type = mri_identify(in_name_only);

  if(in_volume_type == MRI_VOLUME_TYPE_UNKNOWN)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADFILE, "unknown file type for file %s", in_name_only);
    if(in_like_flag)
      MRIfree(&mri_in_like);
    exit(1);
  }

  if(roi_flag && in_volume_type != GENESIS_FILE)
  {
    errno = 0;
    ErrorPrintf(ERROR_BADPARM, "rois must be in GE format");
    if(in_like_flag) MRIfree(&mri_in_like);
    exit(1);
  }

  printf("reading from %s...\n", in_name_only);

  if(in_volume_type == OTL_FILE)
  {

    if(!in_like_flag && !in_n_k_flag)
    {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "parcellation read: must specify"
		  "a volume depth with either in_like or in_k_count");
      exit(1);
    }

    if(!color_file_flag)
    {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "parcellation read: must specify a"
		  "color file name");
      if(in_like_flag)
        MRIfree(&mri_in_like);
      exit(1);
    }

    read_parcellation_volume_flag = TRUE;
    if(read_only_flag && (in_info_flag || in_matrix_flag) && !in_stats_flag)
      read_parcellation_volume_flag = FALSE;

    read_otl_flags = 0x00;

    if(read_parcellation_volume_flag)
      read_otl_flags |= READ_OTL_READ_VOLUME_FLAG;

    if(fill_parcellation_flag)
      read_otl_flags |= READ_OTL_FILL_FLAG;
    else
    {
      printf("notice: unfilled parcellaions currently unimplemented\n");
      printf("notice: filling outlines\n");
      read_otl_flags |= READ_OTL_FILL_FLAG;
    }

    if(translate_labels_flag)
      read_otl_flags |= READ_OTL_TRANSLATE_LABELS_FLAG;

    if(zero_outlines_flag)
    {
/*
      printf("notice: zero outlines currently unimplemented\n");
      printf("notice: outlines won't be cleared\n");
*/
      read_otl_flags |= READ_OTL_ZERO_OUTLINES_FLAG;
    }

    if(in_like_flag)
      mri = MRIreadOtl(in_name, mri_in_like->width, mri_in_like->height, 
		       mri_in_like->depth, color_file_name, read_otl_flags);
    else
      mri = MRIreadOtl(in_name, 0, 0, in_n_k, color_file_name, read_otl_flags);

    if(mri == NULL)
    {
      if(in_like_flag)
        MRIfree(&mri_in_like);
      exit(1);
    }

    /* ----- smooth the parcellation if requested ----- */
    if(smooth_parcellation_flag)
    {
      printf("smoothing parcellation...\n");
      mri2 = MRIsmoothParcellation(mri, smooth_parcellation_count);
      if(mri2 == NULL)
      {
        if(in_like_flag)
          MRIfree(&mri_in_like);
        exit(1);
      }
      MRIfree(&mri);
      mri = mri2;
    }

    resample_type_val = RESAMPLE_NEAREST;
    no_scale_flag = TRUE;

  }
  else if(roi_flag)
  {
    if(!in_like_flag && !in_n_k_flag)
    {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "roi read: must specify a volume"
		  "depth with either in_like or in_k_count");
      if(in_like_flag)
        MRIfree(&mri_in_like);
      exit(1);
    }

    if(in_like_flag)
      mri = MRIreadGeRoi(in_name, mri_in_like->depth);
    else
      mri = MRIreadGeRoi(in_name, in_n_k);

    if(mri == NULL)
    {
      if(in_like_flag)
        MRIfree(&mri_in_like);
      exit(1);
    }

    resample_type_val = RESAMPLE_NEAREST;
    no_scale_flag = TRUE;
  }
  else
  {
    if(read_only_flag && (in_info_flag || in_matrix_flag) && !in_stats_flag)
      mri = MRIreadInfo(in_name);
    else
    {
      if(force_in_type_flag){
	//printf("MRIreadType()\n");
        mri = MRIreadType(in_name, in_volume_type);
      }
      else{
        mri = MRIread(in_name);
      }
    }

  }

  if(mri == NULL)
  {
    if(in_like_flag) MRIfree(&mri_in_like);
    exit(1);
  }

  if(unwarp_flag)
    {
      /* if unwarp_flag is true, unwarp the distortions due
   to gradient coil nonlinearities */
      printf("INFO: unwarping ... ");
      mri_unwarped = unwarpGradientNonlinearity(mri, 
            unwarp_gradientType, 
            unwarp_partialUnwarp,
            unwarp_jacobianCorrection,
            unwarp_interpType,
            unwarp_sincInterpHW);
      MRIfree(&mri);
      mri = mri_unwarped;
      printf("done \n ");      
    }

  printf("TR=%2.2f, te=%2.2f, flip angle=%2.2f\n",
         mri->tr, mri->te, DEGREES(mri->flip_angle)) ;
  if(in_volume_type != OTL_FILE)
  {
  if(fill_parcellation_flag)
    printf("fill_parcellation flag ignored on a non-parcellation read\n");
  if(smooth_parcellation_flag)
    printf("smooth_parcellation flag ignored on a non-parcellation read\n");
  }

  /* ----- apply the in_like volume if it's been read ----- */
  if(in_like_flag)
  {
    if(mri->width   != mri_in_like->width ||
       mri->height  != mri_in_like->height ||
       mri->depth   != mri_in_like->depth ||
       mri->nframes != mri_in_like->nframes)
    {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "volume sizes do not match\n");
      ErrorPrintf(ERROR_BADPARM, "%s: (width, height, depth, frames) = (%d, %d, %d, %d)\n", in_name, mri->width, mri->height, mri->depth, mri->nframes);
      ErrorPrintf(ERROR_BADPARM, "%s: (width, height, depth, frames) = (%d, %d, %d, %d)\n", in_like_name, mri_in_like->width, mri_in_like->height, mri_in_like->depth, mri_in_like->nframes);
      MRIfree(&mri);
      MRIfree(&mri_in_like);
      exit(1);
    }

    temp_type = mri->type;

    if(MRIcopyHeader(mri_in_like, mri) != NO_ERROR)
    {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "error copying information from %s structure to %s structure\n", in_like_name, in_name);
      MRIfree(&mri);
      MRIfree(&mri_in_like);
      exit(1);
    }

    mri->type = temp_type;

    MRIfree(&mri_in_like);

  }

  if(mri->ras_good_flag == 0){
    printf("WARNING: it does not appear that there was sufficient information\n"
     "in the input to assign orientation to the volume... \n");
    if(force_ras_good){
      printf("However, you have specified that the default orientation should\n"
       "be used with by adding --force_ras_good on the command-line.\n");
      mri->ras_good_flag = 1;
    }
    if(in_i_direction_flag || in_j_direction_flag || in_k_direction_flag){
      printf("However, you have specified one or more orientations on the \n"
       "command-line using -i?d or --in-?-direction (?=i,j,k).\n");
      mri->ras_good_flag = 1;
    }
  }

  /* ----- apply command-line parameters ----- */
  if(in_i_size_flag)    mri->xsize = in_i_size;
  if(in_j_size_flag)    mri->ysize = in_j_size;
  if(in_k_size_flag)    mri->zsize = in_k_size;
  if(in_i_direction_flag)
  {
    mri->x_r = in_i_directions[0];
    mri->x_a = in_i_directions[1];
    mri->x_s = in_i_directions[2];
    mri->ras_good_flag = 1;
  }
  if(in_j_direction_flag)
  {
    mri->y_r = in_j_directions[0];
    mri->y_a = in_j_directions[1];
    mri->y_s = in_j_directions[2];
    mri->ras_good_flag = 1;
  }
  if(in_k_direction_flag)
  {
    mri->z_r = in_k_directions[0];
    mri->z_a = in_k_directions[1];
    mri->z_s = in_k_directions[2];
    mri->ras_good_flag = 1;
  }
  if(in_center_flag)
  {
    mri->c_r = in_center[0];
    mri->c_a = in_center[1];
    mri->c_s = in_center[2];
  }
  if(subject_name[0] != '\0')
    strcpy(mri->subject_name, subject_name);

  if(in_tr_flag) mri->tr = in_tr;

  /* ----- correct starts, ends, and fov if necessary ----- */
  if(in_i_size_flag || in_j_size_flag || in_k_size_flag)
  {

    fov_x = mri->xsize * mri->width;
    fov_y = mri->ysize * mri->height;
    fov_z = mri->zsize * mri->depth;

    mri->xend = fov_x / 2.0;
    mri->xstart = -mri->xend;
    mri->yend = fov_y / 2.0;
    mri->ystart = -mri->yend;
    mri->zend = fov_z / 2.0;
    mri->zstart = -mri->zend;

    mri->fov = (fov_x > fov_y ? (fov_x > fov_z ? fov_x : fov_z) : (fov_y > fov_z ? fov_y : fov_z) );

  }

  /* ----- give a warning for non-orthogonal directions ----- */
  i_dot_j = mri->x_r * mri->y_r + mri->x_a * mri->y_a + mri->x_s * mri->y_s;
  i_dot_k = mri->x_r * mri->z_r + mri->x_a * mri->z_a + mri->x_s * mri->z_s;
  j_dot_k = mri->y_r * mri->z_r + mri->y_a * mri->z_a + mri->y_s * mri->z_s;
  if(fabs(i_dot_j) > CLOSE_ENOUGH || fabs(i_dot_k) > CLOSE_ENOUGH || fabs(i_dot_k) > CLOSE_ENOUGH)
  {
    printf("warning: input volume axes are not orthogonal\n");
  }
  printf("i_ras = (%g, %g, %g)\n", mri->x_r, mri->x_a, mri->x_s);
  printf("j_ras = (%g, %g, %g)\n", mri->y_r, mri->y_a, mri->y_s);
  printf("k_ras = (%g, %g, %g)\n", mri->z_r, mri->z_a, mri->z_s);

  /* ----- catch the in info flag ----- */
  if(in_info_flag)
  {
    printf("input structure:\n");
    MRIdump(mri, stdout);
  }

  if(in_matrix_flag)
  {
    MATRIX *i_to_r;
    i_to_r = extract_i_to_r(mri);
    if(i_to_r != NULL)
    {
      printf("input ijk -> ras:\n");
      MatrixPrint(stdout, i_to_r);
      MatrixFree(&i_to_r);
    }
    else
      printf("error getting input matrix\n");
  }

  /* ----- catch the in stats flag ----- */
  if(in_stats_flag)
    MRIprintStats(mri, stdout);

  /* ----- apply a transformation if requested ----- */
  if(transform_flag)
  {

    printf("INFO: Applying transformation from file %s...\n", transform_fname);

    if(!FileExists(transform_fname))
    {
      fprintf(stderr,"ERROR: cannot find transform file %s\n",transform_fname);
      exit(1);
    }

    transform_type = TransformFileNameType(transform_fname);
    if(transform_type == MNI_TRANSFORM_TYPE || 
       transform_type == TRANSFORM_ARRAY_TYPE)
    {
      printf("Reading transform\n");
      lta_transform = LTAread(transform_fname);
      if(lta_transform  == NULL){
    fprintf(stderr, "ERROR: Reading transform from file %s\n", 
      transform_fname);
    exit(1);
  }
      
      printf("Input Matrix --------------------------\n");
      MatrixPrint(stdout,lta_transform->xforms[0].m_L);
      printf("---------------------------------\n");

      if(invert_transform_flag)
      {
  inverse_transform_matrix = MatrixInverse(lta_transform->xforms[0].m_L,
             NULL);
        if(inverse_transform_matrix == NULL)
        {
          fprintf(stderr, "ERROR: inverting transform\n");
    MatrixPrint(stdout,lta_transform->xforms[0].m_L);
          exit(1);
        }

        MatrixFree(&(lta_transform->xforms[0].m_L));
        lta_transform->xforms[0].m_L = inverse_transform_matrix;
      }

      /* Think about calling MRIlinearTransform() here; need vox2vox
   transform. Can create NN version. In theory, LTAtransform()
         can handle multiple transforms, but the inverse assumes only
         one. NN is good for ROI*/

      printf("INFO: resampling input volume \n");
      printf("---------------------------------\n");
      printf("Resampling Matrix input to LTAtransform(): \n");
      MatrixPrint(stdout,lta_transform->xforms[0].m_L);
      printf("---------------------------------\n");

      if (lta_transform->type == LINEAR_RAS_TO_RAS){
  /* What does this do? M = M*V*W, where V is world2vox
     transform for COR, and W is vox2world transform for
     COR. Maybe it changes the center?*/
  printf("INFO: LTAvoxelTransformToCoronalRasTransform()\n");
  LTAvoxelTransformToCoronalRasTransform(lta_transform);
      }

#if 1
      /* LTAtransform() runs either MRIapplyRASlinearTransform() 
   for RAS2RAS or MRIlinearTransform() for Vox2Vox. Since
      the matrix is transformed to Vox2Vox, the LTAtransform line
      should give the same results as the MRIlinearTransfrom().*/
      mri_transformed = LTAtransform(mri, NULL, lta_transform);
      //mri_transformed = MRIlinearTransform(mri, NULL, 
      //           lta_transform->xforms[0].m_L);
#else
      /* This part is experimental. The section that uses NEAREST
   can be used for transforming ROIs that need to be resampled
   given the transform matrix. */
      if (lta_transform->type == LINEAR_RAS_TO_RAS){
  mri_transformed = 
    MRIapplyRASlinearTransform(mri,NULL,
             lta_transform->xforms[0].m_L);
      }
      else{
  printf("INFO: transforming using nearest\n");
  mri_transformed = 
    MRIlinearTransformInterp(mri, NULL,
           lta_transform->xforms[0].m_L,
           SAMPLE_NEAREST);
      }
#endif
      if(mri_transformed == NULL){
        fprintf(stderr, "ERROR: applying transform to volume\n");
        exit(1);
      }

      LTAfree(&lta_transform);
      MRIfree(&mri);
      mri = mri_transformed;
    }

    else if(transform_type == MORPH_3D_TYPE)
    {

      if((m3d_transform = MRI3DreadSmall(transform_fname)) == NULL)
      {
        fprintf(stderr, "error reading transform from file %s\n", 
    transform_fname);
        exit(1);
      }

      if(invert_transform_flag)
        mri_transformed = MRIapplyInverse3DMorph(mri, m3d_transform, NULL);
      else
        mri_transformed = MRIapply3DMorph(mri, m3d_transform, NULL);

      if(mri_transformed == NULL)
      {
        fprintf(stderr, "error applying transform\n");
        exit(1);
      }

      MRI3DmorphFree(&m3d_transform);

      MRIfree(&mri);
      mri = mri_transformed;

    }
    else
    {
      fprintf(stderr, "unknown transform type in file %s\n", transform_fname);
      exit(1);
    }

  }

  if(reslice_like_flag)
  {

    printf("reading template info from volume %s...\n", reslice_like_name);

    template = MRIreadInfo(reslice_like_name);
    if(template == NULL)
    {
      fprintf(stderr, "error reading from volume %s\n", reslice_like_name);
      exit(1);
    }

  }
  else
  {
    template = MRIallocHeader(mri->width, mri->height, mri->depth, mri->type);
    MRIcopyHeader(mri, template);
    if(conform_flag)
    {
      if(out_volume_type == MRI_CORONAL_SLICE_DIRECTORY)
      {
        template->width = template->height = template->depth = 256;
        template->imnr0 = 1;
        template->imnr1 = 256;
        template->type = MRI_UCHAR;
	if (conform_min==TRUE)
	{
	  // find out the min size 
	  minSize = findMinSize(mri);
	  template->thick = minSize;
	  template->ps = minSize;
	  template->xsize = template->ysize = template->zsize = minSize;
	  printf("Data is conformed to %g size for all directin\n", minSize); 
	}
	else
	{
	  template->thick = 1.0;
	  template->ps = 1.0;
	  template->xsize = template->ysize = template->zsize = 1.0;
	}
        template->xstart = template->ystart = template->zstart = -128.0;
        template->xend = template->yend = template->zend = 128.0;
        template->x_r = -1.0;  template->x_a =  0.0;  template->x_s =  0.0;
        template->y_r =  0.0;  template->y_a =  0.0;  template->y_s = -1.0;
        template->z_r =  0.0;  template->z_a =  1.0;  template->z_s =  0.0;
      }
    }
    else if(out_volume_type != MRI_CORONAL_SLICE_DIRECTORY)
        printf("the output volume is not a COR- directory."
         "The --no_conform (-nc) argument is not needed\n");

  }

  /* ----- apply command-line parameters ----- */
  if(out_i_size_flag)
    template->xsize = out_i_size;
  if(out_j_size_flag)
    template->ysize = out_j_size;
  if(out_k_size_flag)
    template->zsize = out_k_size;
  if(out_n_i_flag)
    template->width = out_n_i;
  if(out_n_j_flag)
    template->height = out_n_j;
  if(out_n_k_flag)
    template->depth = out_n_k;
  if(out_i_direction_flag)
  {
    template->x_r = out_i_directions[0];
    template->x_a = out_i_directions[1];
    template->x_s = out_i_directions[2];
  }
  if(out_j_direction_flag)
  {
    template->y_r = out_j_directions[0];
    template->y_a = out_j_directions[1];
    template->y_s = out_j_directions[2];
  }
  if(out_k_direction_flag)
  {
    template->z_r = out_k_directions[0];
    template->z_a = out_k_directions[1];
    template->z_s = out_k_directions[2];
  }
  if(out_center_flag)
  {
    template->c_r = out_center[0];
    template->c_a = out_center[1];
    template->c_s = out_center[2];
  }

  /* ----- correct starts, ends, and fov if necessary ----- */
  if(out_i_size_flag || out_j_size_flag || out_k_size_flag ||
     out_n_i_flag    || out_n_j_flag    || out_n_k_flag)
  {

    fov_x = template->xsize * template->width;
    fov_y = template->ysize * template->height;
    fov_z = template->zsize * template->depth;

    template->xend = fov_x / 2.0;
    template->xstart = -template->xend;
    template->yend = fov_y / 2.0;
    template->ystart = -template->yend;
    template->zend = fov_z / 2.0;
    template->zstart = -template->zend;

    template->fov = (fov_x > fov_y ? (fov_x > fov_z ? fov_x : fov_z) : (fov_y > fov_z ? fov_y : fov_z) );

  }

  /* ----- give a warning for non-orthogonal directions ----- */
  i_dot_j = template->x_r * template->y_r + template->x_a * template->y_a + template->x_s * template->y_s;
  i_dot_k = template->x_r * template->z_r + template->x_a * template->z_a + template->x_s * template->z_s;
  j_dot_k = template->y_r * template->z_r + template->y_a * template->z_a + template->y_s * template->z_s;
  if(i_dot_j != 0.0 || i_dot_k != 0.0 || i_dot_k != 0.0)
  {
    printf("warning: output volume axes are not orthogonal\n");
    printf("i_ras = (%g, %g, %g)\n", template->x_r, template->x_a, template->x_s);
    printf("j_ras = (%g, %g, %g)\n", template->y_r, template->y_a, template->y_s);
    printf("k_ras = (%g, %g, %g)\n", template->z_r, template->z_a, template->z_s);
  }
  if(out_data_type >= 0)
    template->type = out_data_type;

  /* ----- catch the template info flag ----- */
  if(template_info_flag)
  {
    printf("template structure:\n");
    MRIdump(template, stdout);
  }

  /* ----- exit here if read only is desired ----- */
  if(read_only_flag)  exit(0);

  /* ----- change type if necessary ----- */
  if(mri->type != template->type)
  {
    printf("changing data type from %d to %d (noscale = %d)...\n",
	   mri->type,template->type,no_scale_flag);
    mri2 = MRIchangeType(mri, template->type, 0.0, 0.999, no_scale_flag);
    if(mri2 == NULL) {
      printf("ERROR: MRIchangeType\n");
      exit(1);
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- reslice if necessary ----- */
  if(mri->xsize != template->xsize || 
     mri->ysize != template->ysize || 
     mri->zsize != template->zsize ||
     mri->width != template->width || 
     mri->height != template->height || 
     mri->depth != template->depth ||
     mri->x_r != template->x_r || 
     mri->x_a != template->x_a || 
     mri->x_s != template->x_s ||
     mri->y_r != template->y_r || 
     mri->y_a != template->y_a || 
     mri->y_s != template->y_s ||
     mri->z_r != template->z_r || 
     mri->z_a != template->z_a || 
     mri->z_s != template->z_s ||
     mri->c_r != template->c_r || 
     mri->c_a != template->c_a || 
     mri->c_s != template->c_s)
  {
    printf("Reslicing using ");
    switch(resample_type_val){
    case RESAMPLE_INTERPOLATE: printf("trilinear interpolation \n"); break;
    case RESAMPLE_NEAREST:     printf("nearest \n"); break;
    case RESAMPLE_SINC:        printf("sinc \n"); break;
    case RESAMPLE_WEIGHTED:    printf("weighted \n"); break;
    }
    mri2 = MRIresample(mri, template, resample_type_val);
    if(mri2 == NULL) exit(1);
    MRIfree(&mri);
    mri = mri2;
  }


  /* ----- invert contrast if necessary ----- */
  if(invert_val >= 0)
  {
    printf("inverting contrast...\n");
    mri2 = MRIinvertContrast(mri, NULL, invert_val);
    if(mri2 == NULL)
      exit(1);
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- reorder if necessary ----- */
  if(reorder_flag){
    printf("reordering axes...\n");
    mri2 = MRIreorder(mri, NULL, reorder_vals[0], reorder_vals[1], 
          reorder_vals[2]);
    if(mri2 == NULL){
      fprintf(stderr, "error reordering axes\n");
      exit(1);
    }
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- store the gdf file stem ----- */
  strcpy(mri->gdf_image_stem, gdf_image_stem);

  /* ----- catch the out info flag ----- */
  if(out_info_flag){
    printf("output structure:\n");
    MRIdump(mri, stdout);
  }

  /* ----- catch the out matrix flag ----- */
  if(out_matrix_flag)
  {
    MATRIX *i_to_r;
    i_to_r = extract_i_to_r(mri);
    if(i_to_r != NULL)
    {
      printf("output ijk -> ras:\n");
      MatrixPrint(stdout, i_to_r);
      MatrixFree(&i_to_r);
    }
    else
      printf("error getting output matrix\n");
  }

  /* ----- catch the out stats flag ----- */
  if(out_stats_flag) MRIprintStats(mri, stdout);

  /*------ Finally, write the output -----*/
  if(!no_write_flag)
  {
    printf("writing to %s...\n", out_name);
    if(force_out_type_flag){
      if(MRIwriteType(mri, out_name, out_volume_type) != NO_ERROR){
	printf("ERROR: writing %s as %d\n",out_name,out_volume_type);
        exit(1);
      }
    }
    else{
      if(MRIwrite(mri, out_name) != NO_ERROR){
	printf("ERROR: writing %s\n",out_name);
        exit(1);
      }
    }
  }

  exit(0);

} /* end main() */
/*----------------------------------------------------------------------*/

void get_ints(int argc, char *argv[], int *pos, int *vals, int nvals)
{

  char *ep;
  int i;

  if(*pos + nvals >= argc)
  {
    fprintf(stderr, "\n%s: argument %s expects %d integers; only %d arguments after flag\n", Progname, argv[*pos], nvals, argc - *pos - 1);
    usage_message(stdout);
    exit(1);
  }

  for(i = 0;i < nvals;i++)
  {
    if(argv[*pos+i+1][0] == '\0')
    {
      fprintf(stderr, "\n%s: argument to %s flag is null\n", Progname, argv[*pos]);
      usage_message(stdout);
      exit(1);
    }

    vals[i] = (int)strtol(argv[*pos+i+1], &ep, 10);

    if(*ep != '\0')
    {
      fprintf(stderr, "\n%s: error converting \"%s\" to an integer for %s flag\n", Progname, argv[*pos+i+1], argv[*pos]);
      usage_message(stdout);
      exit(1);
    }

  }

  *pos += nvals;

} /* end get_ints() */

void get_floats(int argc, char *argv[], int *pos, float *vals, int nvals)
{

  char *ep;
  int i;

  if(*pos + nvals >= argc)
  {
    fprintf(stderr, "\n%s: argument %s expects %d floats; only %d arguments after flag\n", Progname, argv[*pos], nvals, argc - *pos - 1);
    usage_message(stdout);
    exit(1);
  }

  for(i = 0;i < nvals;i++)
  {
    if(argv[*pos+i+1][0] == '\0')
    {
      fprintf(stderr, "\n%s: argument to %s flag is null\n", Progname, argv[*pos]);
      usage_message(stdout);
      exit(1);
    }

    vals[i] = (float)strtod(argv[*pos+i+1], &ep);

    if(*ep != '\0')
    {
      fprintf(stderr, "\n%s: error converting \"%s\" to an float for %s flag\n", Progname, argv[*pos+i+1], argv[*pos]);
      usage_message(stdout);
      exit(1);
    }

  }

  *pos += nvals;

} /* end get_floats() */

void get_string(int argc, char *argv[], int *pos, char *val)
{

  if(*pos + 1 >= argc)
  {
    fprintf(stderr, "\n%s: argument %s expects an extra argument; none found\n", Progname, argv[*pos]);
    usage_message(stdout);
    exit(1);
  }  

  strcpy(val, argv[*pos+1]);

  (*pos)++;

} /* end get_string() */

void usage_message(FILE *stream)
{

  fprintf(stream, "\n");
  fprintf(stream, "type %s -u for usage\n", Progname);
  fprintf(stream, "\n");

} /* end usage_message() */

void usage(FILE *stream)
{

  fprintf(stream, "\n");
  fprintf(stream, "usage: %s [options] <in volume> <out volume>\n", Progname);
  fprintf(stream, "\n");
  fprintf(stream, "options are:\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -ro, --read_only\n");
  fprintf(stream, "  -nw, --no_write\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -ii, --in_info\n");
  fprintf(stream, "  -oi, --out_info\n");
  fprintf(stream, "  -is, --in_stats\n");
  fprintf(stream, "  -os, --out_stats\n");
  fprintf(stream, "  -im, --in_matrix\n");
  fprintf(stream, "  -om, --out_matrix\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -iis, --in_i_size <size>\n");
  fprintf(stream, "  -ijs, --in_j_size <size>\n");
  fprintf(stream, "  -iks, --in_k_size <size>\n");
  fprintf(stream, "  --force_ras_good : use default when orientation info absent\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -iid, --in_i_direction <R direction> <A direction> <S direction>\n");
  fprintf(stream, "  -ijd, --in_j_direction <R direction> <A direction> <S direction>\n");
  fprintf(stream, "  -ikd, --in_k_direction <R direction> <A direction> <S direction>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -ic, --in_center <R coordinate> <A coordinate> <S coordinate>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -oni, -oic, --out_i_count <count>\n");
  fprintf(stream, "  -onj, -ojc, --out_j_count <count>\n");
  fprintf(stream, "  -onk, -okc, --out_k_count <count>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -ois, --out_i_size <size>\n");
  fprintf(stream, "  -ojs, --out_j_size <size>\n");
  fprintf(stream, "  -oks, --out_k_size <size>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -oid, --out_i_direction <R direction> <A direction> <S direction>\n");
  fprintf(stream, "  -ojd, --out_j_direction <R direction> <A direction> <S direction>\n");
  fprintf(stream, "  -okd, --out_k_direction <R direction> <A direction> <S direction>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -oc, --out_center <R coordinate> <A coordinate> <S coordinate>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -odt, --out_data_type <uchar|short|int|float>\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -rt, --resample_type <interpolate|weighted|nearest|sinc> (default is interpolate)\n");
  fprintf(stream, "\n");
  fprintf(stream, "  --no_scale flag <-ns>: 1 = dont rescale values for COR\n");
  fprintf(stream, "\n");
  fprintf(stream, "\n");
  fprintf(stream, "  -tr TR : TR in seconds\n");
  fprintf(stream, "\n");
  fprintf(stream, "  --unwarp_gradient_nonlinearity \n"
                  "      <sonata | allegra | GE> \n"
                  "      <fullUnwarp | through-plane> \n"
                  "      <JacobianCorrection | noJacobianCorrection> \n"
                  "      <linear | sinc> \n"
            "      <sincInterpHW>  \n"); 
  fprintf(stream, "\n");
  fprintf(stream, "--apply_transform xfmfile (-T or -at)\n");
  fprintf(stream, "--apply_inverse_transform xfmfile (-ait)\n");

  fprintf(stream, 

  "\n"
  "SPECIFYING THE INPUT AND OUTPUT FILE TYPES\n"
  "\n"
  "The file type can be specified in two ways. First, mri_convert will try \n"
  "to figure it out on its own from the format of the file name (eg, files that\n"
  "end in .img are assumed to be in spm analyze format). Second, the user can \n"
  "explicity set the type of file using --in_type and/or --out_type.\n"
  "\n"
  "Legal values for --in_type (-it) and --out_type (-ot) are:\n"
  "\n"
  "  cor           - MGH-NMR COR format\n"
  "  minc          - MNI's Medical Imaging NetCDF format (output may not work)\n"
  "  analyze       - 3D analyze (same as spm)\n"
  "  analyze4d     - 4D analyze \n"
  "  spm           - SPM Analyze format (same as analyze and analyze3d)\n"
  "  ge            - GE Genesis format (input only)\n"
  "  gelx          - GE LX (input only)\n"
  "  lx            - same as gelx\n"
  "  siemens       - Siemens IMA (input only)\n"
  "  dicom         - generic DICOM Format (input only)\n"
  "  siemens_dicom - Siemens DICOM Format (input only)\n"
  "  afni          - AFNI format\n"
  "  brik          - same as afni\n"
  "  bshort        - MGH-NMR bshort format\n"
  "  bfloat        - MGH-NMR bfloat format\n"
  "  sdt           - Varian (?)\n"
  "  outline       - MGH-NMR Outline format\n"
  "  otl           - same as outline\n"
  "  gdf           - GDF volume (requires image stem for output; use -gis)\n"
  "\n"
  "CONVERTING TO SPM-ANALYZE FORMAT \n"
  "\n"
  "Converting to SPM-Analyze format can be done in two ways, depeding upon\n"
  "whether a single frame or multiple frames are desired. For a single frame,\n"
  "simply specify the output file name with a .img extension, and mri_convert \n"
  "will save the first frame into the file.  For multiple frames, specify the \n"
  "base as the output file name and add --out_type spm. This will save each \n"
  "frame as baseXXX.img where XXX is the three-digit, zero-padded frame number.\n"
  "Frame numbers begin at one. By default, the width the of zero padding is 3.\n"
  "This can be controlled with --in_nspmzeropad N where N is the new width.\n"
  "\n"  );

  printf("\n");
  printf("Other options\n");
  printf("\n");
  printf("  -r, --reorder olddim1 olddim2 olddim3\n");
  printf("\n");
  printf("  Reorders axes such that olddim1 is the new column dimension,\n");
  printf("  olddim2 is the new row dimension, olddim3 is the new slice \n");
  printf("  dimension. Example: 2 1 3 will swap rows and cols.\n");
  printf("\n");
  printf("  --invert_contrast threshold\n");
  printf("\n");
  printf("  All voxels in volume greater than threshold are replaced\n");
  printf("  with 255-value. Only makes sense for 8 bit images.\n");
  printf("  Only operates on the first frame.\n");
  printf("\n");
  printf("  -i, --input_volume\n");
  printf("  -o, --output_volume\n");
  printf("  -nc, --no_conform\n");
  printf("  -po, --parse_only\n");
  printf("  -is, --in_stats\n");
  printf("  -os, --out_stats\n");
  printf("  -ro, --read_only\n");
  printf("  -nw, --no_write\n");
  printf("  -sn, --subject_name\n");
  printf("  -rl, --reslice_like\n");
  printf("  -f,  --frame\n");
  printf("  -il, --in_like\n");
  printf("  -roi\n");
  printf("  -fp, --fill_parcellation\n");
  printf("  -sp, --smooth_parcellation\n");
  printf("  -zo, --zero_outlines\n");
  printf("  -cf, --color_file\n");
  printf("  -nt, --no_translate\n");
  printf("  --status (status file for DICOM conversion)\n");
  printf("  --sdcmlist (list of DICOM files for conversion)\n");
  printf("  -ti, --template_info : dump info about template\n");
  printf("  -gis <gdf image file stem>\n");
  printf("\n");
  printf("Notes: \n");
  printf("\n");
  printf("If the user specifies any of the direction cosines, the ras_good_flag is set.\n");
  printf("\n");

} /* end usage() */

float findMinSize(MRI *mri)
{
  float xsize, ysize, zsize;
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
    return (ysize > zsize) ? zsize : ysize;
  else
    return (zsize > xsize) ? xsize : zsize;
  
}
/* EOF */
