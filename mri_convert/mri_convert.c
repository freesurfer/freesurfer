#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "mri.h"
#include "error.h"
#include "mri_identify.h"
#include "utils.h"

void get_ints(int argc, char *argv[], int *pos, int *vals, int nvals);
void get_floats(int argc, char *argv[], int *pos, float *vals, int nvals);
void get_string(int argc, char *argv[], int *pos, char *val);
void usage_message(FILE *stream);
void usage(FILE *stream);

char *Progname;

int main(int argc, char *argv[])
{

  MRI *mri, *mri2, *template;
  int i;
  int reorder_vals[3];
  float invert_val;
  int in_info_flag, out_info_flag;
  int template_info_flag;
  int conform_flag;
  int parse_only_flag;
  int reorder_flag;
  int in_stats_flag, out_stats_flag;
  int read_only_flag, no_write_flag;
  char in_name[STRLEN], out_name[STRLEN];
  int out_volume_type;
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
  float magnitude;
  float i_dot_j, i_dot_k, j_dot_k;
  float in_center[3], out_center[3];
  int in_center_flag, out_center_flag;
  int out_data_type;
  char out_data_type_string[STRLEN];
  int out_n_i, out_n_j, out_n_k;
  int out_n_i_flag, out_n_j_flag, out_n_k_flag;

  /* ----- keep the compiler quiet ----- */
  mri2 = NULL;

  /* ----- get the program name ----- */
  Progname = strrchr(argv[0], '/');
  Progname = (Progname == NULL ? argv[0] : Progname + 1);

  /* ----- pass the command line to mriio ----- */
  mriio_command_line(argc, argv);

  /* ----- catch no arguments here ----- */
  if(argc == 1)
  {
    usage(stderr);
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

  for(i = 1;i < argc;i++)
  {
    if(strcmp(argv[i], "-version2") == 0)
      exit(97);
    if(strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reorder") == 0)
      get_ints(argc, argv, &i, reorder_vals, 3);
    else if(strcmp(argv[i], "--invert_contrast") == 0)
      get_floats(argc, argv, &i, &invert_val, 1);
    else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input_volume") == 0)
      get_string(argc, argv, &i, in_name);
    else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output_volume") == 0)
      get_string(argc, argv, &i, out_name);
    else if(strcmp(argv[i], "-nc") == 0 || strcmp(argv[i], "--no_conform") == 0)
      conform_flag = FALSE;
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
    else if(strcmp(argv[i], "-iis") == 0 || strcmp(argv[i], "--in_i_size") == 0)
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
        usage_message(stderr);
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
        usage_message(stderr);
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
        usage_message(stderr);
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
        usage_message(stderr);
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
        usage_message(stderr);
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
        usage_message(stderr);
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
    else if(strcmp(argv[i], "-odt") == 0 || strcmp(argv[i], "--out_data_type") == 0)
    {
      get_string(argc, argv, &i, out_data_type_string);
      if(strcmp(StrLower(out_data_type_string), "uchar") == 0)
        out_data_type = MRI_UCHAR;
      if(strcmp(StrLower(out_data_type_string), "short") == 0)
        out_data_type = MRI_SHORT;
      if(strcmp(StrLower(out_data_type_string), "int") == 0)
        out_data_type = MRI_INT;
      if(strcmp(StrLower(out_data_type_string), "float") == 0)
        out_data_type = MRI_FLOAT;
      else
      {
        fprintf(stderr, "\n%s: unknown data type \"%s\"\n", Progname, argv[i]);
        usage_message(stderr);
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
      else
      {
        fprintf(stderr, "\n%s: unknown resample type \"%s\"\n", Progname, argv[i]);
        usage_message(stderr);
        exit(1);
      }
    }
    else if(strcmp(argv[i], "-u") == 0 || strcmp(argv[i], "--usage") == 0)
    {
      usage(stdout);
      exit(0);
    }
    else
    {
      if(argv[i][0] == '-')
      {
        fprintf(stderr, "\n%s: unknown flag \"%s\"\n", Progname, argv[i]);
        usage_message(stderr);
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
          usage_message(stderr);
          exit(1);
        }
      }

    }
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
    usage_message(stderr);
    exit(1);
  }

  /* ----- catch missing input or output volume name ----- */
  if(in_name[0] == '\0')
  {
    fprintf(stderr, "\n%s: missing input volume name\n", Progname);
    usage_message(stderr);
    exit(1);
  }

  if(out_name[0] == '\0' && !(read_only_flag || no_write_flag))
  {
    fprintf(stderr, "\n%s: missing output volume name\n", Progname);
    usage_message(stderr);
    exit(1);
  }

  /* ----- catch the parse-only flag ----- */
  if(parse_only_flag)
  {

    printf("input volume name: %s\n", in_name);
    printf("output volume name: %s\n", out_name);
    printf("parse_only_flag = %d\n", parse_only_flag);
    printf("conform_flag = %d\n", conform_flag);
    printf("in_info_flag = %d\n", in_info_flag);
    printf("out_info_flag = %d\n", out_info_flag);

    if(invert_val >= 0)
      printf("inversion, value is %g\n", invert_val);

    if(reorder_flag)
      printf("reordering, values are %d %d %d\n", reorder_vals[0], reorder_vals[1], reorder_vals[2]);

    exit(0);

  }

  /* ----- get the type of the output ----- */
  if(!read_only_flag && !no_write_flag)
  {
    out_volume_type = mri_identify(out_name);
    if(out_volume_type < 0)
    {
      fprintf(stderr, "%s: can't determine type of output volume\n", Progname);
      exit(1);
    }
  }

  /* ----- read the volume ----- */
  printf("reading from %s...\n", in_name);
  if(read_only_flag && in_info_flag && !in_stats_flag)
    mri = MRIreadInfo(in_name);
  else
    mri = MRIread(in_name);
  if(mri == NULL)
    exit(1);

  /* ----- check for ras good flag -- warn if it's not set ----- */
  if(mri->ras_good_flag == 0)
  {
    printf("warning: volume may be incorrectly oriented\n");
    if(mri_identify(in_name) == MRI_CORONAL_SLICE_DIRECTORY)
      printf("(but as a COR- volume, it should be okay...)\n");
  }

  /* ----- apply command-line parameters ----- */
  if(in_i_size_flag)
    mri->xsize = in_i_size;
  if(in_j_size_flag)
    mri->ysize = in_j_size;
  if(in_k_size_flag)
    mri->zsize = in_k_size;
  if(in_i_direction_flag)
  {
    mri->x_r = in_i_directions[0];
    mri->x_a = in_i_directions[1];
    mri->x_s = in_i_directions[2];
  }
  if(in_j_direction_flag)
  {
    mri->y_r = in_j_directions[0];
    mri->y_a = in_j_directions[1];
    mri->y_s = in_j_directions[2];
  }
  if(in_k_direction_flag)
  {
    mri->z_r = in_k_directions[0];
    mri->z_a = in_k_directions[1];
    mri->z_s = in_k_directions[2];
  }
  if(in_center_flag)
  {
    mri->c_r = in_center[0];
    mri->c_a = in_center[1];
    mri->c_s = in_center[2];
  }

  /* ----- give a warning for non-orthogonal directions ----- */
  i_dot_j = mri->x_r * mri->y_r + mri->x_a * mri->y_a + mri->x_s * mri->y_s;
  i_dot_k = mri->x_r * mri->z_r + mri->x_a * mri->z_a + mri->x_s * mri->z_s;
  j_dot_k = mri->y_r * mri->z_r + mri->y_a * mri->z_a + mri->y_s * mri->z_s;
  if(i_dot_j != 0.0 || i_dot_k != 0.0 || i_dot_k != 0.0)
  {
    printf("warning: input volume axes are not orthogonal\n");
    printf("i_ras = (%g, %g, %g)\n", mri->x_r, mri->x_a, mri->x_s);
    printf("j_ras = (%g, %g, %g)\n", mri->y_r, mri->y_a, mri->y_s);
    printf("k_ras = (%g, %g, %g)\n", mri->z_r, mri->z_a, mri->z_s);
  }

  /* ----- catch the in info flag ----- */
  if(in_info_flag)
  {
    printf("input structure:\n");
    MRIdump(mri, stdout);
  }

  /* ----- catch the in stats flag ----- */
  if(in_stats_flag)
    MRIprintStats(mri, stdout);

  if(read_only_flag)
    exit(0);

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
      template->thick = 1.0;
      template->ps = 1.0;
      template->xsize = template->ysize = template->zsize = 1.0;
      template->xstart = template->ystart = template->zstart = -128.0;
      template->xend = template->yend = template->zend = 128.0;
      template->x_r = -1.0;  template->x_a =  0.0;  template->x_s =  0.0;
      template->y_r =  0.0;  template->y_a =  0.0;  template->y_s = -1.0;
      template->z_r =  0.0;  template->z_a =  1.0;  template->z_s =  0.0;
    }
    else
      printf("the output volume is not a COR- directory; the --no_conform (-nc) argument is not needed\n");
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

  /* ----- change type if necessary ----- */
  if(mri->type != template->type)
  {
    printf("changing data type...\n");
    mri2 = MRIchangeType(mri, template->type, 0.0, 0.999);
    if(mri2 == NULL)
      exit(1);
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- catch the template info flag ----- */
  if(template_info_flag)
  {
    printf("template structure:\n");
    MRIdump(template, stdout);
  }

  /* ----- reslice if necessary ----- */
  if(mri->xsize != template->xsize || mri->ysize != template->ysize || mri->zsize != template->zsize ||
     mri->width != template->width || mri->height != template->height || mri->depth != template->depth ||
     mri->x_r != template->x_r || mri->x_a != template->x_a || mri->x_s != template->x_s ||
     mri->y_r != template->y_r || mri->y_a != template->y_a || mri->y_s != template->y_s ||
     mri->z_r != template->z_r || mri->z_a != template->z_a || mri->z_s != template->z_s ||
     mri->c_r != template->c_r || mri->c_a != template->c_a || mri->c_s != template->c_s)
  {
    printf("reslicing...\n");
    mri2 = MRIresample(mri, template, resample_type_val);
    if(mri2 == NULL)
      exit(1);
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
  if(reorder_flag)
  {
    printf("reordering axes...\n");
    MRIreorder(mri, mri2, reorder_vals[0], reorder_vals[1], reorder_vals[2]);
    if(mri2 == NULL)
      exit(1);
    MRIfree(&mri);
    mri = mri2;
  }

  /* ----- catch the out info flag ----- */
  if(out_info_flag)
  {
    printf("output structure:\n");
    MRIdump(mri, stdout);
  }

  /* ----- catch the out stats flag ----- */
  if(out_stats_flag)
    MRIprintStats(mri, stdout);

  if(!no_write_flag)
  {
    printf("writing to %s...\n", out_name);
    MRIwrite(mri, out_name);
  }

  exit(0);

} /* end main() */

void get_ints(int argc, char *argv[], int *pos, int *vals, int nvals)
{

  char *ep;
  int i;

  if(*pos + nvals >= argc)
  {
    fprintf(stderr, "\n%s: argument %s expects %d integers; only %d arguments after flag\n", Progname, argv[*pos], nvals, argc - *pos - 1);
    usage_message(stderr);
    exit(1);
  }

  for(i = 0;i < nvals;i++)
  {
    if(argv[*pos+i+1][0] == '\0')
    {
      fprintf(stderr, "\n%s: argument to %s flag is null\n", Progname, argv[*pos]);
      usage_message(stderr);
      exit(1);
    }

    vals[i] = (int)strtol(argv[*pos+i+1], &ep, 10);

    if(*ep != '\0')
    {
      fprintf(stderr, "\n%s: error converting \"%s\" to an integer for %s flag\n", Progname, argv[*pos+i+1], argv[*pos]);
      usage_message(stderr);
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
    fprintf(stderr, "\n%s: argument %s expects %d integers; only %d arguments after flag\n", Progname, argv[*pos], nvals, argc - *pos - 1);
    usage_message(stderr);
    exit(1);
  }

  for(i = 0;i < nvals;i++)
  {
    if(argv[*pos+i+1][0] == '\0')
    {
      fprintf(stderr, "\n%s: argument to %s flag is null\n", Progname, argv[*pos]);
      usage_message(stderr);
      exit(1);
    }

    vals[i] = (float)strtod(argv[*pos+i+1], &ep);

    if(*ep != '\0')
    {
      fprintf(stderr, "\n%s: error converting \"%s\" to an integer for %s flag\n", Progname, argv[*pos+i+1], argv[*pos]);
      usage_message(stderr);
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
    usage_message(stderr);
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
  fprintf(stream, "\n");
  fprintf(stream, "  -iis, --in_i_size <size>\n");
  fprintf(stream, "  -ijs, --in_j_size <size>\n");
  fprintf(stream, "  -iks, --in_k_size <size>\n");
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

} /* end usage() */

/* EOF */
