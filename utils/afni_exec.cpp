/*
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <fmt/format.h>

#include "error.h"
#include "machine.h"
#include "mghendian.h"
#include "mri.h"
#include "utils.h"

#include "utils/afni.hpp"

/* ----- flags for keeping track of what we've gotten from the header ----- */
#define AFNI_ALL_REQUIRED 0x0000003f

#define ORIENT_SPECIFIC_FLAG    0x00000001
#define BRICK_TYPES_FLAG        0x00000002
#define DATASET_DIMENSIONS_FLAG 0x00000004
#define DELTA_FLAG              0x00000008
#define ORIGIN_FLAG             0x00000010
#define BYTEORDER_STRING_FLAG   0x00000020

namespace fs::utils {

auto afni_orientations() -> std::vector<std::vector<float>> & {
  static std::vector<std::vector<float>> afni_orientations = {
      {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
      {0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}};
  return afni_orientations;
}

void AFinit(AF &pAF) {
  pAF.dataset_rank[0]       = 0;
  pAF.dataset_rank[1]       = 0;
  pAF.dataset_dimensions[0] = 0;
  pAF.dataset_dimensions[1] = 0;
  pAF.dataset_dimensions[2] = 0;
  strcpy(pAF.typestring, "");
  pAF.scene_data[0]      = 0;
  pAF.scene_data[1]      = 0;
  pAF.scene_data[2]      = 0;
  pAF.orient_specific[0] = 0;
  pAF.orient_specific[1] = 0;
  pAF.orient_specific[2] = 0;
  pAF.origin[0]          = 0.;
  pAF.origin[1]          = 0.;
  pAF.origin[2]          = 0.;
  pAF.delta[0]           = 0.;
  pAF.delta[1]           = 0.;
  pAF.delta[2]           = 0.;
  pAF.numchars           = 0;
  pAF.idcode_string      = nullptr;
  pAF.numstats           = 0;
  pAF.brick_stats        = nullptr;
  pAF.numtypes           = 0;
  pAF.brick_types        = nullptr;
  pAF.numfacs            = 0;
  pAF.brick_float_facs   = nullptr;
}

void AFclean(AF *pAF) {
  free(pAF->idcode_string);
  pAF->idcode_string = nullptr;
  free(pAF->brick_stats);
  pAF->brick_stats = nullptr;
  free(pAF->brick_types);
  pAF->brick_types = nullptr;
  free(pAF->brick_float_facs);
  pAF->brick_float_facs = nullptr;
}

void printAFNIHeader(AF &pAF) {
  std::vector<std::string> types = {"byte", "short", "float", "complex"};

  int i = 0;
  fmt::print(
      "AFNI Header Information ============================================\n");
  fmt::print("DATASET_RANK      : spatial dims {}, sub-bricks {}\n",
             pAF.dataset_rank[0], pAF.dataset_rank[1]);
  fmt::print("DATASET_DIMENSIONS: ({}, {}, {})\n", pAF.dataset_dimensions[0],
             pAF.dataset_dimensions[1], pAF.dataset_dimensions[2]);
  fmt::print("TYPESTRING        : {}\n", pAF.typestring);
  fmt::print("SCENE_DATA        : view type {}, func type {}, verify {}\n",
             pAF.scene_data[0], pAF.scene_data[1], pAF.scene_data[2]);
  fmt::print("ORIGIN            : ({:f}, {:f}, {:f})\n", pAF.origin[0],
             pAF.origin[1], pAF.origin[2]);
  fmt::print("DELTA             : ({:f}, {:f}, {:f})\n", pAF.delta[0],
             pAF.delta[1], pAF.delta[2]);
  fmt::print("IDCODE_STRING     : {:s}\n", pAF.idcode_string);
  fmt::print("BYTEORDER_STRING  : {:s}\n", pAF.byteorder_string);
  for (i = 0; i < pAF.numstats; i = i + 2) {
    fmt::print("BRICK_STATS       : min {:f}\t max {:f}\n", pAF.brick_stats[i],
               pAF.brick_stats[i + 1]);
  }
  for (i = 0; i < pAF.numtypes; ++i) {
    fmt::print("BRICK_TYPES       : {:s}\n", types[pAF.brick_types[i]]);
  }
  for (i = 0; i < pAF.numfacs; ++i) {
    fmt::print("BRICK_FLOAT_FACS  : {:f}\n", pAF.brick_float_facs[i]);
  }
  fmt::print(
      "====================================================================\n");
}

static auto get_afni_int(FILE *fp, int count, char *name) -> int * {
  int * buf = nullptr;
  int   i   = 0;
  char  line[STRLEN];
  char *c          = nullptr;
  char *e          = nullptr;
  char  blank_flag = 0;

  buf = static_cast<int *>(malloc(count * sizeof(int)));

  for (i = 0; i < count;) {
    if ((fgets(line, STRLEN, fp) == nullptr) && (ferror(fp) != 0)) {
      ErrorPrintf(ERROR_BADFILE, "failed reading file");
    }

    blank_flag = 1;
    for (c = line; *c != '\0'; c++) {
      if (isspace(*c) == 0) {
        blank_flag = 0;
      }
    }

    if (feof(fp) != 0) {
      free(buf);
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM,
                         "get_afni_int(): hit EOF while reading %d ints for %s",
                         count, name));
    }

    if (blank_flag != 0) {
      free(buf);
      errno = 0;
      ErrorReturn(
          NULL,
          (ERROR_BADPARM,
           "get_afni_int(): hit a blank line while reading %d ints for %s",
           count, name));
    }

    if (strncmp(line, "type", 4) == 0) {
      free(buf);
      errno = 0;
      ErrorReturn(
          NULL, (ERROR_BADPARM,
                 "get_afni_int(): hit a type line while reading %d ints for %s",
                 count, name));
    }

    for (c = line; *c != '\0';) {
      for (; (isspace(*c) != 0) && *c != '\0'; c++) {
        ;
      }
      if (*c != '\0') {
        buf[i] = strtol(c, &e, 10);
        c      = e;
        i++;
      }
    }
  }

  return (buf);

} /* end get_afni_int() */

static auto get_afni_float(FILE *fp, int count, char *name) -> float * {
  float *buf = nullptr;
  int    i   = 0;
  char   line[STRLEN];
  char * c          = nullptr;
  char * e          = nullptr;
  char   blank_flag = 0;

  buf = static_cast<float *>(malloc(count * sizeof(float)));

  for (i = 0; i < count;) {
    if ((fgets(line, STRLEN, fp) == nullptr) && (ferror(fp) != 0)) {
      ErrorPrintf(ERROR_BADFILE, "failed reading file");
    }

    blank_flag = 1;
    for (c = line; *c != '\0'; c++) {
      if (isspace(*c) == 0) {
        blank_flag = 0;
      }
    }

    if (feof(fp) != 0) {
      free(buf);
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_BADPARM,
                   "get_afni_float(): hit EOF while reading %d floats for %s",
                   count, name));
    }

    if (blank_flag != 0) {
      free(buf);
      errno = 0;
      ErrorReturn(
          NULL,
          (ERROR_BADPARM,
           "get_afni_float(): hit a blank line while reading %d floats for %s",
           count, name));
    }

    if (strncmp(line, "type", 4) == 0) {
      free(buf);
      errno = 0;
      ErrorReturn(
          NULL,
          (ERROR_BADPARM,
           "get_afni_float(): hit a type line while reading %d floats for %s",
           count, name));
    }

    for (c = line; *c != '\0';) {
      for (; (isspace(*c) != 0) && *c != '\0'; c++) {
        ;
      }
      if (*c != '\0') {
        buf[i] = static_cast<float>(strtod(c, &e));
        c      = e;
        i++;
      }
    }
  }

  return (buf);

} /* end get_afni_float() */

static auto get_afni_string(FILE *fp, int count, char *name) -> char * {
  char *buf = nullptr;
  int   i   = 0;
  char  c   = 0;

  buf = static_cast<char *>(malloc(count + 1));

  c = fgetc(fp);

  if (c != '\'') {
    free(buf);
    errno = 0;
    ErrorReturn(
        NULL,
        (ERROR_BADPARM,
         "get_afni_string(): afni header string %s does not start with \"'\"",
         name));
  }

  for (i = 0; i < count; i++) {
    if (feof(fp) != 0) {
      free(buf);
      errno = 0;
      ErrorReturn(NULL, (ERROR_BADPARM,
                         "get_afni_string(): end of file reached at %d of %d "
                         "bytes of string %s",
                         i + 1, count, name));
    }

    c = fgetc(fp);

    if (i == count - 1 && c != '~') {
      errno = 0;
      ErrorPrintf(ERROR_BADPARM, "warning: string %s does not end with \"~\"",
                  name);
    }

    buf[i] = (c == '~' ? '\0' : c);
  }

  buf[count] = '\0';

  for (c = fgetc(fp); c != '\n' && (feof(fp) == 0); c = fgetc(fp)) {
    ;
  }

  return (buf);

} /* end get_afni_string() */

auto readAFNIHeader(FILE *fp, AF *pAF) -> int {
  char   line[STRLEN];
  char   line2[STRLEN];
  int    i = 0;
  int    j = 0;
  char   type[STRLEN];
  char   name[STRLEN];
  int    count  = 0;
  char * s      = nullptr;
  float *f      = nullptr;
  int *  ip     = nullptr;
  long   gotten = 0;

  int val = 0;

  fseek(fp, 0, SEEK_SET);

  // AFNI header attribute spec https://github.com/afni/afni/blob/master/doc/README/README.attributes
  while (true) // !feof(fp))
  {
    if ((fgets(line, STRLEN, fp) == nullptr) && (ferror(fp) != 0)) {
      ErrorPrintf(ERROR_BADFILE, "failed reading file");
    }
    if (feof(fp) != 0) { // wow.  we read too many.  get out
      break;
    }

    if (feof(fp) == 0) {
      i = -1;
      j = -1;
      do {
        i++;
        j++;
        for (; (isspace(line[j]) != 0) && line[j] != '\0'; j++) {
          ;
        }
        line2[i] = line[j];
      } while (line[j] != '\0');

      if (strlen(line2) > 0) {
        if (line2[strlen(line2) - 1] == '\n') {
          line2[strlen(line2) - 1] = '\0';
        }
      }

      if (strncmp(line2, "type=", 5) == 0) {
        strcpy(type, &line2[5]);
      }
      if (strncmp(line2, "name=", 5) == 0) {
        strcpy(name, &line2[5]);
      }
      if (strncmp(line2, "count=", 5) == 0) {
        count = atoi(&line2[6]);
        s     = nullptr;
        f     = nullptr;
        ip    = nullptr;
        // depending on the attribute, get the array
        if (strncmp(type, "string-attribute", 16) == 0) {
          s = get_afni_string(fp, count, name);
          if (s == nullptr) {
            return (0);
          }
        } else if (strncmp(type, "float-attribute", 15) == 0) {
          f = get_afni_float(fp, count, name);
          if (f == nullptr) {
            return (0);
          }
        } else if (strncmp(type, "integer-attribute", 17) == 0) {
          ip = get_afni_int(fp, count, name);
          if (ip == nullptr) {
            return (0);
          }
        } else {
          errno = 0;
          ErrorReturn(
              0, (ERROR_BADPARM, "read_afni_header(): unknown type %s", type));
        }
        // now compare name with the attribute
        // mandatory attributes
        if (strcmp(name, "DATASET_RANK") == 0) {
          pAF->dataset_rank[0] = ip[0];
          pAF->dataset_rank[1] = ip[1];
        } else if (strcmp(name, "DATASET_DIMENSIONS") == 0) {
          // check errors for getting non-integer
          if (strncmp(type, "integer-attribute", 17) != 0) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): variable %s listed as %s "
                            "(expecting integer-attribute)",
                            name, type));
          }
          // check errors for getting less than 3 counts
          if (count < 3) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): not enough variables in %s "
                            "(need 3, have %d",
                            name, count));
          }
          pAF->dataset_dimensions[0] = ip[0];
          pAF->dataset_dimensions[1] = ip[1];
          pAF->dataset_dimensions[2] = ip[2];
          // mark got the info on dimension
          gotten = gotten | DATASET_DIMENSIONS_FLAG;
        } else if (strcmp(name, "TYPESTRING") == 0) {
          strcpy(pAF->typestring, s); // s is src
        } else if (strcmp(name, "SCENE_DATA") == 0) {
          pAF->scene_data[0] = ip[0];
          pAF->scene_data[1] = ip[1];
          pAF->scene_data[2] = ip[2];
        } else if (strcmp(name, "ORIENT_SPECIFIC") == 0) {
          // error check to make sure 3 integers
          if (strncmp(type, "integer-attribute", 17) != 0) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): variable %s listed as %s "
                            "(expecting integer-attribute)",
                            name, type));
          }
          if (count < 3) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): not enough variables in %s "
                            "(need 3, have %d",
                            name, count));
          }
          if (ip[0] < 0 || ip[0] > 5 || ip[1] < 0 || ip[1] > 5 || ip[2] < 0 ||
              ip[2] > 5) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): %s variables should be 0 to "
                            "5, inclusive (are %d, %d, %d here)",
                            name, ip[0], ip[1], ip[2]));
          }
          //
          pAF->orient_specific[0] = ip[0];
          pAF->orient_specific[1] = ip[1];
          pAF->orient_specific[2] = ip[2];

          gotten = gotten | ORIENT_SPECIFIC_FLAG;
        } else if (strcmp(name, "ORIGIN") == 0) {
          if (count < 3) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): not enough variables in %s "
                            "(need 3, have %d",
                            name, count));
          }
          if (strncmp(type, "float-attribute", 15) != 0) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): variable %s listed as %s "
                            "(expecting float-attribute)",
                            name, type));
          }
          pAF->origin[0] = f[0];
          pAF->origin[1] = f[1];
          pAF->origin[2] = f[2];

          gotten = gotten | ORIGIN_FLAG;
        } else if (strcmp(name, "DELTA") == 0) {
          if (strncmp(type, "float-attribute", 15) != 0) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): variable %s listed as %s "
                            "(expecting float-attribute)",
                            name, type));
          }
          if (count < 3) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): not enough variables in %s "
                            "(need 3, have %d",
                            name, count));
          }
          pAF->delta[0] = f[0];
          pAF->delta[1] = f[1];
          pAF->delta[2] = f[2];
          gotten        = gotten | DELTA_FLAG;
        }
        ////////////////////////////////////////////////////////////////////
        //// Almost Mandatory Attributes
        ////////////////////////////////////////////////////////////////////
        else if (strcmp(name, "IDCODE_STRING") == 0) {
          pAF->idcode_string =
              static_cast<char *>(malloc(sizeof(char) * count));
          strcpy(pAF->idcode_string, s);
        } else if (strcmp(name, "BYTEORDER_STRING") == 0) {
          if (strncmp(type, "string-attribute", 16) != 0) {
            errno = 0;
            ErrorReturn(0, (ERROR_BADPARM,
                            "read_afni_header(): variable %s listed as %s "
                            "(expecting string-attribute)",
                            name, type));
          }
          strcpy(pAF->byteorder_string, s);
          gotten = gotten | BYTEORDER_STRING_FLAG;
        } else if (strcmp(name, "BRICK_STATS") == 0) {
          pAF->numstats = count;
          pAF->brick_stats =
              static_cast<float *>(malloc(sizeof(float) * count));
          for (i = 0; i < count; ++i) {
            pAF->brick_stats[i] = f[i];
          }
        } else if (strcmp(name, "BRICK_TYPES") == 0) {
          pAF->numtypes    = count;
          pAF->brick_types = static_cast<int *>(malloc(sizeof(int) * count));
          for (i = 0; i < count; ++i) {
            pAF->brick_types[i] = ip[i];
          }
          // verify all same type
          if (pAF->brick_types != nullptr) {
            val = pAF->brick_types[0];
            // we can support only one type nframs
            for (i = 1; i < count; ++i) {
              if (pAF->brick_types[i] != val)
                ErrorReturn(0, (ERROR_BADPARM,
                                "multiple data types are not supported"));
            }
          }
          gotten = gotten | BRICK_TYPES_FLAG;
        }
        // new adddition
        else if (strcmp(name, "BRICK_FLOAT_FACS") == 0) {
          pAF->numfacs = count;
          pAF->brick_float_facs =
              static_cast<float *>(malloc(sizeof(float) * count));
          for (i = 0; i < count; ++i) {
            pAF->brick_float_facs[i] = f[i];
          }
        } else /* ignore unknown variables */
        {
        }

        if (s != nullptr) {
          free(s);
        }
        if (f != nullptr) {
          free(f);
        }
        if (ip != nullptr) {
          free(ip);
        }
      }
    }
  }
  fclose(fp);

  if ((gotten & AFNI_ALL_REQUIRED) != AFNI_ALL_REQUIRED) {
    errno = 0;
    ErrorPrintf(ERROR_BADFILE, "missing fields in afni header file");

    if ((gotten & ORIENT_SPECIFIC_FLAG) == 0) {
      ErrorPrintf(ERROR_BADFILE, "  ORIENT_SPECIFIC missing");
    }
    if ((gotten & BRICK_TYPES_FLAG) == 0) {
      ErrorPrintf(ERROR_BADFILE, "  BRICK_TYPES missing");
    }
    if ((gotten & DATASET_DIMENSIONS_FLAG) == 0) {
      ErrorPrintf(ERROR_BADFILE, "  DATASET_DIMENSIONS missing");
    }
    if ((gotten & DELTA_FLAG) == 0) {
      ErrorPrintf(ERROR_BADFILE, "  DELTA missing");
    }
    if ((gotten & ORIGIN_FLAG) == 0) {
      ErrorPrintf(ERROR_BADFILE, "  ORIGIN missing");
    }
    if ((gotten & BYTEORDER_STRING_FLAG) == 0) {
      ErrorPrintf(ERROR_BADFILE, "  BYTEORDER_STRING missing");
    }

    return (0);
  }
  return 1;
}

// before calling this function initialize sM
auto findMinMaxByte(unsigned char *pS, size_t count, float *fMin, float *fMax)
    -> int {
  size_t i  = 0;
  float  fv = NAN;
  // initialize to the first element
  auto fmin = static_cast<float>(*pS);
  auto fmax = static_cast<float>(*pS);
  for (i = 0; i < count; ++i, pS++) {
    fv = static_cast<float>(*pS);
    if (fv < fmin) {
      fmin = fv;
    }
    if (fv > fmax) {
      fmax = fv;
    }
  }
  *fMin = fmin;
  *fMax = fmax;
  return 0;
}

auto findMinMaxShort(short *pS, size_t count, float *fMin, float *fMax) -> int {
  size_t i  = 0;
  float  fv = NAN;
  // initialize to the first element
  auto fmin = static_cast<float>(*pS);
  auto fmax = static_cast<float>(*pS);
  for (i = 0; i < count; ++i, pS++) {
    fv = static_cast<float>(*pS);
    if (fv < fmin) {
      fmin = fv;
    }
    if (fv > fmax) {
      fmax = fv;
    }
  }
  *fMin = fmin;
  *fMax = fmax;
  return 0;
}

auto findMinMaxFloat(float *pF, size_t count, float *fMin, float *fMax) -> int {
  size_t i = 0;
  // initialize to the first element
  float fmin = *pF;
  float fmax = *pF;
  for (i = 0; i < count; ++i, pF++) {
    if (*pF < fmin) {
      fmin = *pF;
    }
    if (*pF > fmax) {
      fmax = *pF;
    }
  }
  *fMin = fmin;
  *fMax = fmax;
  return 0;
}

/*------------------------------------------------------*/
auto afniRead(const char *fname, int read_volume) -> MRI * {
  FILE *         fp = nullptr;
  char           header_fname[STRLEN];
  char *         c                = nullptr;
  MRI *          mri              = nullptr;
  MRI *          header           = nullptr;
  int            big_endian_flag  = 0;
  long           brik_file_length = 0;
  long           nvoxels          = 0;
  int            bytes_per_voxel  = 0;
  int            i                = 0;
  int            j                = 0;
  int            k                = 0;
  int            swap_flag        = 0;
  AF             af;
  float          scaling = 1.;
  void *         pmem    = nullptr;
  float *        pf      = nullptr;
  short *        ps      = nullptr;
  unsigned char *pc      = nullptr;

  float det         = NAN;
  float xfov        = NAN;
  float yfov        = NAN;
  float zfov        = NAN;
  float fMin        = 0.;
  float fMax        = 0.;
  float flMin       = 0.;
  float flMax       = 0.;
  int   initialized = 0;
  int   frame       = 0;
  // int bytes; // set but not used

  strcpy(header_fname, fname);
  c = strrchr(header_fname, '.');

  if (c == nullptr) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "afniRead(): bad file name %s", fname));
  }

  if (strcmp(c, ".BRIK") != 0) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM, "afniRead(): bad file name %s", fname));
  }

  sprintf(c, ".HEAD");

  if ((fp = fopen(header_fname, "re")) == nullptr) {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "afniRead(): error opening file %s",
                       header_fname));
  }

  // initialize AFNI structure
  AFinit(af);

  // read header file
  if (readAFNIHeader(fp, &af) == 0) {
    return (nullptr);
  }

  printAFNIHeader(af);

  // well, we don't have time
  if (af.numtypes != 1) // should be the same as af.dataset_rank[1] = subbricks
  {
    errno = 0;
    // ErrorReturn(NULL, (ERROR_UNSUPPORTED, "afniRead(): nframes = %d (only 1 frame supported)", af.numtypes));
    fmt::print(
        "INFO: number of frames dataset_rank[1] = {:d} : numtypes = {:d} \n",
        af.dataset_rank[1], af.numtypes);
  }
  // byteorder_string : required field
  if (strcmp(af.byteorder_string, "MSB_FIRST") == 0) {
    big_endian_flag = 1;
  } else if (strcmp(af.byteorder_string, "LSB_FIRST") == 0) {
    big_endian_flag = 0;
  } else {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "read_afni_header(): unrecognized byte order string %s",
                       af.byteorder_string));
  }
  // brick_types : required field
  if (af.brick_types[0] == 2    // int
      || af.brick_types[0] > 3) // 4 = double, 5 = complex, 6 = rgb
  {
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "afniRead(): unsupported data type %d, Must be 0 "
                       "(byte), 1(short), or 3(float)",
                       af.brick_types[0]));
  }
  //////////////////////////////////////////////////////////////////////////////////
  // now we allocate space for MRI
  // dataset_dimensions : required field
  header =
      MRIallocHeader(af.dataset_dimensions[0], af.dataset_dimensions[1],
                     af.dataset_dimensions[2], MRI_UCHAR, af.dataset_rank[1]);
  // set number of frames
  header->nframes = af.dataset_rank[1];

  // direction cosines (use orient_specific)
  // orient_specific : required field
  header->x_r = afni_orientations()[af.orient_specific[0]][0];
  header->x_a = afni_orientations()[af.orient_specific[0]][1];
  header->x_s = afni_orientations()[af.orient_specific[0]][2];
  header->y_r = afni_orientations()[af.orient_specific[1]][0];
  header->y_a = afni_orientations()[af.orient_specific[1]][1];
  header->y_s = afni_orientations()[af.orient_specific[1]][2];
  header->z_r = afni_orientations()[af.orient_specific[2]][0];
  header->z_a = afni_orientations()[af.orient_specific[2]][1];
  header->z_s = afni_orientations()[af.orient_specific[2]][2];

  /* --- quick determinant check --- */
  det = +header->x_r * (header->y_a * header->z_s - header->z_a * header->y_s) -
        header->x_a * (header->y_r * header->z_s - header->z_r * header->y_s) +
        header->x_s * (header->y_r * header->z_a - header->z_r * header->y_a);

  if (det == 0) {
    MRIfree(&header);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "read_afni_header(): error in orientations %d, %d, %d "
                       "(direction cosine matrix has determinant zero)",
                       af.orient_specific[0], af.orient_specific[1],
                       af.orient_specific[2]));
  }
  // sizes use delta
  // delta : required field
  header->xsize = af.delta[0];
  header->ysize = af.delta[1];
  header->zsize = af.delta[2];

  // uses origin
  // origin : required field
  header->c_r =
      header->x_r *
          (header->xsize * (header->width - 1.0) / 2.0 + af.origin[0]) +
      header->y_r *
          (header->ysize * (header->height - 1.0) / 2.0 + af.origin[1]) +
      header->z_r *
          (header->zsize * (header->depth - 1.0) / 2.0 + af.origin[2]);

  header->c_a =
      header->x_a *
          (header->xsize * (header->width - 1.0) / 2.0 + af.origin[0]) +
      header->y_a *
          (header->ysize * (header->height - 1.0) / 2.0 + af.origin[1]) +
      header->z_a *
          (header->zsize * (header->depth - 1.0) / 2.0 + af.origin[2]);

  header->c_s =
      header->x_s *
          (header->xsize * (header->width - 1.0) / 2.0 + af.origin[0]) +
      header->y_s *
          (header->ysize * (header->height - 1.0) / 2.0 + af.origin[1]) +
      header->z_s *
          (header->zsize * (header->depth - 1.0) / 2.0 + af.origin[2]);

  header->ras_good_flag = 1;

  if (header->xsize < 0) {
    header->xsize = -header->xsize;
  }

  if (header->ysize < 0) {
    header->ysize = -header->ysize;
  }

  if (header->zsize < 0) {
    header->zsize = -header->zsize;
  }

  header->imnr0 = 1;
  header->imnr1 = header->depth;

  header->ps    = header->xsize;
  header->thick = header->zsize;

  header->xend   = (header->width / 2.0) * header->xsize;
  header->xstart = -header->xend;
  header->yend   = (header->height / 2.0) * header->ysize;
  header->ystart = -header->yend;
  header->zend   = (header->depth / 2.0) * header->zsize;
  header->zstart = -header->zend;

  xfov = header->xend - header->xstart;
  yfov = header->yend - header->ystart;
  zfov = header->zend - header->zstart;

  header->fov =
      (xfov > yfov ? (xfov > zfov ? xfov : zfov) : (yfov > zfov ? yfov : zfov));

#if (BYTE_ORDER == LITTLE_ENDIAN)
  //#ifdef Linux
  swap_flag = big_endian_flag;
#else
  swap_flag = !big_endian_flag;
#endif

  if ((fp = fopen(fname, "re")) == nullptr) {
    MRIfree(&header);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE, "afniRead(): error opening file %s",
                       header_fname));
  }

  fseek(fp, 0, SEEK_END);
  brik_file_length = ftell(fp);
  fseek(fp, 0, SEEK_SET);

  // number of voxels consecutive

  nvoxels = header->width * header->height * header->depth * header->nframes;

  if ((brik_file_length % nvoxels) != 0) {
    fclose(fp);
    MRIfree(&header);
    errno = 0;
    ErrorReturn(NULL, (ERROR_BADFILE,
                       "afniRead(): BRIK file length (%d) is not divisible by "
                       "the number of voxels (%d)",
                       brik_file_length, nvoxels));
  }
  // bytes_per_voxel = brik_file_length / nvoxels; // this assumes one frame
  bytes_per_voxel =
      af.brick_types[0] + 1; // 0(byte)-> 1, 1(short) -> 2, 3(float)->4

  // bytes = header->width*header->height*header->depth*bytes_per_voxel;

  // this check is for nframes = 1 case which is the one we are supporting
  if (bytes_per_voxel != brik_file_length / nvoxels) {
    fclose(fp);
    MRIfree(&header);
    errno = 0;
    ErrorReturn(NULL, (ERROR_UNSUPPORTED,
                       "afniRead(): type info stored in header does not agree "
                       "with the file size: %d != %d/%d",
                       bytes_per_voxel, brik_file_length / nvoxels));
  }

  // if brick_float_facs != 0 then we scale values to be float
  if (af.numfacs != 0 && af.brick_float_facs[0] != 0.) {
    header->type = MRI_FLOAT;
    scaling      = af.brick_float_facs[0];
  } else {
    if (bytes_per_voxel == 1) {
      header->type = MRI_UCHAR;
    } else if (bytes_per_voxel == 2) {
      header->type = MRI_SHORT;
    } else if (bytes_per_voxel == 4) {
      header->type = MRI_FLOAT;
    } else {
      fclose(fp);
      MRIfree(&header);
      errno = 0;
      ErrorReturn(NULL,
                  (ERROR_UNSUPPORTED,
                   "afniRead(): don't know what to do with %d bytes per voxel",
                   bytes_per_voxel));
    }
  }

  ///////////////////////////////////////////////////////////////////////
  if (read_volume != 0) {
    // mri = MRIalloc(header->width, header->height, header->depth, header->type);
    mri = MRIallocSequence(header->width, header->height, header->depth,
                           header->type, header->nframes);
    MRIcopyHeader(header, mri);

    for (frame = 0; frame < header->nframes; ++frame) {
      initialized = 0;
      for (k = 0; k < mri->depth; k++) {
        for (j = 0; j < mri->height; j++) {
          if (af.brick_float_facs[frame] != 0.0F) {
            scaling = af.brick_float_facs[frame];
          } else {
            scaling = 1.;
          }
          {
            pmem = malloc(bytes_per_voxel * mri->width);
            if (pmem != nullptr) {
              if (static_cast<int>(fread(pmem, bytes_per_voxel, mri->width,
                                         fp)) != mri->width) {
                fclose(fp);
                MRIfree(&header);
                errno = 0;
                ErrorReturn(NULL,
                            (ERROR_BADFILE,
                             "afniRead(): error reading from file %s", fname));
              }
              // swap bytes
              if (swap_flag != 0) {
                if (bytes_per_voxel == 2) // short
                {
                  std::vector<uint8_t> tmp(bytes_per_voxel * mri->width);
                  swab(pmem, tmp.data(), static_cast<size_t>(mri->width * 2));
                  memcpy(pmem, tmp.data(), tmp.size());
                } else if (bytes_per_voxel == 4) // float
                {
                  pf = static_cast<float *>(pmem);
                  for (i = 0; i < mri->width; i++, pf++) {
                    *pf = swapFloat(*pf);
                  }
                }
              }
              // now scaling
              if (bytes_per_voxel == 1) // byte
              {
                pc = static_cast<unsigned char *>(pmem);
                for (i = 0; i < mri->width; i++) {
                  if (scaling == 1.) {
                    MRIseq_vox(mri, i, j, k, frame) = *pc;
                  } else {
                    MRIFseq_vox(mri, i, j, k, frame) =
                        (static_cast<float>(*pc)) * scaling;
                  }
                  ++pc;
                }
                findMinMaxByte(static_cast<unsigned char *>(pmem), mri->width,
                               &flMin, &flMax);
              }
              if (bytes_per_voxel == 2) // short
              {
                ps = static_cast<short *>(pmem);
                for (i = 0; i < mri->width; i++) {
                  // if (*ps != 0)
                  //   fmt::print("{:d} ", *ps);
                  if (scaling == 1.) {
                    MRISseq_vox(mri, i, j, k, frame) = *ps;
                  } else {
                    MRIFseq_vox(mri, i, j, k, frame) =
                        (static_cast<float>(*ps)) * scaling;
                  }
                  ++ps;
                }
                findMinMaxShort(static_cast<short *>(pmem), mri->width, &flMin,
                                &flMax);
              } else if (bytes_per_voxel == 4) // float
              {
                pf = static_cast<float *>(pmem);
                for (i = 0; i < mri->width; i++) {
                  MRIFseq_vox(mri, i, j, k, frame) = (*pf) * scaling;
                  ++pf;
                }
                findMinMaxFloat(static_cast<float *>(pmem), mri->width, &flMin,
                                &flMax);
              }
              free(pmem);
              //
              if (initialized == 0) {
                fMin        = flMin;
                fMax        = flMax;
                initialized = 1;
              } else {
                if (flMin < fMin) {
                  fMin = flMin;
                }
                if (flMax > fMax) {
                  fMax = flMax;
                }
                // fmt::print("\n fmin ={:f}, fmax = {:f}, local min = {:f}, max = {:f}\n", fMin, fMax, flMin, flMax);
              }
            } else {
              fclose(fp);
              MRIfree(&header);
              errno = 0;
              ErrorReturn(
                  NULL, (ERROR_BADFILE,
                         "afniRead(): could not allocate memory for reading %s",
                         fname));
            }
          }
        } // height
        exec_progress_callback(k, mri->depth, frame, header->nframes);
      } // depth
      // valid only for nframs == 1
      {
        fmt::print("BRICK_STATS min = {:f} <--> actual min = {:f}\n",
                   af.brick_stats[0 + 2 * frame], fMin * scaling);
        fmt::print("BRICK_STATS max = {:f} <--> actual max = {:f}\n",
                   af.brick_stats[1 + 2 * frame], fMax * scaling);
      }
    }      // nframes
  } else { // not reading volume
    mri = MRIcopy(header, nullptr);
  }

  strcpy(mri->fname, fname);

  fclose(fp);

  MRIfree(&header);
  AFclean(&af);

  return (mri);

} /* end afniRead() */

auto afniWrite(MRI *mri, const char *fname) -> int {
  char  header_fname[STRLEN];
  FILE *fp = nullptr;
  int   i  = 0;
  int   j  = 0;
  int   k  = 0;
  int   orient_specific[3];
  int   bytes_per_voxel = 0;
  float max             = NAN;
  float min             = NAN;
  int   dest_type       = 0;
  short s               = 0;
  float f               = NAN;
  char *c               = nullptr;

  errno = 0;
  // TODO(aboualiaa): implement then remove this error
  ErrorReturn(ERROR_UNSUPPORTED,
              (ERROR_UNSUPPORTED, "AFNI BRIK write unsupported"));

  /* ----- keep compiler quiet ----- */
  bytes_per_voxel = 0;
  dest_type       = -1;

  if (mri->nframes != 1) {
    errno = 0;
    ErrorReturn(ERROR_UNSUPPORTED, (ERROR_UNSUPPORTED,
                                    "afniRead(): writing of anything but 1 "
                                    "frame unsupported (mri->nframes = %d)",
                                    mri->nframes));
  }

  orient_specific[0] = -1;
  orient_specific[1] = -1;
  orient_specific[2] = -1;

  for (i = 0; i < 6; i++) {
    if (mri->x_r == afni_orientations()[i][0] &&
        mri->x_a == afni_orientations()[i][1] &&
        mri->x_s == afni_orientations()[i][2]) {
      orient_specific[0] = i;
    }
    if (mri->y_r == afni_orientations()[i][0] &&
        mri->y_a == afni_orientations()[i][1] &&
        mri->y_s == afni_orientations()[i][2]) {
      orient_specific[1] = i;
    }
    if (mri->z_r == afni_orientations()[i][0] &&
        mri->z_a == afni_orientations()[i][1] &&
        mri->z_s == afni_orientations()[i][2]) {
      orient_specific[2] = i;
    }
  }

  if (orient_specific[0] == -1 || orient_specific[1] == -1 ||
      orient_specific[2] == -1) {
    errno = 0;
    ErrorPrintf(ERROR_UNSUPPORTED,
                "afniWrite(): oblique volume writing to AFNI unsupported");
    ErrorPrintf(ERROR_UNSUPPORTED, "x_(r, a, s) = (%g, %g, %g)", mri->x_r,
                mri->x_a, mri->x_s);
    ErrorPrintf(ERROR_UNSUPPORTED, "y_(r, a, s) = (%g, %g, %g)", mri->y_r,
                mri->y_a, mri->y_s);
    ErrorPrintf(ERROR_UNSUPPORTED, "z_(r, a, s) = (%g, %g, %g)", mri->z_r,
                mri->z_a, mri->z_s);
    return (ERROR_UNSUPPORTED);
  }

  if (mri->type == MRI_INT || mri->type == MRI_LONG) {
    MRIlimits(mri, &min, &max);
    if (min > -32768.0 && min < 32768.0 && max > -32768.0 && max < 32768.0) {
      dest_type = MRI_SHORT;
    } else {
      dest_type = MRI_FLOAT;
    }
  } else if (mri->type == MRI_UCHAR) {
    bytes_per_voxel = 1;
  } else if (mri->type == MRI_SHORT) {
    bytes_per_voxel = 2;
  } else if (mri->type == MRI_FLOAT) {
    bytes_per_voxel = 4;
  } else {
    errno = 0;
    ErrorReturn(
        ERROR_UNSUPPORTED,
        (ERROR_UNSUPPORTED, "afniRead(): unsupported data type %d", mri->type));
  }

  strcpy(header_fname, fname);
  c = strrchr(header_fname, '.');

  if (c == nullptr) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "afniRead(): bad file name %s", fname));
  }

  if (strcmp(c, ".BRIK") != 0) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "afniRead(): bad file name %s", fname));
  }

  sprintf(c, ".HEAD");

  if ((fp = fopen(header_fname, "we")) == nullptr) {
    errno = 0;
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, "afniWrite(): can't open file %s for writing",
                 header_fname));
  }

  fprintf(fp, "\n");

  fprintf(fp, "type = integer-attribute\n");
  fprintf(fp, "name = ORIENT_SPECIFIC\n");
  fprintf(fp, "count = 3\n");
  fprintf(fp, " 1\n");

  fprintf(fp, "\n");

  fprintf(fp, "type = integer-attribute\n");
  fprintf(fp, "name = BRICK_TYPES\n");
  fprintf(fp, "count = 1\n");
  fprintf(fp, " 1\n");

  fprintf(fp, "\n");

  fprintf(fp, "type = integer-attribute\n");
  fprintf(fp, "name = DATASET_DIMENSIONS\n");
  fprintf(fp, "count = 3\n");
  fprintf(fp, " %d %d %d\n", mri->width, mri->height, mri->depth);

  fprintf(fp, "\n");

  fprintf(fp, "type = float-attribute\n");
  fprintf(fp, "name = DELTA\n");
  fprintf(fp, "count = 3\n");
  fprintf(fp, " %g %g %g\n", mri->xsize, mri->ysize, mri->zsize);

  fprintf(fp, "\n");

  fprintf(fp, "type = float-attribute\n");
  fprintf(fp, "name = ORIGIN\n");
  fprintf(fp, "count = 3\n");
  fprintf(fp, " %g %g %g\n", -(mri->width - 1.0) / 2.0 * mri->xsize,
          -(mri->height - 1.0) / 2.0 * mri->ysize,
          -(mri->depth - 1.0) / 2.0 * mri->zsize);

  fprintf(fp, "\n");

  fprintf(fp, "type = \n");
  fprintf(fp, "name = BYTEORDER_STRING\n");
  fprintf(fp, "count = 10\n");
#if (BYTE_ORDER == LITTLE_ENDIAN)
  //#ifdef Linux
  fprintf(fp, " LSB_FIRST~\n");
#else
  fprintf(fp, " MSB_FIRST~\n");
#endif

  fclose(fp);

  if ((fp = fopen(fname, "we")) == nullptr) {
    errno = 0;
    ErrorReturn(
        ERROR_BADFILE,
        (ERROR_BADFILE, "afniWrite(): can't open file %s for writing", fname));
  }

  for (k = 0; k < mri->depth; k++) {
    for (j = 0; j < mri->height; j++) {
      if (mri->type == MRI_INT || mri->type == MRI_LONG) {
        for (i = 0; i < mri->width; i++) {
          if (dest_type == MRI_SHORT) {
            if (mri->type == MRI_INT) {
              s = static_cast<short> MRIIvox(mri, i, j, k);
            }
            if (mri->type == MRI_LONG) {
              s = static_cast<short> MRILvox(mri, i, j, k);
            }
            fwrite(&s, sizeof(short), 1, fp);
          }
          if (dest_type == MRI_FLOAT) {
            if (mri->type == MRI_INT) {
              f = static_cast<float> MRIIvox(mri, i, j, k);
            }
            if (mri->type == MRI_LONG) {
              f = static_cast<float> MRILvox(mri, i, j, k);
            }
            fwrite(&f, sizeof(float), 1, fp);
          }
        }
      }

      else {
        fwrite(mri->slices[k][j], bytes_per_voxel, mri->width, fp);
      }
    }
  }

  fclose(fp);

  return (0);

} /* end afniWrite() */
} // namespace fs::utils

// expose to rest of freesurfer
// TODO(aboualiaa): delete once migration is complete
auto afniWrite(MRI *mri, const char *fname) -> int {
  fs::utils::afniWrite(mri, fname);
}

auto afniRead(const char *fname, int read_volume) -> MRI * {
  return fs::utils::afniRead(fname, read_volume);
}