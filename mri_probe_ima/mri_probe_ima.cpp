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

#ifndef Darwin
#include <malloc.h>
#endif
#include "diag.h"
#include "imautils.h"
#include "mri2.h"
#include "version.h"

int main(int argc, char *argv[]);

const char *Progname = NULL;

static int  parse_commandline(int argc, char **argv);
static void check_options();
static void print_usage();
static void usage_exit();
static void print_help();
static void print_version();
static void argnerr(char *option, int n);
static int  singledash(char *flag);
static int  stringmatch(const char *s1, const char *s2);

char *      imafile      = NULL;
char *      typestring   = NULL;
int         type         = -1;
int         typesize     = 1;
int         offset       = -1;
int         stringlen    = 1;
const char *key          = NULL;
int         keyno        = -1;
int         dumpfileinfo = 0;
int         debug, verbose;
FILE *      fp;
char *      attrname;
int         getattr = 0;

const char *bstem = "img";
short *     pixeldata;
int         npixels;

#define TMPSTRLEN 10000
static char tmpstr[TMPSTRLEN];

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  void *       pVal;
  IMAFILEINFO *ifi;
  int          nargs;

  nargs = handleVersionOption(argc, argv, "mri_probe_ima");
  if (nargs && argc - nargs == 1)
    exit(0);
  argc -= nargs;

  tmpstr[0] = 'a'; /* to stop compiler warning */

  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL);
  DiagInit(nullptr, nullptr, nullptr);

  if (argc == 0) {
    usage_exit();
  }

  parse_commandline(argc, argv);
  check_options();

  if (getattr != 0) {
    if (stringmatch(attrname, "isima") != 0) {
      printf("%d\n", imaIsSiemensIMA(imafile));
      return (0);
    }
  }

  if (imaIsSiemensIMA(imafile) == 0) {
    printf("WARNING: %s does not seem to be a Siemens IMA file\n", imafile);
  }

  MkImaDictionary();

  if (dumpfileinfo != 0) {
    ifi = imaLoadFileInfo(imafile);
    imaDumpFileInfo(stdout, ifi);
    return (0);
  }

  if (getattr != 0) {
    ifi = imaLoadFileInfo(imafile);
    if (stringmatch(attrname, "studydate") != 0) {
      printf("%s\n", ifi->StudyDate);
      return (0);
    }
    if (stringmatch(attrname, "studytime") != 0) {
      printf("%s\n", ifi->StudyTime);
      return (0);
    }
    if (stringmatch(attrname, "voldim") != 0) {
      printf("%d %d %d\n", ifi->VolDim[0], ifi->VolDim[1], ifi->VolDim[2]);
      return (0);
    }
    if (stringmatch(attrname, "volres") != 0) {
      printf("%g %g %g\n", ifi->VolRes[0], ifi->VolRes[1], ifi->VolRes[2]);
      return (0);
    }
    if (stringmatch(attrname, "nframes") != 0) {
      printf("%d\n", ifi->NFrames);
      return (0);
    }
    if (stringmatch(attrname, "tr") != 0) {
      printf("%g\n", ifi->RepetitionTime);
      return (0);
    }
    if (stringmatch(attrname, "pulseseq") != 0) {
      printf("%s\n", ifi->PulseSequence);
      return (0);
    }
    if (stringmatch(attrname, "patname") != 0) {
      printf("%s\n", ifi->PatientName);
      return (0);
    }
    if (stringmatch(attrname, "patdob") != 0) {
      printf("%s\n", ifi->PatientDOB);
      return (0);
    }
    if (stringmatch(attrname, "patgender") != 0) {
      printf("%s\n", ifi->PatientGender);
      return (0);
    }
    if (stringmatch(attrname, "ismosaic") != 0) {
      printf("%d\n", ifi->IsMosaic);
      return (0);
    }
    if (stringmatch(attrname, "nfilesact") != 0) {
      printf("%d\n", ifi->NFilesInSeries);
      return (0);
    }
    if (stringmatch(attrname, "nfilesexp") != 0) {
      printf("%d\n", ifi->NFilesInSeriesExp);
      return (0);
    }
    if (stringmatch(attrname, "nfilesperframe") != 0) {
      printf("%d\n", ifi->NFilesPerFrame);
      return (0);
    }
    if (stringmatch(attrname, "error") != 0) {
      printf("%d\n", ifi->ErrorFlag);
      return (0);
    }
    if (stringmatch(attrname, "pixeldata") != 0) {
      npixels   = ifi->NImageRows * ifi->NImageCols;
      pixeldata = imaReadPixelData(ifi, nullptr);
      if (pixeldata == nullptr) {
        printf("ERROR: could not read pixel data\n");
        exit(1);
      }

      sprintf(tmpstr, "%s_000.bshort", bstem);
      fp = fopen(tmpstr, "we");
      if (fp == nullptr) {
        printf("ERROR: cannot open %s for writing\n", tmpstr);
        exit(1);
      }
      fwrite(pixeldata, sizeof(short), npixels, fp);
      fclose(fp);

      sprintf(tmpstr, "%s_000.hdr", bstem);
      fp = fopen(tmpstr, "we");
      if (fp == nullptr) {
        printf("ERROR: cannot open %s for writing\n", tmpstr);
        exit(1);
      }
      fprintf(fp, "%d %d 1 %d\n", ifi->NImageRows, ifi->NImageCols, Arch486());
      fclose(fp);

      return (0);
    }

    printf("ERROR: attribute %s not recognized\n", attrname);
    return (1);
  }

  fp = fopen(imafile, "re");
  if (fp == nullptr) {
    printf("ERROR: could not open %s\n", imafile);
    exit(1);
  }

  if (keyno > -1) {
    key = ImaDictionary[keyno].key;
  }
  if (key != nullptr) {
    pVal = imaLoadValFromKey(fp, key, nullptr);
    type = imaTypeFromKey(key);
    imaPrintVal(stdout, type, pVal);
    return (0);
  }

  if (offset > -1) {
    type     = imaTypeFromString(typestring);
    typesize = imaTypeSize[type];
    if (debug != 0) {
      printf("type = %s (%d), offset = %d\n", typestring, type, offset);
    }
    pVal = imaLoadVal(fp, offset, typesize, stringlen, nullptr);
    imaPrintVal(stdout, type, pVal);
    printf("\n");
    return (0);
  }

  DumpImaDictionaryVal(stdout, imafile);

  fclose(fp);
  return (0);
}
/*---------------------------------------------------------------*/
/* -----//////////////////<<<<<<<<>>>>>>>>>>\\\\\\\\\\\\\\\\---- */
/*---------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int    nargc;
  int    nargsused;
  char **pargv;
  char * option;

  if (argc < 1) {
    usage_exit();
  }

  nargc = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug != 0) {
      printf("%d %s\n", nargc, option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (strcasecmp(option, "--help") == 0) {
      print_help();
    } else if (strcasecmp(option, "--version") == 0) {
      print_version();
    } else if (strcasecmp(option, "--debug") == 0) {
      debug = 1;
    } else if (strcasecmp(option, "--verbose") == 0) {
      verbose = 1;
    } else if (strcasecmp(option, "--fileinfo") == 0) {
      dumpfileinfo = 1;

      /* -------- source volume inputs ------ */
    } else if (strcmp(option, "--i") == 0) {
      if (nargc < 1) {
        argnerr(option, 1);
      }
      imafile   = pargv[0];
      nargsused = 1;
    } else if (strcmp(option, "--attr") == 0) {
      if (nargc < 1) {
        argnerr(option, 1);
      }
      attrname  = pargv[0];
      getattr   = 1;
      nargsused = 1;
    } else if (strcmp(option, "--ob") == 0) {
      if (nargc < 1) {
        argnerr(option, 1);
      }
      bstem     = pargv[0];
      nargsused = 1;
    } else if (strcmp(option, "--dictionary") == 0) {
      DumpImaDictionary(stdout);
      exit(0);
    } else if (strcmp(option, "--key") == 0) {
      if (nargc < 1) {
        argnerr(option, 1);
      }
      key       = pargv[0];
      nargsused = 1;
    } else if (strcmp(option, "--keyno") == 0) {
      if (nargc < 1) {
        argnerr(option, 1);
      }
      sscanf(pargv[0], "%d", &keyno);
      nargsused = 1;
    } else if (strcmp(option, "--o") == 0) {
      if (nargc < 2) {
        argnerr(option, 2);
      }
      sscanf(pargv[0], "%d", &offset);
      typestring = pargv[1];
      nargsused  = 2;
      if (strcmp(typestring, "string") == 0) {
        if (nargc < 3) {
          printf("ERROR: type string needs length argument\n");
          exit(1);
        }
        sscanf(pargv[2], "%d", &stringlen);
        if (stringlen < 1) {
          printf("ERROR: string length = %d, must be >= 1\n", stringlen);
          exit(1);
        }
        nargsused = 3;
      }
    } else {
      fprintf(stderr, "ERROR: Option %s unknown\n", option);
      if (singledash(option) != 0) {
        fprintf(stderr, "       Did you really mean -%s ?\n", option);
      }
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return (0);
}
/* ------------------------------------------------------ */
static void usage_exit() {
  print_usage();
  exit(1);
}
/* --------------------------------------------- */
static void print_usage() {
  fprintf(stdout, "USAGE: %s \n", Progname);
  fprintf(stdout, "\n");
  fprintf(stdout, "   --i imafile           : path to ima file \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --key keystring       : string from dictionary\n");
  fprintf(stdout,
          "   --o offset type <len> : decimal offset, data type, and len\n");
  fprintf(stdout,
          "   --attr attrname       : print value of file info attribute\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --fileinfo            : print interpreted file info\n");
  fprintf(stdout,
          "   --dictionary          : print dictionary (no file needed)\n");
#if 0
  fprintf(stdout, "   --view          : view the image  \n");
  fprintf(stdout, "   --o file        : dump binary pixel data into file\n");
#endif
  fprintf(stdout,
          "   --ob stem             : dump binary pixel data into bshort\n");
  fprintf(stdout, "   --help                : how to use this program \n");
  fprintf(stdout, "   --version             : print version and exit\n");

  fprintf(stdout, "\n");
}
/* --------------------------------------------- */
static void print_help() {
  printf("\n");
  print_usage();
  printf("\n");

  printf("This program allows the user to query a Siemens IMA file, \n"
         "and can be used to print out a single value from the IMA header \n"
         "or to dump lots of info.\n");

  printf("This program allows the user to query a Siemens IMA file, \n"
         "and can be used to print out a single value from the IMA header \n"
         "or to dump lots of info.\n");

  printf(
      "\n"
      "DESCRIPTION\n"
      "\n"
      "A single value can be queried in one of three ways: (1) offset/type, \n"
      "(2) key string, (3) attribute name. In the offset method, the user "
      "supplies \n"
      "an offset which is the number of bytes into the header; the user also "
      "supplies  \n"
      "a type string indicating the data type (short, int, long, float, "
      "double, string). \n"
      "If a string type is specified, the user must also supply the string "
      "length. In  \n"
      "the key method, the user supplies a key string from a dictionary (see "
      "below). \n"
      "In the attribute method, the user supplies an attribute name. The "
      "difference \n"
      "between an attribute and a key is that an attribute may have been "
      "interpreted \n"
      "in some way (eg, the raw (key-based) row resolution for mosaics is "
      "wrong, \n"
      "but the attribute is correct). \n"
      "\n"
      "A dump of header information can be obtained in two ways. First, if "
      "only \n"
      "a file name is given, then the dictionary with corresponding values "
      "will  \n"
      "be dumped. If the --fileinfo flag is added, then values that have been  "
      "\n"
      "interpreted are printed out. \n"
      "\n"
      "ARGUMENTS\n"
      "\n"
      "  --i imafile\n"
      "\n"
      "      Path to the IMA file to be probe. If this is the only option, \n"
      "      the dictionary with corresponding values is printed out. See "
      "also\n"
      "      --dictionary.\n"
      "\n"
      "  --o offset type <stringlen>\n"
      "\n"
      "      offset is the number of bytes from the beginning of the file. "
      "type \n"
      "      is the type of data. Valid values are short, int, long, float, "
      "double, \n"
      "      and string. If string is used, the length of the string must be "
      "supplied.\n"
      "\n"
      "  --d key \n"
      "\n"
      "      key is a string as found in the dictionary (see DICTIONARY "
      "below).\n"
      "\n"
      "  --attr attrname\n"
      "\n"
      "     Name of an attribute.\n"
      "\n"
      "  --fileinfo \n"
      "\n"
      "     Dump the interpreted file information.\n"
      "\n"
      "  --dictionary\n"
      "\n"
      "     Dump the dictionary (no ima file name need be supplied). Each "
      "entry has\n"
      "     six columns: entry number, key string, offset, type string, number "
      "of \n"
      "     bytes in the type, and the string length. If only the imafile is "
      "supplied,\n"
      "     the dictionary will be printed out with the value as a seventh "
      "column.\n"
      "\n"
      "  --ob bstem\n"
      "\n"
      "     Save pixel data to bstem_000.bshort instead of img_000.bshort\n"
      "\n"
      "DICTIONARY\n"
      "\n"
      "The dictionary is a list of character strings (keys) that describe an \n"
      "aspect of a value in the IMA header along with its offset, data type,  "
      "\n"
      "and, if a string, the length of the string. The dictionary is unique to "
      "\n"
      "to this program and does not represent anything official from Siemens \n"
      "or anyone else. The key names, offsets, and data types were gleaned "
      "from \n"
      "other programs and reverse engineering. \n"
      "\n"
      "ATTRIBUTES\n"
      "\n"
      "  isima     : returns 1 if the file is a Siemens IMA file, 0 otherwise "
      "\n"
      "  studydate : date of the scan  (YYYYMMDD)\n"
      "  studytime : time of the scan  (HHMMSS)\n"
      "  voldim    : number of columns, rows, and slices in the volume. \n"
      "  volres    : spacing between columns, rows, and slices in the volume\n"
      "  nframes   : number of frames (expected)\n"
      "  ismosaic  : returns 1 if slices are mosaiced, 0 otherwise\n"
      "  nfilesact : actual number of files found for this series\n"
      "  nfilesexp : number of files expected for this series\n"
      "  error     : 1 if actual number of files does not equal expected, 0 "
      "otherwise\n"
      "  tr        : repetition time (sec)\n"
      "  pulseseq  : pulse sequence name\n"
      "  patname   : patient name \n"
      "  patdob    : patient date of birth (YYYYMMDD)\n"
      "  patgender : patient gender\n"
      "  pixeldata : stores pixel data as a 2D image in img_000.bshort\n"
      "\n"
      "AUTHOR\n"
      "\n"
      "Written by Douglas N. Greve.\n"
      "\n"
      "BUG REPORTING\n"
      "\n"
      "Send bug reports to analysis-bugs@nmr.mgh.harvard.edu. Make sure to "
      "include\n"
      "enough information to replicate the problem. This includes, but is not "
      "limited\n"
      "to, the version and command-line.\n"
      "\n"
      "BUGS\n"
      "\n"
      "This program does not make use of any 'official' documentation from "
      "Siemens. \n"
      "As such, some of the results may be wrong. Some of the header elements "
      "may or\n"
      "may not be filled with legitimate values, depending upon how the pulse "
      "sequence\n"
      "was programmed.\n"
      "\n"
      "The G28_Pre_PixelSize_Row and G28_Pre_PixelSize_Column values are \n"
      "incorrect for mosaics. This is actually the way Siemens stores \n"
      "these values and is not an bug in  this program.\n"
      "\n"
      "G20_Rel_Study (3200, long) appears to be the Series Number. There is\n"
      "another variable (G20_Rel_Series, 3204, long), but this has garbage in\n"
      "it, at least in the MGH-NMR files.\n"
      "\n"
      "\n"
      "VERSION\n"
      "\n");
  printf("   %s\n\n", getVersion().c_str());

  exit(1);
}
/* --------------------------------------------- */
static void print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str());
  exit(1);
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n == 1) {
    fprintf(stderr, "ERROR: %s flag needs %d argument\n", option, n);
  } else {
    fprintf(stderr, "ERROR: %s flag needs %d arguments\n", option, n);
  }
  exit(-1);
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) {
    return (0);
  }

  if (flag[0] == '-' && flag[1] != '-') {
    return (1);
  }
  return (0);
}
/* --------------------------------------------- */
static void check_options() {

  if (imafile == nullptr) {
    printf("ERROR: no file name supplied\n");
    exit(1);
  }

  return;

  if (key != nullptr && offset > -1) {
    printf("ERROR: cannot specify key and offset\n");
    exit(1);
  }

  if (key == nullptr && offset < 0) {
    printf("ERROR: must specify either key or offset\n");
    exit(1);
  }

  if (offset > -1 && typestring == nullptr) {
    fprintf(stderr, "ERROR: no data type supplied\n");
    exit(1);
  }
}
/*-------------------------------------------------------------*/
static int stringmatch(const char *s1, const char *s2) {
  if (strcmp(s1, s2) == 0)
    return (1);
  return (0);
}
