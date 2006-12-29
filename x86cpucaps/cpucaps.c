/*
 *
 * 'cpucaps', sample code for libx86cpucaps
 * by Osamu Kayasono <jacobi@jcom.home.ne.jp>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "x86cpucaps.h"
#ifndef _WIN32_
#include <getopt.h>
#endif

#define OUT_ALL       0
#define OUT_KERNELOPT 1
#define OUT_GCCTARGET 2
#define OUT_GCCSIMD   3
#define OUT_GCCALL    4
#define OUT_MODEL     5

struct x86cpucaps cpucaps;
struct simdcaps sdcaps;

#define DEBUG 1 // set to 1, verbose

void print_usage() {
  printf("Usage: %s [OPTIONS]\n\n", PACKAGE);
  printf("    -h, --help                   this help\n");
  printf("    -w, --withgccversion=GCCVER  specify gcc version. \n");
  printf("                                 'cpucaps -w 2.953', if you use gcc-2.95.3\n");
  printf("    -g, --outgccopt              print optimal gcc target only\n");
  printf("    -s, --outsimdopt             print optimal gcc SIMD options only\n");
  printf("    -k, --outkernelopt           print optimal ProcessorType for\n");
  printf("                                 kernel building only\n");
  printf("    -m, --outmodel               print cpu model\n");
  printf("    -v, --version                print version and exit\n\n");

  return;
}

void print_version() {
  printf("%s version %s\n", PACKAGE, VERSION);
  printf("Copyright (C) 2002 Osamu Kayasono\n");
  printf("This is free software. There is NO WARRANTY,\n");
  printf("to the extent permitted by law.\n");
  printf("You may redistribute copies of %s\n", PACKAGE);
  printf("under the terms of the GNU General Public License.\n");
  printf("For more information about these matters,\n");
  printf("see the files named COPYING.\n");

  return;
}

void print_simdcaps(int i) {
  switch (i) {
  case HAS_SSE2:
    printf("SSE2\n");
    break;
  case HAS_SSE:
    printf("SSE\n");
    break;
  case HAS_3DNOWEXT:
    printf("3DNow! extensions\n");
    break;
  case HAS_3DNOW:
    printf("3DNow!\n");
    break;
    /*    case HAS_MMXEXT:
          printf("MMX extensions\n");
          break; */
  case HAS_MMX:
    printf("MMX\n");
    break;
  default:
    printf("none\n");
  }
  return;
}

void print_cpuinfo(int outmode) {
  switch (outmode) {
  case OUT_KERNELOPT:
    printf("%s\n", cpucaps.kernelopt);
    break;
  case OUT_GCCTARGET:
  case OUT_GCCSIMD:
  case OUT_GCCALL:
    if (outmode != OUT_GCCSIMD ) printf("%s ", cpucaps.gcctarget);
    if (outmode != OUT_GCCTARGET ) printf("%s", cpucaps.gccsimdopt);
    printf("\n");
    break;
  case OUT_MODEL:
    printf("%s\n", cpucaps.cpu_name);
    break;
  default:
    printf("\n");
    printf("CPU info                           : %d - %d - %d (%s)\n",
           cpucaps.cpu_family, cpucaps.cpu_model, cpucaps.cpu_stepping,
           cpucaps.vendor_name);
    printf("CPU Model Name                     : %s\n", cpucaps.cpu_name);
    printf("Recommended Kernel building option : %s\n", cpucaps.kernelopt);
    printf("Recommended gcc (%.4f) target    : %s %s\n", cpucaps.gccver, cpucaps.gcctarget, cpucaps.gccsimdopt);

    printf("checking Intel SIMD capability     : ");
    print_simdcaps(cpucaps.intel_simd);

    printf("checking AMD 3DNow! capability     : ");
    print_simdcaps(cpucaps.amd_simd);
    printf("\n");

#ifndef _WIN32_
    if (sdcaps.has_sse == TRUE)
      sdcaps.has_sse= x86cpucaps_check_sse_supported(sdcaps.has_sse,DEBUG);
#endif

    printf("SIMD capabilities checking results\n");
    printf("   SSE2:%d, SSE:%d, MMXext:%d, MMX:%d,  3DNow!Ex:%d, 3DNow!:%d\n",
           sdcaps.has_sse2,sdcaps.has_sse,sdcaps.has_mmxext, sdcaps.has_mmx,
           sdcaps.has_3dnowext,sdcaps.has_3dnow);

  }

  return;
}

float getgccver(char *gccver) {
  const float defver = 3.11;
  float ver;
  char *tmpver;
  int i;

  tmpver = (char *)malloc((strlen(gccver)+1)*sizeof(char));
  /* can't allocale, use default */
  if (tmpver == NULL) return (float) defver;

  memset(tmpver, 0, strlen(gccver)+1);

  for (i=0;i<strlen(gccver);i++) {
    if ( (gccver[i] >= '1') && (gccver[i] <= '9')) break;
  }
  /* 'gccver' is not numeric, use default */
  if ( i == strlen(gccver)-1 ) return (float) defver;

  strncpy(tmpver, gccver, strlen(gccver));
  memmove(&tmpver[0], &tmpver[i], strlen(gccver)-i+1);

  ver = (float) atof(tmpver);
  free(tmpver);

  return ver;
}


int main(int argc, char *argv[]) {
  int c= 0;
  int outmode = OUT_ALL;
  int cpu_id;
  char *tmpstr;

  cpucaps.gccver = 3.11; /* default GCC version, 3.1.1 */
  cpucaps.vendor_id = x86cpucaps_vendor(cpucaps.vendor_name);

#ifndef _WIN32_
  while (1) {
    /* int this_option_optind = optind ? optind : 1; */
    int option_index = 0;
    static struct option long_options[] = {
                                            {"help", 0, NULL, 'h'
                                            },
                                            {"version", 0, NULL, 'v'},
                                            {"withgccversion", 1, NULL, 'w'},
                                            {"outgccopt", 0, NULL, 'g'},
                                            {"outsimdopt", 0, NULL, 's'},
                                            {"outkernelopt", 0, NULL, 'k'}
                                          };

    c = getopt_long (argc, argv, "hvw:gskm", long_options, &option_index);
    if (c == EOF ) break;

    switch (c) {
    case 'w':
      if (optarg) cpucaps.gccver = getgccver(optarg);
      break;

    case 'g':
      if (outmode == OUT_GCCSIMD)
        outmode = OUT_GCCALL;
      else
        outmode = OUT_GCCTARGET;
      break;

    case 'k':
      outmode = OUT_KERNELOPT;
      break;

    case 'm':
      outmode = OUT_MODEL;
      break;

    case 's':
      if (outmode == OUT_GCCTARGET)
        outmode = OUT_GCCALL;
      else
        outmode = OUT_GCCSIMD;
      break;

    case 'v':
      print_version();
      return 0;
      break;

    case 'h':
    case '?':
    default:
      print_usage();
      return 0;
      break;
    }
  }
#endif /* _WIN32_ */

  cpucaps.cpu_family   = x86cpucaps_cpumodel(GET_FAMILY);
  cpucaps.cpu_model    = x86cpucaps_cpumodel(GET_MODEL);
  cpucaps.cpu_stepping = x86cpucaps_cpumodel(GET_STEPPING);
  /*   cpu_id = x86cpucaps_cpumodel(GET_ALLID); */

  cpu_id = (cpucaps.cpu_family << 8) + (cpucaps.cpu_model << 4) + cpucaps.cpu_stepping;

  tmpstr =  x86cpucaps_getcpuname(cpucaps.vendor_id, cpu_id);
  strncpy(cpucaps.cpu_name, tmpstr, LEN_CPUNAME);
  free(tmpstr);

  tmpstr = x86cpucaps_getkernelopt(cpucaps.vendor_id, cpu_id);
  strncpy(cpucaps.kernelopt, tmpstr, LEN_KERNELOPT);
  free(tmpstr);

  tmpstr = x86cpucaps_getgcctarget(cpucaps.gccver, cpucaps.vendor_id, cpu_id);
  strncpy(cpucaps.gcctarget, tmpstr, LEN_GCCTARGET);
  free(tmpstr);

  tmpstr = x86cpucaps_getgccsimdopt(cpucaps.gccver, cpucaps.vendor_id, cpu_id);
  strncpy(cpucaps.gccsimdopt, tmpstr, LEN_GCCSIMDOPT);
  free(tmpstr);

  cpucaps.intel_simd = x86cpucaps_simd(GET_INTELSIMD);
  cpucaps.amd_simd = x86cpucaps_simd(GET_AMDSIMD);

  x86cpucaps_simdall(&sdcaps,0);

  print_cpuinfo(outmode);

  return 0;
}
