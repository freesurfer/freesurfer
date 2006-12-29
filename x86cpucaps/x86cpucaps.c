/*
 *
 * x86cpucaps
 * by Osamu Kayasono <jacobi@jcom.home.ne.jp>
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "x86cpucaps.h"

#define B8(y)  ( ((0x##y##L & 0x00000001L) >> 0) \
               | ((0x##y##L & 0x00000010L) >> 3) \
               | ((0x##y##L & 0x00000100L) >> 6) \
               | ((0x##y##L & 0x00001000L) >> 9) \
               | ((0x##y##L & 0x00010000L) >> 12) \
               | ((0x##y##L & 0x00100000L) >> 15) \
               | ((0x##y##L & 0x01000000L) >> 18) \
               | ((0x##y##L & 0x10000000L) >> 21) )

#define B16(h,l)            (B8(h)<<8 | B8(l))
#define B32(hh, hl, lh, ll) (B16(hh, hl)<<16 | B16(lh, ll))

/* see "Identifying Supported Features" in
       AMD Processor Recognition Application Note, p10-11 */
#define FLAG_MMX      B32(00000000,10000000,00000000,00000000) /* bit 23 */
#define FLAG_SSE      B32(00000010,00000000,00000000,00000000) /* bit 25 */
#define FLAG_SSE2     B32(00000100,00000000,00000000,00000000) /* bit 26 */
#define FLAG_3DNOW    B32(10000000,00000000,00000000,00000000) /* bit 31 */
#define FLAG_3DNOWEXT B32(01000000,00000000,00000000,00000000) /* bit 30 */

#define MASK_FAMILY   B32(00000000,00000000,00001111,00000000)
#define MASK_MODEL    B32(00000000,00000000,00000000,11110000)
#define MASK_STEPPING B32(00000000,00000000,00000000,00001111)
#define MASK_ALL   ( MASK_FAMILY | MASK_MODEL | MASK_STEPPING )

#define MASK_1CHAR B32(00000000,00000000,00000000,11111111)
#define MASK_2CHAR B32(00000000,00000000,11111111,00000000)
#define MASK_3CHAR B32(00000000,11111111,00000000,00000000)
#define MASK_4CHAR B32(11111111,00000000,00000000,00000000)

#define GET_CPUNAME    0
#define GET_KERNELOPT  1
#define GET_GCCTARGET  2
#define GET_GCCSIMDOPT 3

char *x86cpucaps_detectcpu(int, float, int, int);

/* cpuid, from kernel source ( linux/include/asm-i386/processor.h ) */
inline void cpuid(int op, int *eax, int *ebx, int *ecx, int *edx) {
  __asm__("cpuid"
        : "=a" (*eax),
          "=b" (*ebx),
          "=c" (*ecx),
          "=d" (*edx)
              : "a" (op));
}


/*
 * check all SIMD capabilities
 *   struct simdcaps : checking results
 *   i : i=0 (/w check os support), i=1 (no check)
  */
int x86cpucaps_simdall(struct simdcaps *simd, int i) {
  int ex[4], f;

  memset(&(*simd),0,sizeof(struct simdcaps));

  /* check CPU has CPUID */
  cpuid(0,&ex[0],&ex[1],&ex[2],&ex[3]);
  if ( ex[0] < 1) return 1;

  cpuid(0x80000001,&ex[0],&ex[1],&ex[2],&ex[3]);
  if ( (ex[3] & FLAG_3DNOW ) == FLAG_3DNOW) simd->has_3dnow = TRUE;
  if ( (ex[3] & FLAG_3DNOWEXT ) == FLAG_3DNOWEXT) simd->has_3dnowext = TRUE;

  cpuid(1,&ex[0],&ex[1],&ex[2],&ex[3]);
//printf("Feature Flags : %08x\n", ex[3]);
  if ( (ex[3] & FLAG_MMX  ) == FLAG_MMX ) simd->has_mmx = TRUE;
  if ( (ex[3] & FLAG_SSE  ) == FLAG_SSE ) simd->has_sse = TRUE;
  if ( (ex[3] & FLAG_SSE2 ) == FLAG_SSE2) simd->has_sse2 = TRUE;

  /* SSE CPU supports mmxext too */
  if (simd->has_sse == TRUE) simd->has_mmxext = TRUE;

#ifndef _WIN32_
  if ((i == 0 ) && (simd->has_sse == TRUE)) {
    f = x86cpucaps_check_sse_supported(simd->has_sse, 0);
    if ( f == 0 ) { /* OS not support SSE, disabling to be safe */
      simd->has_sse  = FALSE;
      simd->has_sse2 = FALSE;
    }
  }
#endif /* _WIN32_ */

  return 0;
}


/*
 * check SIMD capability
 *  ( i=0:Intel SIMD(SSE2/SSE/MMX), i=1:AMD 3DNow! )
 */
int x86cpucaps_simd(int i) {
  int f = 0, k = 0;
  struct simdcaps sdall;

  memset(&sdall,0,sizeof(struct simdcaps));
  if ( i == GET_INTELSIMD_NOCHECK) k = 1;

  /* check CPU has CPUID */
  if (x86cpucaps_simdall(&sdall,k) == 1 ) return f;

  switch (i) {
  case GET_AMDSIMD:
    if (sdall.has_3dnow == TRUE ) f = HAS_3DNOW;
    if (sdall.has_3dnowext == TRUE ) f = HAS_3DNOWEXT;
    break;
  default:
    if (sdall.has_mmx  == TRUE ) f = HAS_MMX;
    if (sdall.has_sse  == TRUE ) f = HAS_SSE;
    if (sdall.has_sse2 == TRUE ) f = HAS_SSE2;
    break;
  }

  return f;
}

/*
 * check CPU Family-Model-Stepping
 *  ( i=1:Family, i=2:Model, i=3:Stepping )
 */
int x86cpucaps_cpumodel(int i) {
  int ex[4];
  int f = 0;

  /* check CPU has CPUID */
  cpuid(0,&ex[0],&ex[1],&ex[2],&ex[3]);
  if ( ex[0] < 1) return f;

  cpuid(1,&ex[0],&ex[1],&ex[2],&ex[3]);
  switch (i) {
  case 1:
    f = (ex[0] & MASK_FAMILY) >> 8;
    break;
  case 2:
    f = (ex[0] & MASK_MODEL) >> 4;
    break;
  case 3:
    f = ex[0] & MASK_STEPPING;
    break;
  default:
    f = ex[0] & MASK_ALL;
    break;
  }

  return f;
}

/*
 * check CPU Vendor
 */
int x86cpucaps_vendor(char *vendorname) {
  int ex[4];
  int f = 0;
  char vendorstr[LEN_VENDORNAME];

  /* check CPU has CPUID */
  cpuid(0,&ex[0],&ex[1],&ex[2],&ex[3]);
  if ( ex[0] < 1) return f;

  /* read Vendor Strings */
  vendorstr[0] =  ex[1] & MASK_1CHAR;
  vendorstr[1] = (ex[1] & MASK_2CHAR) >>  8;
  vendorstr[2] = (ex[1] & MASK_3CHAR) >> 16;
  vendorstr[3] = (ex[1] & MASK_4CHAR) >> 24;

  vendorstr[4] =  ex[3] & MASK_1CHAR;
  vendorstr[5] = (ex[3] & MASK_2CHAR) >>  8;
  vendorstr[6] = (ex[3] & MASK_3CHAR) >> 16;
  vendorstr[7] = (ex[3] & MASK_4CHAR) >> 24;

  vendorstr[8] =  ex[2] & MASK_1CHAR;
  vendorstr[9] = (ex[2] & MASK_2CHAR) >>  8;
  vendorstr[10]= (ex[2] & MASK_3CHAR) >> 16;
  vendorstr[11]= (ex[2] & MASK_4CHAR) >> 24;

  vendorstr[12]= '\0';

  /* check VendorName */
  if ( strcmp(vendorstr, "GenuineIntel") == 0 )
    f = VENDOR_INTEL;
  else if ( strcmp(vendorstr, "AuthenticAMD") == 0 )
    f = VENDOR_AMD;
  else if ( strcmp(vendorstr, "CyrixInstead") == 0 )
    f = VENDOR_CYRIX;
  else if ( strcmp(vendorstr, "CentaurHauls") == 0 )
    f = VENDOR_CENTAUR;
  else if ( strcmp(vendorstr, "GenuineTMx86") == 0 )
    f = VENDOR_TRANSMETA;

  /*
   * " UMC UMC UMC", UMC processor
   * "NexGenDriven", NexGen processor
   * "RiseRiseRise", Rise Technology processor
   */

  strncpy(vendorname, vendorstr, LEN_VENDORNAME);

  return f;
}


/*
 * define CPU name
 */

char *x86cpucaps_getcpuname(int vendor_id, int cpu_id) {
  return x86cpucaps_detectcpu(GET_CPUNAME, (float) 0, vendor_id, cpu_id);
}

/*
 * define kernel build option
 */
char *x86cpucaps_getkernelopt(int vendor_id, int cpu_id) {
  return x86cpucaps_detectcpu(GET_KERNELOPT, (float) 0, vendor_id, cpu_id);

}


/*
 * define gcc target
 */
char *x86cpucaps_getgcctarget(float gccver, int vendor_id, int cpu_id) {
  return x86cpucaps_detectcpu(GET_GCCTARGET, gccver, vendor_id, cpu_id);
}

/*
 * define gcc SIMD opt
 */
char *x86cpucaps_getgccsimdopt(float gccver, int vendor_id, int cpu_id) {
  return x86cpucaps_detectcpu(GET_GCCSIMDOPT, gccver, vendor_id, cpu_id);
}


char *x86cpucaps_detectcpu(int r, float gccver, int vendor_id, int cpu_id) {
  int i;
  int family, model, stepping;
  char cpuname[LEN_CPUNAME], kernelopt[LEN_KERNELOPT],
  gcctarget[LEN_GCCTARGET], gccsimdopt[LEN_GCCSIMDOPT] = "";
  char kdefopt[] = "M586";
  char gdefarch[] = "i586";
  char *retstr;

  family   = (cpu_id & MASK_FAMILY) >> 8;
  model    = (cpu_id & MASK_MODEL)  >> 4;
  stepping = cpu_id & MASK_STEPPING;

  /* it's very ugly... */
  switch (vendor_id) {

    /* Intel CPUs */
  case VENDOR_INTEL:
    if ( family == 4 ) {
      strncpy(cpuname, "i486 series", LEN_CPUNAME);
      strncpy(kernelopt, "M486", LEN_KERNELOPT);
      strncpy(gcctarget, "i486", LEN_GCCSIMDOPT);
    } else if ( family == 5 ) {
      if ( model < 4) {
        strncpy(cpuname, "Pentium Classic", LEN_CPUNAME);
        strncpy(kernelopt, "M586TSC", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium", LEN_GCCSIMDOPT);
      } else {
        strncpy(cpuname, "Pentium MMX", LEN_CPUNAME);
        strncpy(kernelopt, "M586MMX", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium-mmx", LEN_GCCSIMDOPT);
      }
    } else if ( family == 6 ) {
      if ( model <= 1 ) {
        strncpy(cpuname, "Pentium Pro", LEN_CPUNAME);
        strncpy(kernelopt, "M686", LEN_KERNELOPT);
        strncpy(gcctarget, "pentiumpro", LEN_GCCSIMDOPT);
      } else if ( model < 7 ) {
        strncpy(cpuname, "Pentium II/Pentium II Xeon/Celeron", LEN_CPUNAME);
        strncpy(kernelopt, "M686", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium2", LEN_GCCSIMDOPT);
      } else if ( model == 7) {
        strncpy(cpuname, "Pentium III Katmai", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUMIII", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium3", LEN_GCCSIMDOPT);
      } else if ( model == 8 ) {
        strncpy(cpuname, "Pentium III Coppermine", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUMIII", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium3", LEN_GCCSIMDOPT);
      } else if ( model == 9 ) {  // two version of Pentium M
        strncpy(cpuname, "Pentium M Banias", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUMIII", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium3", LEN_GCCSIMDOPT);
      } else if ( model == 10) {
        strncpy(cpuname, "Pentium III Large L2 cache", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUMIII", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium3", LEN_GCCSIMDOPT);
      } else if ( model == 11) {
        strncpy(cpuname, "Pentium III Tualatin", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUMIII", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium3", LEN_GCCSIMDOPT);
      } else if ( model == 13) {
        strncpy(cpuname, "Pentium M Dothan", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUM4", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium4", LEN_GCCSIMDOPT);
      } else {
        strncpy(cpuname, "Pentium III/Pentium III Xeon/Celeron", LEN_CPUNAME);
        strncpy(kernelopt, "MPENTIUMIII", LEN_KERNELOPT);
        strncpy(gcctarget, "pentium3", LEN_GCCSIMDOPT);
      }
    } else if ( family > 6 ) { /* family == 15 */
      // interestingly Xeon returns 15-2-7 (the same as P4-M)
      // Pentium M uses the same optimization
      strncpy(cpuname, "Pentium 4", LEN_CPUNAME);
      strncpy(kernelopt, "MPENTIUM4", LEN_KERNELOPT);
      strncpy(gcctarget, "pentium4", LEN_GCCSIMDOPT);
    } else {
      strncpy(cpuname, "i386 class", LEN_CPUNAME);
      strncpy(kernelopt, "M386", LEN_KERNELOPT);
      strncpy(gcctarget, "i386", LEN_GCCSIMDOPT);
    }
    break;

    /* AMD CPUs */
  case VENDOR_AMD:
    if ( family == 4 ) {
      if ( model <= 9 ) {
        strncpy(cpuname, "AMD i80486 series", LEN_CPUNAME);
        strncpy(kernelopt, "M486", LEN_KERNELOPT);
        strncpy(gcctarget, "i486", LEN_GCCSIMDOPT);
      } else {
        strncpy(cpuname, "AMD 5x86", LEN_CPUNAME);
        strncpy(kernelopt, "M586", LEN_KERNELOPT);
        strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
      }
    } else if ( family == 5 ) {
      if ( model <= 3 ) {
        strncpy(cpuname, "AMD K5", LEN_CPUNAME);
        strncpy(kernelopt, "M586", LEN_KERNELOPT);
        strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
      } else {
        if ( model <= 7 ) {
          strncpy(cpuname, "AMD K6", LEN_CPUNAME);
          strncpy(kernelopt, "MK6", LEN_KERNELOPT);
          strncpy(gcctarget, "k6", LEN_GCCSIMDOPT);
        } else if ( model == 8 ) {
          strncpy(cpuname, "AMD K6-2", LEN_CPUNAME);
          strncpy(kernelopt, "MK6", LEN_KERNELOPT);
          strncpy(gcctarget, "k6-2", LEN_GCCSIMDOPT);
        } else if ( model == 9 ) {
          strncpy(cpuname, "AMD K6-III", LEN_CPUNAME);
          strncpy(kernelopt, "MK6", LEN_KERNELOPT);
          strncpy(gcctarget, "k6-3", LEN_GCCSIMDOPT);
        } else {
          strncpy(cpuname, "AMD K6-2+/III+", LEN_CPUNAME);
          strncpy(kernelopt, "MK6", LEN_KERNELOPT);
          strncpy(gcctarget, "k6-3", LEN_GCCSIMDOPT);
        }
      }
    } else if ( family == 6 ) {
      if ( model == 1 ) {
        strncpy(cpuname, "AMD Athlon (K7)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon", LEN_GCCSIMDOPT);
      } else if ( model == 2 ) {
        strncpy(cpuname, "AMD Athlon (K75)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon", LEN_GCCSIMDOPT);
      } else if ( model == 3 ) {
        strncpy(cpuname, "AMD Duron (Spitfire)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon", LEN_GCCSIMDOPT);
      } else if ( model == 4 ) {
        strncpy(cpuname, "AMD Athlon (Thunderbird)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon-tbird", LEN_GCCSIMDOPT);
      } else if ( model == 6 ) {
        strncpy(cpuname, "AMD Athlon XP/MP/4 (Palomino)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon-xp", LEN_GCCSIMDOPT);
      } else if ( model == 7 ) {
        strncpy(cpuname, "AMD Duron (Morgan)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon-xp", LEN_GCCSIMDOPT);
      } else if ( model == 8 ) {
        strncpy(cpuname, "AMD Athlon XP/MP (Thoroughbred)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon-xp", LEN_GCCSIMDOPT);
      } else if ( model == 10 ) {
        strncpy(cpuname, "AMD Athlon XP/MP (Barton)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon-xp", LEN_GCCSIMDOPT);
      } else {
        strncpy(cpuname, "AMD Athlon (unknown)", LEN_CPUNAME);
        strncpy(kernelopt, "MK7", LEN_KERNELOPT);
        strncpy(gcctarget, "athlon-xp", LEN_GCCSIMDOPT);
      }
    } else if ( family == 15 ) {
      if ( model == 4) {
        strncpy(cpuname, "AMD Athlon64", LEN_CPUNAME);
      } else if ( model == 5) {
        strncpy(cpuname, "AMD Opteron/FX", LEN_CPUNAME);
      }
      strncpy(kernelopt, "MK7", LEN_KERNELOPT);
      // AMD now deprecates 3dnow and supports SSE and SSE2
      strncpy(gcctarget, "pentium4", LEN_GCCSIMDOPT);
    } else {
      strncpy(cpuname, "unknown", LEN_CPUNAME);
      strncpy(kernelopt, kdefopt, LEN_KERNELOPT);
      strncpy(gcctarget, "%s", LEN_GCCSIMDOPT);
    }
    break;

    /* VIA Cyrix CPUs */
  case VENDOR_CYRIX:
    if ( family == 4 ) {
      strncpy(cpuname, "Cyrix 5x86", LEN_CPUNAME);
      strncpy(kernelopt, "M586", LEN_KERNELOPT);
      strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
    } else if ( family == 5) {
      strncpy(cpuname, "Cyrix M1 (6x86)", LEN_CPUNAME);
      strncpy(kernelopt, "M586", LEN_KERNELOPT);
      strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
    } else if ( family == 6) {
      if ( model == 0 ) {
        strncpy(cpuname, "Cyrix M2 (6x86MX)", LEN_CPUNAME);
        strncpy(kernelopt, "M586", LEN_KERNELOPT);
        strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
      } else if ( model <= 5 ) {
        strncpy(cpuname, "VIA Cyrix III (M2 core)", LEN_CPUNAME);
        strncpy(kernelopt, "M686", LEN_KERNELOPT);
        strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
      } else if ( model == 6 ) {
        strncpy(cpuname, "VIA Cyrix III (WinChip C5A)", LEN_CPUNAME);
        strncpy(kernelopt, "MCYRIXIII", LEN_KERNELOPT);
        strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
      } else if ( model == 7 ) {
        strncpy(cpuname, "VIA Cyrix III (WinChip C5B/C)", LEN_CPUNAME);
        strncpy(kernelopt, "MCYRIXIII", LEN_KERNELOPT);
        strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
      } else {
        strncpy(cpuname, "VIA Cyrix III (WinChip C5C-T)", LEN_CPUNAME);
        strncpy(kernelopt, "MCYRIXIII", LEN_KERNELOPT);
        strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
      }
    } else {
      strncpy(cpuname, "unknown", LEN_CPUNAME);
      strncpy(kernelopt, kdefopt, LEN_KERNELOPT);
      strncpy(gcctarget, gdefarch, LEN_GCCSIMDOPT);
    }
    break;

    /* Centaur CPUs */
  case VENDOR_CENTAUR:
    if ( family == 5) {
      if ( model <= 4) {
        strncpy(cpuname, "Centaur WinChip C6", LEN_CPUNAME);
        strncpy(kernelopt, "MWINCHIPC6", LEN_KERNELOPT);
        strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
      } else if ( model <= 8 ) {
        strncpy(cpuname, "Centaur WinChip 2", LEN_CPUNAME);
        strncpy(kernelopt, "MWINCHIP2", LEN_KERNELOPT);
        strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
      } else {
        strncpy(cpuname, "Centaur WinChip 2A", LEN_CPUNAME);
        strncpy(kernelopt, "MWINCHIP3D", LEN_KERNELOPT);
        strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
      }
      break;
    } else {
      strncpy(cpuname, "unknown", LEN_CPUNAME);
      strncpy(kernelopt, kdefopt, LEN_KERNELOPT);
      strncpy(gcctarget, gdefarch, LEN_GCCSIMDOPT);
    }
    break;

    /* Transmeta CPUs */
  case VENDOR_TRANSMETA:
    strncpy(cpuname, "Transmeta Crusoe TM3x00/5x00", LEN_CPUNAME);
    strncpy(kernelopt, "MCRUSOE", LEN_KERNELOPT);
    strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
    break;

    /* Others */
  default:
    strncpy(cpuname, "unknown", LEN_CPUNAME);
    strncpy(kernelopt, kdefopt, LEN_KERNELOPT);
    strncpy(gcctarget, gdefarch, LEN_GCCSIMDOPT);
    break;
  }

  /* some targets not supported by older gcc */
  if ( gccver < (float)3.4 ) {
    if (strstr(gcctarget, "k8") != NULL ) strncpy(gcctarget, "athlon-xp", LEN_GCCSIMDOPT);
  }
  if ( gccver < (float)3.1 ) {
    if (strstr(gcctarget, "athlon-") != NULL )
      strncpy(gcctarget, "athlon", LEN_GCCSIMDOPT);
    else if (strstr(gcctarget, "k6-") != NULL )
      strncpy(gcctarget, "k6", LEN_GCCSIMDOPT);
    else if (strstr(gcctarget, "pentium-mmx") != NULL )
      strncpy(gcctarget, "pentium", LEN_GCCSIMDOPT);
    else if ((strstr(gcctarget, "pentium2") != NULL ) ||
             (strstr(gcctarget, "pentium3") != NULL ) ||
             (strstr(gcctarget, "pentium4") != NULL ) )
      strncpy(gcctarget, "pentiumpro", LEN_GCCSIMDOPT);

    if ( gccver < (float)3.0) {
      if (strstr(gcctarget, "athlon") != NULL )
        strncpy(gcctarget, "pentiumpro", LEN_GCCSIMDOPT);
      else if (strstr(gcctarget, "k6") != NULL )
        strncpy(gcctarget, "pentium", LEN_GCCSIMDOPT);
    }

    if ( gccver < (float)2.9) {
      if (strstr(gcctarget, "pentiumpro") != NULL )
        strncpy(gcctarget, "i686", LEN_GCCSIMDOPT);
      else if (strstr(gcctarget, "pentium") != NULL )
        strncpy(gcctarget, "i586", LEN_GCCSIMDOPT);
    }
  }

  /* added SIMD options */
  if (gccver >= 3.1) {
    i = x86cpucaps_simd(GET_AMDSIMD);
    if (( i == HAS_3DNOW ) || ( i == HAS_3DNOWEXT))
      strncpy(gccsimdopt, "-m3dnow", LEN_GCCSIMDOPT);
    else {
      i = x86cpucaps_simd(GET_INTELSIMD);
      switch (i) {
      case HAS_SSE2:
        strncpy(gccsimdopt, "-msse2 -mfpmath=sse", LEN_GCCSIMDOPT);
        break;
      case HAS_SSE:
        strncpy(gccsimdopt, "-msse", LEN_GCCSIMDOPT);
        break;
      case HAS_MMX:
        strncpy(gccsimdopt, "-mmmx", LEN_GCCSIMDOPT);
        break;
      }
    }
  }

  /* out result */
  switch (r) {
  case GET_CPUNAME:
    if ((retstr=(char *)malloc((LEN_CPUNAME+1)*sizeof(char))) == NULL)
      return NULL;
    strncpy(retstr, cpuname, LEN_CPUNAME);
    break;
  case GET_GCCTARGET:
    if ((retstr=(char *)malloc((LEN_GCCTARGET+1)*sizeof(char))) == NULL)
      return NULL;
    snprintf(retstr, LEN_GCCTARGET, "-march=%s", gcctarget);
    break;
  case GET_GCCSIMDOPT:
    if ((retstr=(char *)malloc((LEN_GCCSIMDOPT+1)*sizeof(char))) == NULL)
      return NULL;
    strncpy(retstr, gccsimdopt, LEN_GCCSIMDOPT);
    break;
  default:
    if ((retstr=(char *)malloc((LEN_KERNELOPT+1)*sizeof(char))) == NULL)
      return NULL;
    snprintf(retstr, LEN_KERNELOPT, "CONFIG_%s", kernelopt);
    break;
  }

  return retstr;
}
