#define FALSE 0
#define TRUE  1

#define LEN_VENDORNAME 13
#define LEN_CPUNAME    40
#define LEN_KERNELOPT  32
#define LEN_GCCSIMDOPT     32
#define LEN_GCCTARGET    40 /* -march=GCCSIMDOPT */

struct simdcaps {
  unsigned int has_sse2;
  unsigned int has_sse;
  unsigned int has_3dnowext;
  unsigned int has_3dnow;
  unsigned int has_mmxext;
  unsigned int has_mmx;
};

struct x86cpucaps {
  int vendor_id;
  char vendor_name[LEN_VENDORNAME];
  int cpu_family;
  int cpu_model;
  int cpu_stepping;
  char cpu_name[LEN_CPUNAME];
  char kernelopt[LEN_KERNELOPT];
  float gccver;
  char gcctarget[LEN_GCCTARGET];
  char gccsimdopt[LEN_GCCSIMDOPT];
  int intel_simd;
  int amd_simd;
};

/*
 * x86cpucaps_simd : check SIMD capability
 * i=0: check Intel SIMD (SSE2/SSE/MMX)
 * i=1: check AMD 3DNow!
 */
extern int x86cpucaps_simdall(struct simdcaps *simd, int i);
extern int x86cpucaps_check_sse_supported(int i, int verbose);
extern int x86cpucaps_simd(int i);

#define GET_INTELSIMD           0
#define GET_AMDSIMD             1
#define GET_INTELSIMD_NOCHECK   2

/* result value of 'x86cpucaps_simd' */
#define HAS_NOSIMD      0x00
#define HAS_MMX  0x01
/* #define HAS_MMXEXT  0x02 */
#define HAS_SSE  0x03
#define HAS_SSE2 0x04
#define HAS_3DNOW 0x11
#define HAS_3DNOWEXT 0x12

/*
 * check CPU Family-Model-Stepping
 *  ( i=0: Family * 256 + Model * 16 + Stepping,
 *    i=1:Family, i=2:Model, i=3:Stepping )
 */
extern int x86cpucaps_cpumodel(int i);

#define GET_ALLID    0
#define GET_FAMILY   1
#define GET_MODEL    2
#define GET_STEPPING 3

/*
 * get Vendor ID
 */
extern int x86cpucaps_vendor(char *vendorname);

/* result value of 'x86cpucaps_vendor' */
#define VENDOR_INTEL     1
#define VENDOR_AMD       2
#define VENDOR_CYRIX     3
#define VENDOR_CENTAUR   4
#define VENDOR_TRANSMETA 5
#define VENDOR_OTHERS    0

extern char *x86cpucaps_getcpuname(int vendor_id, int cpu_id);

extern char *x86cpucaps_getkernelopt(int vendor_id, int cpu_id);

extern char *x86cpucaps_getgcctarget(float gccver, int vendor_id, int cpu_id);

extern char *x86cpucaps_getgccsimdopt(float gccver, int vendor_id, int cpu_id);
