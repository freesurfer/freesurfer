#include <string>
#include <iostream>
#include <fstream>

#ifdef __cplusplus
extern "C"
{
#endif

#include "error.h"
#include "utils.h"
#include "macros.h"
#include "mri.h"
#include "version.h"
#include "transform.h"
#include "gcamorph.h"

#ifdef __cplusplus
}
#endif

using namespace std;

static GCAM *concatenate(GCA_MORPH *gcam1, GCA_MORPH *gcam2);
static int GCAMsampleOrig(const GCA_MORPH *gcam, float x, float y, float z, float *pxd, float *pyd, float *pzd);
static int boundsCheckf(float x, float y, float z, int width, int height, int depth);
static bool parseCommandLine(int argc, char *argv[], string &warp1name, string &warp2name, string &outname);

static int parseOption(int argc, char *argv[]);
static void usage(void);


int main(int argc, char *argv[])
{
  string warp1name, warp2name, outname;

  if (!parseCommandLine(argc, argv, warp1name, warp2name, outname))
  { exit(1); }

  // read warp 1
  cout << "reading " << warp1name << endl;
  GCA_MORPH *gcam1 = GCAMread(warp1name.c_str());
  if (!gcam1) exit(1);

  // read warp 2
  cout << "reading " << warp2name << endl;
  GCA_MORPH *gcam2 = GCAMread(warp2name.c_str());
  if (!gcam2) exit(1);

  // begin concatenation
  GCAM *gcam_concat = concatenate(gcam1, gcam2);
  if (!gcam_concat) {
    ErrorExit(ERROR_BADPARM, "could not concatenate!");
  }

  // write warp:
  printf("warp concatenation successful!\n");
  GCAMwrite(gcam_concat, outname.c_str());

  // clean up & exit:
  GCAMfree(&gcam1);
  GCAMfree(&gcam2);
  GCAMfree(&gcam_concat);
  exit(0);
}


static GCAM *concatenate(GCA_MORPH *gcam1, GCA_MORPH *gcam2)
{
  int x, y, z, width, height, depth, space;
  int out_of_gcam1;
  float xd, yd, zd, xdd, ydd, zdd, oxdd, oydd, ozdd;
  oxdd = oydd = ozdd = 0;
  GCA_MORPH_NODE *gcamnC, *gcamn2;

  printf("allocating gcam...(%d, %d, %d)\n", gcam2->width, gcam2->height,
                                                                  gcam2->depth);
  GCAM *gcamC = GCAMalloc(gcam2->width, gcam2->height, gcam2->depth);

  // copy gcam geometry information
  gcamC->image = gcam1->image;
  gcamC->atlas = gcam2->atlas;
  gcamC->spacing = space = gcam2->spacing;
  width = gcamC->width;
  height = gcamC->height;
  depth = gcamC->depth;

  for (x = 0 ; x < width ; x++){
    for (y = 0 ; y < height ; y++){
      for (z = 0 ; z < depth ; z++){

        gcamnC = &gcamC->nodes[x][y][z];
        gcamn2 = &gcam2->nodes[x][y][z];

        xd = gcamn2->x;
        yd = gcamn2->y;
        zd = gcamn2->z;

        gcamnC->xn = x;
        gcamnC->yn = y;
        gcamnC->zn = z;

        out_of_gcam1 = GCAMsampleMorph(gcam1, xd, yd, zd, &xdd, &ydd, &zdd);

        if (!out_of_gcam1) {
          GCAMsampleOrig(gcam1, (float)(x*space), (float)(y*space), (float)(z*space), &oxdd, &oydd, &ozdd);

          gcamnC->origx = oxdd;
          gcamnC->origy = oydd;
          gcamnC->origz = ozdd;
          gcamnC->x = xdd;
          gcamnC->y = ydd;
          gcamnC->z = zdd;
        }
        else {
          // This node is out of bounds with GCAM1 - its net value cannot be
          // estimated. To prevent it from being included in interpolation, 
          // it must be forced invalid by setting [xyz] and orig[xyz] to 0
          gcamnC->origx = 0;
          gcamnC->origy = 0;
          gcamnC->origz = 0;
          gcamnC->x = 0;
          gcamnC->y = 0;
          gcamnC->z = 0;
        }
      }
    }
  }

  return gcamC;
}


static int GCAMsampleOrig(const GCA_MORPH *gcam, float x, float y, float z, float *pxd, float *pyd, float *pzd)
{
  int            xm, xp, ym, yp, zm, zp, width, height, depth ;
  float          xmd, ymd, zmd, xpd, ypd, zpd ;  // d's are distances
  int            errCode = NO_ERROR;

  // x, y, z are in node coords
  x /= gcam->spacing;
  y /= gcam->spacing;
  z /= gcam->spacing;
  width = gcam->width;
  height = gcam->height;
  depth = gcam->depth;

  if ((errCode = boundsCheckf(x, y, z, width, height, depth)) != NO_ERROR) return errCode;

  xm = MAX((int)x, 0);
  xp = MIN(width-1, xm+1);
  ym = MAX((int)y, 0);
  yp = MIN(height-1, ym+1);
  zm = MAX((int)z, 0);
  zp = MIN(depth-1, zm+1);

  xmd = x - (float)xm;
  ymd = y - (float)ym;
  zmd = z - (float)zm;
  xpd = (1.0f - xmd);
  ypd = (1.0f - ymd);
  zpd = (1.0f - zmd);

  if (
    (gcam->nodes[xm][ym][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xm][ym][zp].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xm][yp][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xm][yp][zp].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][ym][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][ym][zp].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][yp][zm].invalid == GCAM_POSITION_INVALID) ||
    (gcam->nodes[xp][yp][zp].invalid == GCAM_POSITION_INVALID))
  {
    return ERROR_BADPARM;
  }

  *pxd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].origx +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].origx +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].origx +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].origx +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].origx +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].origx +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].origx +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].origx ;
  *pyd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].origy +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].origy +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].origy +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].origy +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].origy +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].origy +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].origy +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].origy ;
  *pzd =
    xpd * ypd * zpd * gcam->nodes[xm][ym][zm].origz +
    xpd * ypd * zmd * gcam->nodes[xm][ym][zp].origz +
    xpd * ymd * zpd * gcam->nodes[xm][yp][zm].origz +
    xpd * ymd * zmd * gcam->nodes[xm][yp][zp].origz +
    xmd * ypd * zpd * gcam->nodes[xp][ym][zm].origz +
    xmd * ypd * zmd * gcam->nodes[xp][ym][zp].origz +
    xmd * ymd * zpd * gcam->nodes[xp][yp][zm].origz +
    xmd * ymd * zmd * gcam->nodes[xp][yp][zp].origz ;

  return NO_ERROR;
}


static int boundsCheckf(float x, float y, float z, int width, int height, int depth)
{
  if (x >= width) {
    return ERROR_BADPARM;
  }
  else if (y >= height) {
    return ERROR_BADPARM;
  }
  else if (z >= depth) {
    return ERROR_BADPARM;
  }
  else if (x < 0) {
    return ERROR_BADPARM;
  }
  else if (y < 0) {
    return ERROR_BADPARM;
  }
  else if (z < 0) {
    return ERROR_BADPARM;
  }
  else {
    return NO_ERROR;
  }
}


static bool parseCommandLine(int argc, char *argv[], string &warp1name, string &warp2name, string &outname)
{
  int nargs, index;

  argc -= 1;
  argv += 1;
  index = 0;
  for (; argc > 0; argc--, argv++) {
    // parse options
    if (ISOPTION(*argv[0])) {
      nargs = parseOption(argc, argv);
      argc -= nargs;
      argv += nargs;
    }
    else {
      if (index == 0) {
        warp1name = string(argv[0]);
      }
      else if (index == 1) {
        warp2name = string(argv[0]);
      }
      else if (index == 2) {
        outname = string(argv[0]);
      }
      else if (index > 2) {
        printf("error: too many arguments!\n");
        return false;
      }
      index++;
    }
  }
  if (index <= 2) {
    printf("error: not enough arguments!\n");
    return false;
  }

  return true;
}


static int parseOption(int argc, char *argv[])
{
  int nargs = 0;
  char *option;

  // remove option dashes
  option = argv[0] + 1;
  if (option[0] == '-') {
    option = option + 1;
  }
  StrUpper(option);

  if (!strcmp(option, "HELP") || !strcmp(option, "H")) {
    usage();
    exit(1);
  }
  else {
    printf("%s is not a valid option\n", option);
    exit(1);
  }

  return nargs;
}


#include "mri_warp_concat.help.xml.h"
static void usage(void)
{
  outputHelpXml(mri_warp_concat_help_xml, mri_warp_concat_help_xml_len);
}
