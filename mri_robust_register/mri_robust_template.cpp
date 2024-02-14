/**
 * @brief combine multiple volumes by mean or median
 *
 * Creation of robust template of several volumes together with
 * their linear registration
 */

/*
 * Original Author: Martin Reuter
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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

// for rand:
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "MultiRegistration.h"

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

using namespace std;

#define SAT -1 // leave blank, either passed by user or --satit
//#define SAT 4.685 // this is suggested for gaussian noise
//#define SAT 20
#define SSAMPLE -1

// unsigned int printmemusage ()
// {
//   char buf[30];
//   snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
//   FILE* pf = fopen(buf, "r");
//   unsigned size = 0; //       total program size
//   if (pf) {
//     //unsigned resident;//   resident set size
//     //unsigned share;//      shared pages
//     //unsigned text;//       text (code)
//     //unsigned lib;//        library
//     //unsigned data;//       data/stack
//     //unsigned dt;//         dirty pages (unused in Linux 2.6)
//     //fscanf(pf, "%u" /* %u %u %u %u %u"*/, &size/*, &resident, &share, &text, &lib, &data*/);
//     fscanf(pf, "%u" , &size);
//     //DOMSGCAT(MSTATS, std::setprecision(4) << size / (1024.0) << "MB mem used");
//     cout <<  size / (1024.0) << " MB mem used" << endl;
//     fclose(pf);
//   }
// // while (1)
// // {
// //     if ('n' == getchar())
// //        break;
// // }
//   return size;
// }

struct Parameters
{
  vector<std::string> mov;
  vector<std::string> masks;
  string mean;
  vector<string> iltas;
  vector<string> nltas;
  vector<string> nweights;
  vector<string> mapmovhdr;
  bool fixvoxel;
  bool floattype;
  bool lta_vox2vox;
  bool affine;
  bool iscale;
  bool iscaleonly;
  bool transonly;
  bool leastsquares;
  bool nomulti;
  int iterate;
  double epsit;
  double sat;
  vector<string> nwarps;
  int debug;
  int average;
  int inittp;
  bool noit;
  bool quick;
  int subsamplesize;
  bool fixtp;
  bool satit;
  string conform;
  bool doubleprec;
  bool oneminusweights;
  vector<string> iscalein;
  vector<string> iscaleout;
  int finalinterp;
  int highit;
  unsigned int seed;
  bool crascenter;
  int pairiterate;
  double pairepsit;
  float  resthresh;
  double frobnormthresh;
  bool AllowDiffVoxSize;
};

// Initializations:
static struct Parameters P =
{ vector<string>(0), vector<string>(0), "", vector<string>(0), vector<string>(0), vector<string>(
    0), vector<string>(0), false, false, false, false, false, false, false, false, false,
    5, -1.0, SAT, vector<string>(0), 0, 1, -1, false, false, SSAMPLE, false, false, "", false,
  true, vector<string>(0), vector<string>(0), SAMPLE_CUBIC_BSPLINE, -1, 0 , false, 5, 0.01, 0.01, 0.0001,false};

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[], Parameters & P);

const char *Progname = NULL;

int getRandomNumber(int start, int end, unsigned int & seed)
// return n in [start,end]
// also return seed if it was -1
{

  if (seed == 0)
  {
    seed = time(NULL);
  }

  // initialize random seed:
  srand(seed);

  // generate random number: 
  int range = end - start + 1;

  return rand() % range + start;
}

int main(int argc, char *argv[])
{
  {
    cout << getVersion() << endl << endl;
    // set the environment variable
//  setenv("SURFER_FRONTDOOR","",1) ;
    // to store mri as chunk in memory:
//  setenv("FS_USE_MRI_CHUNK","",1) ;

    // Default initialization
    int nargs = handleVersionOption(argc, argv, "mri_robust_template");
    if (nargs && argc - nargs == 1)
    {
      exit(0);
    }
    argc -= nargs;
    Progname = argv[0];
    argc--;
    argv++;
    ErrorInit(NULL, NULL, NULL);
//  DiagInit(NULL, NULL, NULL) ;

    if (!parseCommandLine(argc, argv, P))
    {
      printUsage();
      exit(1);
    }
//  if (P.outdir[P.outdir.length()-1] != '/') P.outdir += "/";

    // Timer
    Timer start;

    int msec, minutes, seconds;
    start.reset();
    ///////////////////////////////////////////////////////////////

    MultiRegistration MR;

    // set parameters
    size_t pos = P.mean.rfind("/");    // position of "." in str
    if (pos != string::npos)
    {
      MR.setOutdir(P.mean.substr(0, pos - 1));
    }
    else
    {
      MR.setOutdir("./");
    }
    MR.setDebug(P.debug);
    MR.setRigid(!P.affine);
    MR.setTransonly(P.transonly);
    if (P.iscale) std::cout << "Setting iscale ..." << std::endl;
    MR.setIscale(P.iscale);
    if (P.iscaleonly) std::cout << "Setting iscaleonly ..." << std::endl;
    MR.setIscaleOnly(P.iscaleonly);
    MR.setNoMulti(P.nomulti);
    MR.setRobust(!P.leastsquares);
    MR.setSaturation(P.sat);
    MR.setSatit(P.satit);

    MR.setFixVoxel(P.fixvoxel);
    MR.setKeepType(!P.floattype);
    MR.setAverage(P.average);
    MR.setDoublePrec(P.doubleprec);
    MR.setSubsamplesize(P.subsamplesize);
    MR.setHighit(P.highit);
    if (P.nweights.size() > 0)
      MR.setBackupWeights(true);
    MR.useCRAS(P.crascenter);
    MR.setResthresh(P.resthresh);
    MR.SetFrobnormthresh(P.frobnormthresh);
    
    // init MultiRegistration and load movables
    //int nnin = (int) P.mov.size();
    //assert (P.mov.size() >1);
    //assert (MR.loadMovables(P.mov)==nin);
    int nin = MR.loadMovables(P.mov,P.masks,P.AllowDiffVoxSize);
    if (nin <=1)
    {
      std::cerr << "Could not load movables!" << std::endl;
      exit(1);
    }

    // load initial ltas if set:
    if (P.iltas.size() > 0)
    {
      const int numLoaded = MR.loadLTAs(P.iltas);
      if (numLoaded != nin)
      {
        std::cerr << "Error: num LTA files != num movables!" << std::endl;
	exit(1);
      }
    }

    // load initial iscales if set:
    if (P.iscalein.size() > 0)
    {
      const int numLoaded = MR.loadIntensities(P.iscalein);
      if (numLoaded != nin) 
      {
        std::cerr << "Error: num intensity files != num movables!" << std::endl;
	exit(1);
      }
    }

    // Randomly pick target (default):
    // Only if no inittp is passed and if no init transforms are given
    // In those cases we either fixed the target (--inittp)
    // or we won't need to compute ixforms when they are passed (--ixforms)
    if (P.inittp < 0 && P.iltas.size() == 0)
    {
      // compute seed based on input images if not set:
      if (P.seed == 0)
        P.seed = MR.getSeed();
      // inittp randomly
      P.inittp = getRandomNumber(1, P.mov.size(), P.seed);
      cout << "Will use TP " << P.inittp << " as random initial target (seed "
          << P.seed << " )." << endl;
    }

    if (P.noit) // no registration to mean space, only initial reg. and averaging
    {
      // if no initial xforms are given, use initialization to median space
      //   by registering everything first to inittp
      //   res 0: up to highest resolution, eps 0.01: accurate
      if (P.iltas.size() == 0 && P.inittp > 0)
      {
        MR.initialXforms(P.inittp, P.fixtp, 0, 5, 0.01);
      }

      // create template:
      MR.setSampleType(P.finalinterp);
      MR.mapAndAverageMov(0);
    }
    // run registrations
    else if (nin == 2)
    {

      if (P.iltas.size() == 0 && P.inittp > 0)
      {
        MR.initialXforms(P.inittp, P.fixtp, 0, 5, 0.01);
      }

      // we are done here, since with 2TP the registration
      // to the first mean space is allready accurate

      // create template:
      MR.setSampleType(P.finalinterp);
      MR.mapAndAverageMov(0);

//    // here default params are adjusted for just 2 images (if not passed):
//    if (P.iterate == -1) P.iterate = 5;
//    if (P.epsit <= 0)    P.epsit   = 0.01;
//    MR.halfWayTemplate(0,P.iterate,P.epsit,P.lta_vox2vox);
    }
    else
    {
      // if no initial xforms are given, use initialization to median space
      //   by registering everything first to inittp
      //   a) res 1: up to second highest resolution, eps 0.05: not too accurate
      //   b) res 0: up to highest res., eps 0.01 accurate reg.
      //   turns out accurate b) performs better and saves us
      //   from more global iterations (reg to template) later
      // remains open if subsampling speeds up things w/o increasing iterations later
      if (P.iltas.size() == 0 && P.inittp > 0)
      {
        MR.initialXforms(P.inittp, P.fixtp, 0, 5, 0.01);
      }
      //MR.initialXforms(P.inittp,1,5,0.05);

      // here default is adjusted for several images (and real mean/median target):
      if (P.iterate == -1)
      {
        P.iterate = 6;
      }
      if (P.epsit <= 0)
      {
        P.epsit = 0.03;
      }

      // P.iterate and P.epsit are used for terminating the template iterations
      // while 5 and 0.01 are default for the individual pairwise registrations
      MR.computeTemplate(P.iterate, P.epsit, P.pairiterate, P.pairepsit);
      // if final interp not set in multireg, re-average with new type
      if (P.finalinterp != MR.getSampleType())
      {
        MR.setSampleType(P.finalinterp);
        MR.mapAndAverageMov(0);
      }
    }

    cout << "Writing final template: " << P.mean << endl;
    MR.writeMean(P.mean);
    if (P.conform != "")
    {
      MR.writeConformMean(P.conform);
    }

    // output transforms and warps
    cout << "Writing final transforms (warps etc.)..." << endl;
    if (P.nltas.size() > 0)
    {
      MR.writeLTAs(P.nltas, P.lta_vox2vox, P.mean);
    }
    if (P.nwarps.size() > 0)
    {
      MR.writeWarps(P.nwarps);
    }
    if (P.mapmovhdr.size() > 0)
    {
      MR.writeMapMovHdr(P.mapmovhdr);
    }
    if (P.iscaleout.size() > 0)
    {
      MR.writeIntensities(P.iscaleout);
    }
    if (!P.leastsquares && P.nweights.size() > 0)
    {
      MR.writeWeights(P.nweights, P.oneminusweights);
    }

    MR.clear();

    ///////////////////////////////////////////////////////////////
    msec = start.milliseconds();
    seconds = nint((float) msec / 1000.0f);
    minutes = seconds / 60;
    seconds = seconds % 60;
    cout << "registration took " << minutes << " minutes and " << seconds
        << " seconds." << endl;
    //if (diag_fp) fclose(diag_fp) ;

    cout << endl << " Thank you for using RobustTemplate! " << endl;
    cout << " If you find it useful and use it for a publication, please cite: "
        << endl << endl;
    cout
        << " Within-Subject Template Estimation for Unbiased Longitudinal Image Analysis"
        << endl;
    cout
        << " M. Reuter, N.J. Schmansky, H.D. Rosas, B. Fischl. NeuroImage 2012."
        << endl;
    cout << " http://dx.doi.org/10.1016/j.neuroimage.2012.02.084" << endl;
    cout << " http://reuter.mit.edu/papers/reuter-long12.pdf" << endl << endl;

  }
  exit(0);
  return (0);
}

/*----------------------------------------------------------------------
 ----------------------------------------------------------------------*/
#include "mri_robust_template.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_robust_template_help_xml, mri_robust_template_help_xml_len);
}

/*!
 \fn int parseNextCommand(int argc, char **argv)
 \brief Parses the command-line for next command
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       number of used arguments for this command
 */
static int parseNextCommand(int argc, char *argv[], Parameters & P)
{
  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  StrUpper(option);

  //cout << " option: " << option << endl;

  if (!strcmp(option, "MOV"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.mov.push_back(string(argv[nargs]));
        cout << "--mov: Using " << P.mov.back() << " as movable/source volume."
            << endl;
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "    Total: " << nargs << " input volumes" << endl;
  }
//  else if (!strcmp(option, "OUTDIR") )
//  {
//     P.outdir = string(argv[1]);
//     nargs = 1;
//     cout << "--outdir: Using "<< P.outdir << " as output directory." << endl;
//  }
  else if (!strcmp(option, "TEMPLATE"))
  {
    P.mean = string(argv[1]);
    nargs = 1;
    cout << "--template: Using " << P.mean << " as template output volume."
        << endl;
  }
  else if (!strcmp(option, "LTA"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.nltas.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nltas.back() << " as LTA." << endl;
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--lta: Will output LTA transforms" << endl;
  }
  else if (!strcmp(option, "MASKS"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.masks.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nltas.back() << " as LTA." << endl;
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--masks: Will use masks on inputs." << endl;
  }
  else if (!strcmp(option, "ISCALEOUT"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.iscaleout.push_back(string(argv[nargs]));
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    P.iscale = true;
    cout << "--iscaleout: Will perform intensity scaling and output results"
        << endl;
  }
  else if (!strcmp(option, "ISCALEIN"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.iscalein.push_back(string(argv[nargs]));
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--iscalein: Will use init intensity scales" << endl;
  }
  else if (!strcmp(option, "IXFORMS"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;

        P.iltas.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nltas.back() << " as LTA." << endl;
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--ixforms: Will use init XFORMS." << endl;
  }
  else if (!strcmp(option, "AVERAGE"))
  {
    P.average = atoi(argv[1]);
    nargs = 1;
    cout << "--average: Using method " << P.average
        << " for template computation." << endl;
  }
  else if (!strcmp(option, "VOX2VOX"))
  {
    P.lta_vox2vox = true;
    cout << "--vox2vox: Output transforms as VOX2VOX. " << endl;
  }
  else if (!strcmp(option, "AFFINE") || !strcmp(option, "A"))
  {
    P.affine = true;
    cout << "--affine: Enabling affine transform!" << endl;
  }
  else if (!strcmp(option, "ISCALE") || !strcmp(option, "I"))
  {
    P.iscale = true;
    cout << "--iscale: Enabling intensity scaling!" << endl;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    P.transonly = true;
    cout << "--transonly: Using only translation!" << endl;
  }
  else if (!strcmp(option, "ISCALEONLY"))
  {
    P.iscale = true;
    P.iscaleonly = true;
    P.transonly = false;
    P.affine = false;
    cout << "--iscaleonly: Computing only global scaling!" << endl;
  }
  else if (!strcmp(option, "NOMULTI"))
  {
    P.nomulti = true;
    cout << "--nomulti: process highest resolution only!" << endl;
  }
  else if (!strcmp(option, "LEASTSQUARES") || !strcmp(option, "L"))
  {
    P.leastsquares = true;
    cout << "--leastsquares: Using standard least squares (non-robust)!"
        << endl;
  }
  else if (!strcmp(option, "MAXIT"))
  {
    P.iterate = atoi(argv[1]);
    nargs = 1;
    cout << "--maxit: Performing maximal " << P.iterate
        << " iterations on each resolution." << endl;
  }
  else if (!strcmp(option, "HIGHIT"))
  {
    P.highit = atoi(argv[1]);
    nargs = 1;
    cout << "--highit: Performing maximal " << P.highit
        << " iterations on highest resolution" << endl;
  }
  else if (!strcmp(option, "EPSIT"))
  {
    P.epsit = atof(argv[1]);
    nargs = 1;
    cout << "--epsit: Stop iterations when change is less than " << P.epsit
        << " . " << endl;
  }
  else if (!strcmp(option, "PAIRMAXIT"))
  {
    P.pairiterate = atoi(argv[1]);
    nargs = 1;
    cout << "--pairmaxit: Performing maximal " << P.pairiterate
        << " iterations on each resolution for individual pairwise registrations." << endl;
  }
  else if (!strcmp(option, "PAIREPSIT"))
  {
    P.pairepsit = atof(argv[1]);
    nargs = 1;
    cout << "--pairepsit: Stop individual pairwise iterations when change is less than " << P.pairepsit
        << " . " << endl;
  }
  else if (!strcmp(option, "SAT"))
  {
    P.sat = atof(argv[1]);
    nargs = 1;
    cout << "--sat: Using saturation " << P.sat << " in M-estimator!" << endl;
  }
  else if (!strcmp(option, "SUBSAMPLE"))
  {
    P.subsamplesize = atoi(argv[1]);
    nargs = 1;
    if (P.subsamplesize >= 0)
    {
      cout << "--subsample: Will subsample if size is larger than "
          << P.subsamplesize << " on all axes!" << endl;
    }
    else
    {
      cout << "--subsample -1: Will not subsample on any scale!" << endl;
    }
  }
  else if (!strcmp(option, "DEBUG"))
  {
    P.debug = 1;
    nargs = 0;
    cout << "--debug: Will output debug info and files!" << endl;
  }
  else if (!strcmp(option, "NOIT"))
  {
    P.noit = true;
    nargs = 0;
    cout << "--noit: Will output only first template (no iterations)!" << endl;
  }
  else if (!strcmp(option, "FIXTP"))
  {
    P.fixtp = true;
    nargs = 0;
    cout << "--fixtp: Will map everything to init TP!" << endl;
  }
  else if (!strcmp(option, "ALLOW-DIFF-VOX-SIZE"))
  {
    P.AllowDiffVoxSize = true;
    nargs = 0;
    cout << "--allow-diff-vox-size " << endl;
  }
  else if (!strcmp(option, "SATIT"))
  {
    P.satit = true;
    nargs = 0;
    cout << "--satit: Will estimate SAT iteratively!" << endl;
  }
  else if (!strcmp(option, "DOUBLEPREC"))
  {
    P.doubleprec = true;
    nargs = 0;
    cout
        << "--doubleprec: Will perform algorithm with double precision (higher mem usage)!"
        << endl;
  }
  else if (!strcmp(option, "WEIGHTS"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.nweights.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nweights.back() << " as weights volume." << endl;
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--weights: Will output weights in target space" << endl;
  }
  else if (!strcmp(option, "WARP") || !strcmp(option, "MAPMOV"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.nwarps.push_back(string(argv[nargs]));
        //cout << "Using "<< P.nwarps.back() << " as weights volume." << endl;
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--mapmov: Will save mapped movables/sources!" << endl;
  }
  else if (!strcmp(option, "MAPMOVHDR"))
  {
    nargs = 0;
    do
    {
      option = argv[nargs + 1];
      if (option[0] != '-')
      {
        nargs++;
        P.mapmovhdr.push_back(string(argv[nargs]));
      }
    } while (nargs + 1 < argc && option[0] != '-');
    assert(nargs > 0);
    cout << "--mapmovhdr: Will save header-adjusted movables!" << endl;
  }
  else if (!strcmp(option, "TEST"))
  {
    cout << "--test: TEST-MODE " << endl;
    RegRobust R;
    R.testRobust(argv[2], atoi(argv[1]));
    nargs = 2;
    exit(0);
  }
  else if (!strcmp(option, "CONFORM"))
  {
    P.conform = argv[1];
    nargs = 1;
    cout << "--conform: Will output conform template (256^3 and 1mm voxels)!"
        << endl;
  }
//   else if (!strcmp(option, "CONFORM") )
//   {
//     P.fixvoxel = true;
//     nargs = 0 ;
//     cout << "Will conform images to 256^3 and voxels to 1mm!" << endl;
//   }
  else if (!strcmp(option, "FLOATTYPE"))
  {
    P.floattype = true;
    nargs = 0;
    cout << "--floattype: Use float images internally (independent of input)!"
        << endl;
  }
  else if (!strcmp(option, "FINALNEAREST"))
  {
    P.finalinterp = SAMPLE_NEAREST;
    nargs = 0;
    cout
        << "--finalnearest: Use nearest neighbor interpolation for final average!"
        << endl;
  }
  else if (!strcmp(option, "INITTP"))
  {
    P.inittp = atoi(argv[1]);
    nargs = 1;
    if (P.inittp == 0)
    {
      cout
          << "--inittp 0: No initialization, construct first mean from original TPs"
          << endl;
    }
    else
    {
      cout << "--inittp: Using TP " << P.inittp
          << " as target for initialization" << endl;
    }
  }
  else if (!strcmp(option, "ONEMINUSW"))
  {
    P.oneminusweights = false;
    nargs = 0;
    cout
        << "--oneminusw: Will output 1-weights (zero=outlier), as in earlier versions!"
        << endl;
  }
  else if (!strcmp(option, "SEED"))
  {
    P.seed = atoi(argv[1]);
    nargs = 1;
    cout << "--seed: Will use random seed " << P.seed << endl;
  }
  else if (!strcmp(option, "CRAS"))
  {
    P.crascenter = true;
    nargs = 0;
    cout << "--cras: Will center template at avgerage CRAS!" << endl;
  }
  else if (!stricmp(option, "HELP") || !stricmp(option, "USAGE")
      || !stricmp(option, "h") || !stricmp(option, "u"))
  {
    printUsage();
    exit(0);
  }
  else if (!strcmp(option, "RES-THRESH"))
  {
    P.resthresh = atof(argv[1]);
    nargs = 1;
    cout << "--res-thresh: Voxsize threshold " << P.resthresh << endl;
  }
  else if (!strcmp(option, "FROBNORM-THRESH"))
  {
    P.frobnormthresh = atof(argv[1]);
    nargs = 1;
    cout << "--frobnorm-thresh: Matrix frobenius norm threshold " << P.frobnormthresh << endl;
  }
  else if (!strcmp(option, "THREADS"))
  {
    nargs = 1;
#ifdef _OPENMP
    // Does not speed things up much
    omp_set_num_threads(atoi(argv[1]));
#endif
  }
  else
  {
    cerr << endl << endl << "ERROR: Option " << argv[0] << " unknown !!! "
        << endl << endl;
    printUsage();
    exit(1);
  }

  fflush(stdout);

  return (nargs);
}

/*!
 \fn int parseCommandLine(int argc, char **argv)
 \brief Parses the command-line
 \param   argc  number of command line arguments
 \param   argv  pointer to a character pointer
 \param      P  reference to parameters
 \returns       if all necessary parameters were set
 */
static bool parseCommandLine(int argc, char *argv[], Parameters & P)
{
  int nargs;
  int n = argc;

  if (n == 0)
  {
    printUsage();
    exit(1);
  }

  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv, P);
    argc -= nargs;
    argv += nargs;
  }

  bool ntest = true;

  if (P.mov.size() < 1)
  {
    ntest = false;
    cerr << "ERROR: Specify at least 2 inputs with --mov (or one 4D input)!"
        << endl;
    exit(1);
  }
  if (P.mean == "")
  {
    ntest = false;
    cerr << "ERROR: Specify output image with --template !" << endl;
    exit(1);
  }
  if (P.nwarps.size() > 0 && P.mov.size() != P.nwarps.size())
  {
    ntest = false;
    cerr
        << "ERROR: Number of filenames for --warp should equal number of inputs!"
        << endl;
    exit(1);
  }
  if (P.mapmovhdr.size()>0 && P.mov.size()!=P.mapmovhdr.size())
  {
    ntest = false;
      std::cerr
      << "ERROR: Number of filenames for --mapmovhdr should equal number of inputs!"
      << std::endl;
    exit(1);
  }
  if (P.nltas.size() > 0 && P.mov.size() != P.nltas.size())
  {
    ntest = false;
    cerr
        << "ERROR: Number of filenames for --lta should equal number of inputs!"
        << endl;
    exit(1);
  }
  if (P.iltas.size() > 0 && P.mov.size() != P.iltas.size())
  {
    ntest = false;
    cerr
        << "ERROR: Number of filenames for --ixforms should equal number of inputs!"
        << endl;
    exit(1);
  }
  if (P.nweights.size() > 0 && P.mov.size() != P.nweights.size())
  {
    ntest = false;
    cerr
        << "ERROR: Number of filenames for --weights should equal number of inputs!"
        << endl;
    exit(1);
  }
  if (P.iscalein.size() > 0 && P.iltas.size() != P.iscalein.size())
  {
    ntest = false;
    cerr
        << "ERROR: Number of filenames for --iscalein should agree with number of init LTAs (--ixforms)!"
        << endl;
    exit(1);
  }
  if (P.iscaleout.size() > 0 && P.mov.size() != P.iscaleout.size())
  {
    ntest = false;
    cerr
        << "ERROR: Number of filenames for --iscaleout should equal number of inputs!"
        << endl;
    exit(1);
  }
  if (n > 0 && !P.satit && P.sat <= 0 && !(P.noit && P.iltas.size() > 0))
  {
    ntest = false;
    cerr
        << "ERROR: Specify either --satit or set sensitivity manually with --sat <real> !"
        << endl;
    exit(1);
  }

  return (ntest);
}
