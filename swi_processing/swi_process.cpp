/**
 * @brief Processing for Susceptibility Weighted Imaging Pipeline
 *
 */
/*
 * Original Author: Krish Subramaniam
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

#include "fsenv.h"
#include "mri.h"
#include "diag.h"
#include "DICOMRead.h"
#include <stdexcept>
#include <vector>
#include <iterator>
#include <algorithm>
#include "cmd_line_interface.h"
const char *Progname;

// static function declarations
// forward declaration
struct IoParams;
using namespace std;
#include <limits>       // std::numeric_limits

/*
  IO structure which store the command line args
*/
struct IoParams
{
  std::string magfile;
  std::string phasefile;
  std::string swioutput;
  float stddev;
  std::string phasemethod;
  int phasemult;
  float phasecutoff;
  float phasercutoff;
  float sigmoid_a;
  float sigmoid_b;
  int miplevel;

  IoParams();
  void parse(int ac, char* av[]);
};

void do_phase_mask(IoParams& params, MRI* mag, MRI* pha, MRI* out)
{

  float mval, pval, outval, _ov =0.0, v;

  if (params.phasemethod.compare("negative") == 0)
  {
    std::cout << "Performing phase mask multiplication with cutoff " << params.phasecutoff << " ..\n";
    std::cout << "Number of phase multiplications " << params.phasemult << " ..\n";
    for ( int i=0; i < out->width; i++)
      for ( int j=0; j < out->height; j++)
        for ( int k=0; k < out->depth; k++)
        {
          mval = MRIgetVoxVal(mag, i, j, k, 0);
          pval = MRIgetVoxVal(pha, i, j, k, 0);
          // do phase cut off
          if ( pval > 0.0 )
          {
            _ov = 1.0;
          }
          else if ( pval <= params.phasecutoff )
          {
            _ov = 0.0;
          }
          else
          {
            _ov = 1.0 - (pval / params.phasecutoff );
          }

          // do phase multiplications
          outval=1.0;
          for (int phm=0; phm < params.phasemult; phm++)
          {
            outval = outval * _ov;
          }
          v = outval * mval;
          MRIsetVoxVal(out, i, j, k, 0, v);
        }
  }

  if (params.phasemethod.compare("positive") == 0)
  {
    std::cout << "Performing phase mask multiplication with right cutoff " << params.phasercutoff << " ..\n";
    std::cout << "Number of phase multiplications " << params.phasemult << " ..\n";
    for ( int i=0; i < out->width; i++)
      for ( int j=0; j < out->height; j++)
        for ( int k=0; k < out->depth; k++)
        {
          mval = MRIgetVoxVal(mag, i, j, k, 0);
          pval = MRIgetVoxVal(pha, i, j, k, 0);
          // do phase cut off
          if ( pval < 0.0 )
          {
            _ov = 1.0;
          }
          else if ( pval >= params.phasercutoff )
          {
            _ov = 0.0;
          }
          else
          {
            _ov =  1.0 - (pval / params.phasercutoff ) ;
          }

          // do phase multiplications
          outval=1.0;
          for (int phm=0; phm < params.phasemult; phm++)
          {
            outval = outval * _ov;
          }
          v = outval * mval;
          MRIsetVoxVal(out, i, j, k, 0, v);
        }
  }

  if (params.phasemethod.compare("symmetric") == 0)
  {
    std::cout << "Performing phase mask multiplication with cutoff " << params.phasecutoff << " ..\n";
    std::cout << "Performing phase mask multiplication with right cutoff " << params.phasercutoff << " ..\n";
    std::cout << "Number of phase multiplications " << params.phasemult << " ..\n";
    if ( fabs( fabs(params.phasecutoff) - fabs(params.phasercutoff)) > 0.0001 )
    {
      std::cerr << "ERROR: symmetric option implies both left cutoff and right cutoff are almost same\n";
      exit(1);
    }
    for ( int i=0; i < out->width; i++)
      for ( int j=0; j < out->height; j++)
        for ( int k=0; k < out->depth; k++)
        {
          mval = MRIgetVoxVal(mag, i, j, k, 0);
          pval = MRIgetVoxVal(pha, i, j, k, 0);
          // do phase cut off
          if ( pval < params.phasecutoff )
          {
            _ov = 0.0;
          }
          else if ( pval >= params.phasercutoff )
          {
            _ov = 0.0;
          }
          else if ( pval >= 0.0 )
          {
            _ov =  1.0 - (pval / params.phasercutoff ) ;
          }
          else
          {
            _ov =  1.0 - (pval / params.phasecutoff ) ;
          }

          // do phase multiplications
          outval=1.0;
          for (int phm=0; phm < params.phasemult; phm++)
          {
            outval = outval * _ov;
          }
          v = outval * mval;
          MRIsetVoxVal(out, i, j, k, 0, v);
        }
  }

  if (params.phasemethod.compare("asymmetric") == 0)
  {
    std::cout << "Performing phase mask multiplication with cutoff " << params.phasecutoff << " ..\n";
    std::cout << "Performing phase mask multiplication with right cutoff " << params.phasercutoff << " ..\n";
    std::cout << "Number of phase multiplications " << params.phasemult << " ..\n";
    for ( int i=0; i < out->width; i++)
      for ( int j=0; j < out->height; j++)
        for ( int k=0; k < out->depth; k++)
        {
          mval = MRIgetVoxVal(mag, i, j, k, 0);
          pval = MRIgetVoxVal(pha, i, j, k, 0);
          // do phase cut off
          if ( pval < params.phasecutoff )
          {
            _ov = 0.0;
          }
          else if ( pval >= params.phasercutoff )
          {
            _ov = 1.0;
          }
          else if ( pval >= 0.0 )
          {
            _ov =  1.0 - (pval / params.phasercutoff ) ;
          }
          else
          {
            _ov =  1.0 - (pval / params.phasecutoff ) ;
          }

          // do phase multiplications
          outval=1.0;
          for (int phm=0; phm < params.phasemult; phm++)
          {
            outval = outval * _ov;
          }
          v = outval * mval;
          MRIsetVoxVal(out, i, j, k, 0, v);
        }
  }

  if (params.phasemethod.compare("sigmoid") == 0)
  {
    std::cout << "Number of phase multiplications " << params.phasemult << " ..\n";
    float a = params.sigmoid_a;
    float b = params.sigmoid_b;
    for ( int i=0; i < out->width; i++)
      for ( int j=0; j < out->height; j++)
        for ( int k=0; k < out->depth; k++)
        {
          mval = MRIgetVoxVal(mag, i, j, k, 0);
          pval = MRIgetVoxVal(pha, i, j, k, 0);
          if ( pval < params.phasecutoff )
          {
            _ov = 0.0;
          }
          else if ( pval >= params.phasercutoff )
          {
            _ov = 1.0;
          }
          else
          {
            _ov = 1.0 / ( 1.0 + exp (-a * ( pval - b ) ) );
          }

          // do phase multiplications
          outval=1.0;
          for (int phm=0; phm < params.phasemult; phm++)
          {
            outval = outval * _ov;
          }
          v = outval * mval;
          MRIsetVoxVal(out, i, j, k, 0, v);
        }
  }


}

int main(int argc, char*argv[])
{
  IoParams params;
  try
  {
    params.parse(argc,argv);
  }
  catch (std::exception& excp)
  {
    std::cerr << "ERROR: Exception caught while parsing the command-line\n"
              << excp.what() << std::endl;
    exit(1);
  }

  std::cout << "Reading the input file(s)\n";
  MRI *mrimag=NULL, *mriphase=NULL;

  mrimag = MRIread(params.magfile.c_str() );
  std::cout << "Read the files..\n";
  if ( mrimag == NULL )
  {
    std::cerr << "ERROR: The magnitude image can't be read. Check the file\n";
    exit(1);
  }
  mriphase = MRIread(params.phasefile.c_str() );
  if ( mriphase == NULL )
  {
    std::cerr << "ERROR: The phase image can't be read. Check the file\n";
    exit(1);
  }

  if  ((params.phasemethod.compare("negative") != 0 ) &&
       (params.phasemethod.compare("positive") != 0 ) &&
       (params.phasemethod.compare("symmetric") != 0 ) &&
       (params.phasemethod.compare("asymmetric") != 0 ) &&
       (params.phasemethod.compare("sigmoid") != 0 ))
  {
    std::cerr << "ERROR: The phase mask method has to be one of negative, positive, symmetric, asymmetric or sigmoid\n";
  }

  float vmin, vmax;
  MRI *outimg = MRIclone(mriphase, NULL);
  // the following are the various intermediary images
  MRI *_hpimg = MRIclone(mriphase, NULL);
  MRI *_tmpimg = MRIclone(mriphase, NULL);

  std::cout << "Performing Gaussian smoothing with stddev " << params.stddev << " ..\n";
  MRIgaussianSmooth(mriphase, params.stddev, 1, _tmpimg);
  //MRIwrite(_tmpimg, "gaussian.mgz");

  std::cout << "Subtracting the smoothed image from the phase image..\n";
  MRIsubtract(mriphase, _tmpimg, _tmpimg);
  /*MRI_fft_highpass(mriphase, _tmpimg, 20);*/
  //MRIwrite(_tmpimg, "highpass.mgz");
  //
  MRIvalRange(_tmpimg, &vmin, &vmax);
  std::cout << "Min and max of the highpass image " << vmin << " " << vmax << "\n";

  if ( params.phasecutoff < vmin )
  {
    params.phasecutoff = vmin;
  }
  if (params.phasercutoff > vmax )
  {
    params.phasercutoff = vmax;
  }

  // perform phase masking operation. Result is in _hpimg
  do_phase_mask(params, mrimag, _tmpimg, _hpimg);

  // minimum Intensity Projection along y
  std::cout << "Performing Minimum Intensity Projection along y direction with levels=" << params.miplevel << " ..\n";
  int tmp_idx;
  std::vector<float> vals(params.miplevel);
  float mip_val;
  for ( int i=0; i < _hpimg->width; i++)
    for ( int j=0; j < _hpimg->height; j++)
      for ( int k=0; k < _hpimg->depth; k++)
      {
        for ( int jc=j; jc < j+params.miplevel; jc++)
        {
          tmp_idx = jc;
          // bounds check
          if ( jc < 0 )
          {
            tmp_idx = 0;
          }
          if ( jc >= _hpimg->height-1 )
          {
            tmp_idx = _hpimg->height - 1;
          }
          vals.push_back( MRIgetVoxVal(_hpimg, i, tmp_idx, k, 0));
        }

        mip_val = *(std::min_element( vals.begin(), vals.end() ));
        vals.clear();
        MRIsetVoxVal(outimg, i, j, k, 0, mip_val);
      }
  std::cout << "Writing the swi processed image\n";
  MRIwrite(outimg, params.swioutput.c_str() );
  // Freeing all
  MRIfree(&_hpimg);
  MRIfree(&mrimag);
  MRIfree(&mriphase);
  MRIfree(&outimg);
  MRIfree(&_tmpimg);
}

// These are the default values for the command line arguments
IoParams::IoParams()
{
  stddev      = 2.0;
  phasemethod = "negative";
  phasecutoff = -std::numeric_limits<float>::max();
  phasercutoff = std::numeric_limits<float>::max();
  sigmoid_a   = 1.0;
  sigmoid_b   = 0.0;
  phasemult   = 4;
  miplevel    = 4;
}

void IoParams::parse(int argc, char* argv[])
{
  CCmdLineInterface cmd(argv[0]);

  cmd.AddOptionString(
    "mag_file",
    &magfile,
    "The magnitude image ( Output from the PRELUDE program)");
  cmd.AddOptionString(
    "phase_file",
    &phasefile,
    "The phase image ( Output from the PRELUDE program)");
  cmd.AddOptionString(
    "swi_output",
    &swioutput,
    "Name of the SWI processed output image");
  cmd.AddOptionFloat(
    "stddev",
    &stddev,
    (char*)"Specify the standard deviation of the Gaussian Smoothing Filter. Default is 2.0");
  cmd.AddOptionString(
    "phase_mask_method",
    &phasemethod,
    "Specify the phase mask method. One of negative, positive, symmetric, asymmetric, sigmoid. Default is negative");
  cmd.AddOptionInt(
    "phase_multiplications",
    &phasemult, 
    "Specify the number of phase multiplications. Default is 4");
  cmd.AddOptionFloat(
    "phase_mask_cutoff",
    &phasecutoff,
    (char*)"Specify the negative phase mask cutoff frequency ( in radianѕ). Default is the minimum value of the phase image.");
  cmd.AddOptionFloat(
    "phase_mask_right_cutoff",
    &phasercutoff,
    (char*)"Specify the positive phase mask cutoff frequency ( in radianѕ). Default is the maximum value of the phase image.");
  cmd.AddOptionFloat(
    "sigmoid_a",
    &sigmoid_a,
    (char*)"Specify 'a' for the sigmoid formula f(phase)=1/(1+exp(-a*(phase-b))). Default is 1.0. Meaningless with phase_method != sigmoid");
  cmd.AddOptionFloat(
    "sigmoid_b",
    &sigmoid_b,
    (char*)"Specify 'b' for the sigmoid formula f(phase)=1/(1+exp(-a*(phase-b))). Default is 0.0. Meaningless with phase_method != sigmoid");
  cmd.AddOptionInt(
    "mip_level",
    &miplevel,
    "Specify the number of levels of mIP across the y direction. Default is 4");
  if ( argc == 1 )
  {
    std::cout << "\n"
              " Process the Susceptibility-weighted images. Make sure the inputs to this program is after the phase unwrapping step using PRELUDE\n";
    cmd.PrintHelp();
    exit(0);
  }

  cmd.Parse(argc, argv);

}


