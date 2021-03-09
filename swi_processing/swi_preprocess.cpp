/**
 * @brief Pre-processing for Susceptibility Weighted Imaging
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
#include "cmd_line_interface.h"
const char *Progname;

// static function declarations
// forward declaration
struct IoParams;
using namespace std;

/*
  IO structure
*/
struct IoParams
{
  std::string scanner;
  std::string gefile;
  std::string philipsfile;
  std::string magfile;
  std::string phasefile;
  std::string outmag;
  std::string outphase;

  IoParams();
  void parse(int ac, char* av[]);
};

int main(int argc, char*argv[])
{
  IoParams params;
  try
  {
    params.parse(argc,argv);
    if ( params.scanner == "siemens" && ( params.gefile != "" || params.magfile == "" || params.phasefile == "" ) )
    {
      throw std::logic_error("Siemens scanner should only have siemens mag and siemens phase options. See help");
    }
    if ( params.scanner == "ge" && ( params.gefile == "" || params.magfile != "" || params.phasefile != "" ) )
    {
      throw std::logic_error("GE scanner should only have a ge-file option in the command line. See help");
    }
  }
  catch (std::exception& excp)
  {
    std::cerr << "ERROR: Exception caught while parsing the command-line\n"
              << excp.what() << std::endl;
    exit(1);
  }

  std::cout << "Reading the input file(s)\n";
  MRI *mrimag=NULL, *mriphase=NULL;

  if ( params.scanner == "siemens" )
  {
    mrimag = DICOMRead2(params.magfile.c_str(), TRUE );
    if ( mrimag == NULL )
    {
      std::cerr << "ERROR: The Siemens magnitude DICOM can't be read. Check the file\n";
      exit(1);
    }
    MRI* _sphase = DICOMRead2(params.phasefile.c_str(), TRUE );
    if ( _sphase == NULL )
    {
      std::cerr << "ERROR: The Siemens phase DICOM can't be read. Check the file\n";
      exit(1);
    }
    //mriphase = MRIcloneBySpace(_sphase, MRI_FLOAT, 1);
    /*float vmin, vmax;
    MRIvalRange(_sphase, &vmin, &vmax);
    printf("%f, %f\n", vmin, vmax);*/
    mriphase = MRIchangeType(_sphase, MRI_FLOAT, -M_PI, M_PI, 0);
    MRIvalScale(mriphase, mriphase, -M_PI, M_PI);


    MRIfree(&_sphase);

  }
  else if ( params.scanner == "ge" )
  {
    MRI *_gedcm = NULL;
    _gedcm      = DICOMRead2(params.gefile.c_str(), TRUE);
    if ( _gedcm == NULL )
    {
      std::cerr << "ERROR: The GE DICOM can't be read. Check the file\n";
      exit(1);
    }
    mrimag = MRIcloneBySpace(_gedcm, MRI_FLOAT, 1);
    mriphase = MRIcloneBySpace(_gedcm, MRI_FLOAT, 1);
    float _r, _i;
    for ( int i=0; i < _gedcm->width; i++)
      for ( int j=0; j < _gedcm->height; j++)
        for ( int k=0; k < _gedcm->depth; k++)
        {
          _r = (float)(MRIgetVoxVal(_gedcm, i, j, k, 1)) ;
          _i = (float)(MRIgetVoxVal(_gedcm, i, j, k, 2)) ;
          MRIsetVoxVal(mrimag, i, j, k, 0, (float)MRIgetVoxVal(_gedcm, i, j, k, 0));
          //MRIsetVoxVal(mrimag, i, j, k, 0, sqrt((_i/1000.0)*(_i/1000.0) + (_r/1000.0)*(_r/1000.0)));
          MRIsetVoxVal(mriphase, i, j, k, 0, atan2(_i, _r) );
        }
    MRIfree(&_gedcm);
  }
  else if ( params.scanner == "philips" )
  {
    MRI *_phildcm = NULL;
    _phildcm      = DICOMRead2(params.philipsfile.c_str(), TRUE);
    if ( _phildcm == NULL )
    {
      std::cerr << "ERROR: The Philips DICOM can't be read. Check the file\n";
      exit(1);
    }
    mrimag = MRIcloneBySpace(_phildcm, MRI_FLOAT, 1);
    MRI *_pphase = MRIcloneBySpace(_phildcm, MRI_FLOAT, 1);
    for ( int i=0; i < _phildcm->width; i++)
      for ( int j=0; j < _phildcm->height; j++)
        for ( int k=0; k < _phildcm->depth; k++)
        {
          MRIsetVoxVal(mrimag, i, j, k, 0, (float)MRIgetVoxVal(_phildcm, i, j, k, 0));
          MRIsetVoxVal(_pphase, i, j, k, 0, (float)MRIgetVoxVal(_phildcm, i, j, k, 1));
        }
    mriphase = MRIchangeType(_pphase, MRI_FLOAT, -M_PI, M_PI, 0);
    MRIvalScale(mriphase, mriphase, -M_PI, M_PI);
    MRIfree(&_phildcm);
    MRIfree(&_pphase);
  }
  else
  {
    std::cerr << "ERROR: You didn't specify the scanner make ( one of ge, siemens and philips )\n";
    exit(1);
  }

  std::cout << "Writing out the Magnitude and Phase files\n";
  MRIwrite(mrimag, params.outmag.c_str());
  MRIwrite(mriphase, params.outphase.c_str());
  MRIfree(&mrimag);
  MRIfree(&mriphase);
}


IoParams::IoParams()
{

}

void IoParams::parse(int argc, char* argv[])
{
  CCmdLineInterface cmd(argv[0]);

  cmd.AddOptionString("scanner", &scanner, "Name of the scanner ( one of ge, siemens or philips )");
  cmd.AddOptionString("ge_file", &gefile, "Name of the input GE file (only compatible with --scanner ge option)");
  cmd.AddOptionString("philips_file", &philipsfile, "Name of the input Philips file (only compatible with --scanner philips option)");
  cmd.AddOptionString("siemens_mag", &magfile, "Name of the input Siemens magnitude file (only compatible with --scanner siemens option)");
  cmd.AddOptionString("siemens_phase", &phasefile, "Name of the input Siemens phase file (only compatible with --scanner siemens option)");
  cmd.AddOptionString("out_magnitude", &outmag, "Name of the output magnitude file (after preprocessing). Make sure it has .nii suffix");
  cmd.AddOptionString("out_phase", &outphase, "Name of the output phase file (after preprocessing). Make sure it has .nii suffix");
  if ( argc == 1 )
  {
    std::cout << "\n"
              " Pre-process the Susceptibility-weighted images. Write out nifti files so that they can be fed into PRELUDE ( Phase Unwarpping Library of FSL).\n";
    cmd.PrintHelp();
    exit(0);
  }

  cmd.Parse(argc, argv);

}

