/*

Gheorghe Postelnicu, 2006

This will load a morph and 2 volumes
- a reference and a moving

as well as a transform.

It will apply the transform to the volume, while caching it, 
then it will export it to the gcam structure.

As a validation, it will use the gcam to morph again the volumes

*/

// STL
#include <iostream>
#include <string>

#include "morph.h"

// FS
#include "argparse.h"
#include "timer.h"

#include "mri.h"
#include "gcamorph.h"


#include "exportGcam.help.xml.h"

// required by FS
const char* Progname;

struct IoParams
{
  std::string strFixed;
  std::string strMoving;
  std::string strGcam;
  std::string strMorph;
  std::string strInterp;

  bool useBoundingBox;
  int  threshold;
  unsigned int zlibBuffer;
  int doTest;

  void parse(int ac, const char** av);
};

char* getChar(const std::string& str)
{
  return const_cast<char*>( str.c_str() );
}

void initOctree2( gmp::VolumeMorph& morph)
{
  using namespace gmp;

  for ( VolumeMorph::TransformContainerType::iterator
        it = morph.m_transforms.begin();
        it != morph.m_transforms.end();
        ++it )
  {
    (*it)->performInit();
  } // next it
}

int
main(int ac, const char** av)
{
  Timer timer;

  IoParams params;
  try
  {
    params.parse(ac,av);
  }
  catch (const std::exception& e)
  {
    std::cerr << " Exception caught while parsing cmd-line: "
    << e.what() << "\n Type --help for command info! " << std::endl;
    exit(1);
  }
  std::cout << " using bounding box: " 
            << (params.useBoundingBox?"yes":"no") << std::endl;

  // load fixed
  MRI* mriFixed = MRIread( getChar( params.strFixed ) );
  if ( !mriFixed )
  {
    std::cerr << " failed reading fixed volumes\n";
    exit(1);
  }

  // load moving
  MRI* mriMoving = MRIread( getChar( params.strMoving ) );
  if ( !mriMoving )
  {
    std::cerr << " failed reading moving volume\n";
    exit(1);
  }

  // load morph
  std::shared_ptr<gmp::VolumeMorph> pmorph(new gmp::VolumeMorph);
  try
  {
    pmorph->load( params.strMorph.c_str(), params.zlibBuffer );
  }
  catch (const char* msg)
  {
    std::cerr << " Failed loading morph: \n"
    << msg << std::endl;
    exit(1);
  }
  pmorph->m_template = mriFixed;
  initOctree2(*pmorph);

  // set the interpolation
  if ( params.strInterp.empty() )
    pmorph->m_interpolationType = SAMPLE_TRILINEAR;
  else if (strcmp(params.strInterp.c_str(), "linear")==0)
    pmorph->m_interpolationType = SAMPLE_TRILINEAR;
  else if (strcmp(params.strInterp.c_str(), "nearest")==0)
    pmorph->m_interpolationType = SAMPLE_NEAREST;
  else pmorph->m_interpolationType = SAMPLE_TRILINEAR;

  if (params.doTest)
  {
    std::cout << " Writing out some tests.\n";
    Timer t1;
    VOL_GEOM vgLike;
    initVolGeom(&vgLike);
    getVolGeom(pmorph->m_template, &vgLike);

    MRI* mriOut  = pmorph->apply_transforms(mriMoving, true, &vgLike);
    std::cout << " morph completed in " << t1.minutes() << " minutes\n";
    //MRIwrite(mriOut, "tmpout1.mgz");
    MRIfree(&mriOut);
  }
  // produce the GCAM
  printf("#VMPC# exportGcam:pre-export VmPeak  %d\n",GetVmPeak());
  GCA_MORPH* gcam = pmorph->exportGcam(mriMoving,
                                       params.useBoundingBox,
                                       params.threshold);
  printf("#VMPC# exportGcam:pre-write VmPeak  %d\n",GetVmPeak());
  GCAMwrite( gcam,
             const_cast<char*>
             (params.strGcam.c_str()));

  // test the presence of the gc structures -- LZ: WHAT DOES THAT DO???
  printf("#VMPC# exportGcam:pre-norm VmPeak  %d\n",GetVmPeak());
  GCAMnormalizeIntensities(gcam, mriFixed);

  if ( params.doTest && (!params.useBoundingBox || 1) )
  {
    std::cout << " applying morph\n"
    << " width = " << gcam->width << std::endl;
    MRI* mriOut = GCAMmorphToAtlas(mriMoving,
                                   gcam,
                                   mriFixed,
                                   0,
                                   pmorph->m_interpolationType);
    std::cout << " AFTER MORPH\n";
    std::cout << " out MRI params = "
    << " width = " << mriOut->width
    << " height = " << mriOut->height
    << " depth = " << mriOut->depth
    << " nframes = " << mriOut->nframes << std::endl;
    //MRIwrite(mriOut, "tmpout2.mgz");
  }
  else
  {
    std::cout << " skipping tmpout2.mgz - using bounding box\n";
  }

  std::cout << " Export performed in " << timer.minutes() << " minutes \n";
  printf("#VMPC# exportGcam VmPeak  %d\n",GetVmPeak());

  return 0;
}

//---------------------------

void
IoParams::parse(int ac, const char** av)
{
  ArgumentParser parser;
  // required
  parser.addArgument("--fixed", 1, String, true);
  parser.addArgument("--moving", 1, String, true);
  parser.addArgument("--morph", 1, String, true);
  parser.addArgument("--out_gcam", 1, String, true);
  // optional
  parser.addArgument("--zlib_buffer", 1, Int);
  parser.addArgument("--bbox_threshold", 1, Int);
  parser.addArgument("--interp", 1, String);
  parser.addArgument("--test");
  // help text
  parser.addHelp(exportGcam_help_xml, exportGcam_help_xml_len);
  parser.parse(ac, av);

  strFixed = parser.retrieve<std::string>("fixed");
  strMoving = parser.retrieve<std::string>("moving");
  strMorph = parser.retrieve<std::string>("morph");
  strGcam = parser.retrieve<std::string>("out_gcam");
  strInterp = parser.retrieve<std::string>("interp");

  doTest = 0;
  if (parser.exists("test")) {
    doTest = 1;
  }

  zlibBuffer = 5;
  if (parser.exists("zlib_buffer")) {
    zlibBuffer = parser.retrieve<int>("zlib_buffer");
  }

  threshold = 0;
  useBoundingBox = false;
  if (parser.exists("bbox_threshold")) {
    threshold = parser.retrieve<int>("bbox_threshold");
    useBoundingBox = true;
  }
}
