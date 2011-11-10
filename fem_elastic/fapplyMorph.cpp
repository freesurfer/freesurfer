
// STL includes
#include <cmath>
#include <fstream>
#include <iostream>

// BOOST
#include <boost/program_options.hpp>
#include <boost/progress.hpp>

// MPI
//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR
//#include <mpi.h>

// OWN
#include "morph.h"

#include "surf_utils.h"
#include "morph_utils.h"

// FreeSurfer
extern "C"
{
#include "mri.h"
};

// required by FreeSurfer
char *Progname;

int tractPointList = 0;

////////////

template<class T>
T sqr(const T& val)
{
  return val*val;
}

std::vector<int> g_vDbgCoords;

//void initOctree( gmp::VolumeMorph& morph);

struct DataItem
{
  enum Type
  {
    surf,
    volume,
    sprobe,
    snormals,
    pointList
  };

  Type m_type;
  std::string strInput;
  std::string strOutput;
  std::string strAttached;

  std::string interpolation;
};


class AbstractFilter
{
public:
  AbstractFilter()
  {}
  std::string strInput;
  std::string strOutput;

  boost::shared_ptr<gmp::VolumeMorph> pmorph;

  virtual void Execute() =0;
  virtual ~AbstractFilter()
  {}
};


class VolumeFilter : public AbstractFilter
{
public:
  void Execute();
  std::string strGcam;
  int m_interpolationType;
};

class SurfaceFilter : public AbstractFilter
{
public:
  std::string strAttached;
  MRI* mriTemplate;
  void Execute();
};

class SurfaceProbeFilter : public AbstractFilter
{
public:
  SurfaceProbeFilter() : m_mrisTemplate(NULL), m_mrisSubject(NULL)
  {}
  std::string strDestinationSurf;
  virtual void Execute();
protected:
  void LoadInputData();

  MRIS*  m_mrisTemplate;
  MRIS*  m_mrisSubject;
  void ConvertSurfaceToVoxelSpace(const VG& vg, MRIS* mris);
};

class SurfaceNormalsProbeFilter : public SurfaceProbeFilter
{
public:
  SurfaceNormalsProbeFilter() : SurfaceProbeFilter()
  {}
  void Execute();
};

class PointListProbeFilter : public AbstractFilter
{
public:
  void Execute();

  std::string strMode;

};

//------------------------------------------------------

struct IoParams
{
  std::vector<DataItem> items;

  std::string strTemplate;
  std::string strTransform;
  std::string strGcam; // option to export gcam -- not yet implemented

  unsigned int zlibBuffer;

  void parse(int ac, char* av[]);
};

//------------------------------------------------------

int
main(int argc,
     char* argv[])
{
  boost::timer t0;

  // cmd-line

  IoParams params;

  try
  {
    params.parse(argc, argv);
  }
  catch (const char* msg)
  {
    std::cerr << " Exception caught while parsing cmd-line\n"
    << msg << std::endl;
    exit(1);
  }

  // load template
  std::cout<<"Template name:";  // TODO: why does not it read nii volumes????
  std::cout<< const_cast<char*>( params.strTemplate.c_str()) << "\n";
  MRI* mriTemplate = MRIread( const_cast<char*>( params.strTemplate.c_str()) );
  if ( !mriTemplate )
  {
    std::cerr << " Failed reading template volume "
    << params.strTemplate << std::endl;
    exit(1);
  }
  std::cout<<"After loading template\n";

  // load template
  std::cout<<"Template name:";  // TODO: why does not it read nii volumes????
  std::cout<< const_cast<char*>( params.strTemplate.c_str()) << "\n";
  mriTemplate = MRIread( const_cast<char*>( params.strTemplate.c_str()) );
  if ( !mriTemplate )
    {
      std::cerr << " Failed reading template volume "
		<< params.strTemplate << std::endl;
      exit(1);
    }
  std::cout<<"After loading template\n";
  

  // load transform
  boost::shared_ptr<gmp::VolumeMorph> pmorph(new gmp::VolumeMorph);
  pmorph->m_template = mriTemplate;

  try
  {
    pmorph->load( params.strTransform.c_str(), params.zlibBuffer );
  }
  catch (const char* msg)
    {
    std::cerr << " Exception caught while loading transform\n"
	      << msg << std::endl;
    exit(1);
    }
  std::cout << " loaded transform\n";
  initOctree(*pmorph);

  typedef std::vector<boost::shared_ptr<AbstractFilter> > FilterContainerType;
  FilterContainerType filterContainer;

  for ( std::vector<DataItem>::const_iterator cit = params.items.begin();
        cit != params.items.end(); ++cit )
  {
    boost::shared_ptr<AbstractFilter> p;
    switch (cit->m_type)
    {
    case DataItem::surf :
      {
        boost::shared_ptr<SurfaceFilter> pTmp(new SurfaceFilter);
        pTmp->strAttached = cit->strAttached;
        pTmp->mriTemplate = mriTemplate;
        p = pTmp;
        break;
      }
    case DataItem::volume :
      {
        p = boost::shared_ptr<AbstractFilter>(new VolumeFilter);
      }
      break;
    case DataItem::sprobe :
      {
        boost::shared_ptr<SurfaceProbeFilter> pTmp(new SurfaceProbeFilter);
        pTmp->strDestinationSurf = cit->strAttached;
        p = pTmp;
      }
      break;
    case DataItem::snormals :
      {
        boost::shared_ptr<SurfaceNormalsProbeFilter> pTmp(new SurfaceNormalsProbeFilter);
        pTmp->strDestinationSurf = cit->strAttached;
        p = pTmp;
      }
      break;
    case DataItem::pointList:
      {
        boost::shared_ptr<PointListProbeFilter> pTmp(new PointListProbeFilter);
        pTmp->strMode = cit->strAttached;
        p = pTmp;
      }
    default:
      ;
    }

    if ( cit->interpolation=="linear" )
      pmorph->m_interpolationType = SAMPLE_TRILINEAR;
    else if ( cit->interpolation=="nearest")
      pmorph->m_interpolationType = SAMPLE_NEAREST;

    p->pmorph = pmorph;
    p->strInput = cit->strInput;
    p->strOutput = cit->strOutput;

    filterContainer.push_back(p);

  } // next cit

  try
  {
    // apply each filter
    for ( FilterContainerType::iterator it = filterContainer.begin();
          it != filterContainer.end(); ++it )
    {
      std::cout << " executing filter on file "
      << (*it)->strInput << std::endl;
      (*it)->Execute();
    } // next it
  }
  catch (const char* msg)
  {
    std::cerr << " Exception caught while applying filters " 
              << msg << std::endl;
    exit(1);
  }

  std::cout << " morph applied in " << t0.elapsed() / 60. << " minutes \n";

  // apply morph to one point for debug if needed
  if ( !g_vDbgCoords.empty() )
  {
    tDblCoords pt, img;
    for (unsigned int ui=0; ui<3; ++ui)
      pt(ui) = g_vDbgCoords[ui];
    img = pmorph->image(pt);
    std::cout << " computing image for point " << pt << std::endl
    << "\t = " << img << std::endl
    << (img.isValid()?"not valid":"") << std::endl;
  }
  return 0;
}

//---------------------

void
IoParams::parse(int ac,
                char* av[])
{
  zlibBuffer = 5;

  namespace po = boost::program_options;
  typedef std::vector<std::string> StringContainerType;
  StringContainerType container;

  po::options_description desc("Allowed Options");

  desc.add_options()
  ("help", " produce help message")
  ("template", po::value<std::string>(), " template volume ")
  ("sourcevol", po::value<std::string>(), " source (moving) volume ")
  ("transform", po::value<std::string>(), " transform file")
  //("gcam", po::value(&strGcam), " if present, will write a gcam at that location" )
  ("zlib_buffer", po::value(&zlibBuffer), " zlib buffer pre-allocation multiplier")
  ("dbg_coords", po::value(&g_vDbgCoords)->multitoken(), " debug coordinates")
  ;

  po::options_description hidden;
  hidden.add_options()
  ("data", po::value<StringContainerType>(), " input files");

  po::positional_options_description p;
  p.add("data", -1);

  po::options_description cmd_line;
  cmd_line.add(desc).add(hidden);

  po::variables_map vm;
  po::store( po::command_line_parser(ac,av).
             options(cmd_line).positional(p).run(), vm);
  po::notify(vm);

  if ( vm.count("help") )
  {
    std::cout << desc << std::endl;
    exit(0);
  }

  if ( !vm.count("template") )
    throw " IoParams - you need to specify a template volume";
  strTemplate = vm["template"].as<std::string>();

  if ( !vm.count("transform") )
    throw " IoParams - you need to specify a transform";
  strTransform = vm["transform"].as<std::string>();

  // process the data vector
  if ( !vm.count("data") )
    throw " IoParams - data missing";
  container = vm["data"].as<StringContainerType>();

  StringContainerType::const_iterator cit = container.begin();

  while ( cit != container.end() )
  {
    // read a data item
    if ( *cit == "surf" )
    {
      if ( ++cit == container.end() )
        throw " Incomplete data item";

      DataItem item;
      item.m_type = DataItem::surf;
      item.strInput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item";
      item.strOutput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item";
      item.strAttached = *cit;

      items.push_back(item);
    }
    else if ( *cit == "sprobe" ||
              *cit == "snormals" )
    {
      if ( ++cit == container.end() )
        throw " Incomplete data item";

      DataItem item;
      if ( *cit == "sprobe" )
        item.m_type = DataItem::sprobe;
      else
        item.m_type = DataItem::snormals;
      item.strInput = *cit;

      if ( ++cit == container.end() ) throw "Incomplete data item";
      item.strAttached = *cit;

      if ( ++cit == container.end() ) throw "Incomplete data item";
      item.strOutput = *cit;

      items.push_back(item);
    }
    else if ( *cit == "vol" )
    {
      if ( ++cit == container.end() )
        throw " Incomplete data item";

      DataItem item;
      item.m_type = DataItem::volume;
      item.strInput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item";
      item.strOutput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item";
      item.interpolation = *cit;

      items.push_back(item);
    }
    else if ( *cit == "point_list" )
    {

      if ( ++cit == container.end() )
        throw " Incomplete data item ";

      DataItem item;
      item.m_type = DataItem::pointList;
      item.strInput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item ";
      item.strOutput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item ";
      item.strAttached = *cit;

      items.push_back( item );
    }
    else if ( *cit == "tract_point_list" )
    {
      tractPointList = 1;
      if ( ++cit == container.end() )
        throw " Incomplete data item ";

      DataItem item;
      item.m_type = DataItem::pointList;
      item.strInput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item ";
      item.strOutput = *cit;

      if ( ++cit == container.end() )
        throw " Incomplete data item ";
      item.strAttached = *cit;

      items.push_back( item );
    }
    else
      throw " Unrecognized data item type";

    ++cit;
  }
}

//----------------------------------------------


void
VolumeFilter::Execute()
{
  // load input volume
  MRI* mri = MRIread( const_cast<char*>
                      ( strInput.c_str() ) );
  if ( !mri )
  {
    std::cerr << " Failed reading input volume "
    << strInput << std::endl;
    exit(1);
  }

  VOL_GEOM vgLike;
  initVolGeom(&vgLike);
  getVolGeom(pmorph->m_template, &vgLike);

  MRI* mriOut  = pmorph->apply_transforms(mri,
                                          false,
                                          &vgLike);

  std::cout << " done morphing - will write next\n";
  MRIwrite(mriOut,
           const_cast<char*>
           (strOutput.c_str()) );
  std::cout << " done\n";

  // free data
  MRIfree(&mriOut);
  MRIfree(&mri);
}

void
SurfaceFilter::Execute()
{
  //std::cout << " before serialize\n";
  pmorph->serialize();
  //std::cout << " before inverting\n";
  //pmorph->invert();

  // load the surface
  MRIS* mris = MRISread
               ( const_cast<char*>( strInput.c_str() ) );
  if ( !mris )
    throw " SurfaceFilter Execute - failed reading surface file";

  // load attached volume
  MRI* mriAttached = MRIread
                     ( const_cast<char*>( strAttached.c_str() ) );
  if ( !mriAttached )
    throw " SurfaceFilter Execute - failed reading attached volume";

  // change from RAS to VOX
  convert_surf_to_vox( mris, mriAttached );

  // apply morph
  MRIS* mrisOut = pmorph->apply_transforms( mris );

  // change from VOX to RAS
  MRIScopyVolGeomFromMRI( mrisOut, mriTemplate );
  convert_vox_to_surf( mrisOut, mriTemplate);

  // save surface
  MRISwrite( mrisOut,
             const_cast<char*>( strOutput.c_str()) );

  // free data
  MRISfree(&mris);
  MRISfree(&mrisOut);
  MRIfree(&mriAttached);
}

void
SurfaceProbeFilter::Execute()
{

  this->LoadInputData();

  // now apply the morph for every vertex and
  //     write the difference
  //
  // if the image of a vertex is not valid, mark it as -.1

  // iterate through vertices
  //
  // assumption - the correspondence between surface vertices
  // is index/position based

  VERTEX* pvtxTemplate = &(this->m_mrisTemplate->vertices[0]);
  VERTEX* pvtxSubject = &(this->m_mrisSubject->vertices[0]);
  const unsigned int nVertices = m_mrisTemplate->nvertices;

  //double ddist;
  Coords3d img;

  for (unsigned int ui = 0;
       ui < nVertices;
       ++ui, ++pvtxTemplate, ++pvtxSubject )
  {
    // apply the morph to the point
    // remember the vertex is now in index-space
    Coords3d pt;
    pt(0) = pvtxTemplate->x;
    pt(1) = pvtxTemplate->y;
    pt(2) = pvtxTemplate->z;

    img = pmorph->image( pt );

    if ( !img.isValid() )
      pvtxTemplate->curv = -.1;
    else
    {
      pvtxTemplate->curv = std::sqrt
                           ( sqr(img(0)-pvtxSubject->x) +
                             sqr(img(1)-pvtxSubject->y) +
                             sqr(img(2)-pvtxSubject->z)
                           );
    }
  } // next ui, pvtxTemplate, pvtxSubject

  // finally, save the curvature file
  MRISwriteCurvature( this->m_mrisTemplate,
                      const_cast<char*>( strOutput.c_str() )
                    );
  MRISfree( &this->m_mrisTemplate );
  MRISfree( &this->m_mrisSubject );
}

void
SurfaceProbeFilter::LoadInputData()
{
  // before loading anything, check that the geometries are valid
  if ( !this->pmorph->vgFixed().valid ||
       !this->pmorph->vgMoving().valid )
  {
    std::cerr << " Morph doesn't have valid geometries \n";
    return;
  }

  // load template surface
  this->m_mrisTemplate = MRISread
                         ( const_cast<char*>( this->strInput.c_str() ) );
  if ( !this->m_mrisTemplate )
    throw "SurfaceProbeFilter Execute - failed to read input surface ";

  // load subject surface
  m_mrisSubject = MRISread
                  ( const_cast<char*>( this->strDestinationSurf.c_str() ) );
  if ( !this->m_mrisSubject )
    throw "SurfaceProbeFilter Execute - failed to read input surface ";

  // test that surfaces have the same number of vertices
  if ( this->m_mrisTemplate->nvertices != this->m_mrisSubject->nvertices )
    throw "SurfaceProbeFilter Execute - surfaces do not have equal number of vertices ";

  // convert surfaces to voxel space
  this->ConvertSurfaceToVoxelSpace( this->pmorph->vgFixed(),
                                    this->m_mrisTemplate );
  this->ConvertSurfaceToVoxelSpace( this->pmorph->vgMoving(),
                                    this->m_mrisSubject );

}

void
SurfaceProbeFilter::ConvertSurfaceToVoxelSpace(const VG& vg,
    MRIS* mris)
{
  MRI* mri = MRIalloc( vg.width,
                       vg.height,
                       vg.depth,
                       MRI_UCHAR );
  useVolGeomToMRI( &vg, mri );

  convert_surf_to_vox( mris, mri );

  MRIfree(&mri);
}

void
SurfaceNormalsProbeFilter::Execute()
{
  this->LoadInputData();

  // compute normals for the template surface
  MRIScomputeNormals( this->m_mrisSubject );

  VERTEX* pvtxTemplate = &(this->m_mrisTemplate->vertices[0]);
  VERTEX* pvtxSubject = &(this->m_mrisSubject->vertices[0]);
  const unsigned nVertices = m_mrisTemplate->nvertices;

  //double ddist;
  Coords3d img, normal;

  for (unsigned int ui=0;
       ui < nVertices;
       ++ui, ++pvtxTemplate, ++pvtxSubject )
  {
    // apply the morph to the point first
    Coords3d pt;
    pt(0) = pvtxTemplate->x;
    pt(1) = pvtxTemplate->y;
    pt(2) = pvtxTemplate->z;

    img = this->pmorph->image( pt );

    if ( !img.isValid() )
      pvtxTemplate->curv = -.1;
    else
    {
      pt(0) = pvtxSubject->x;
      pt(1) = pvtxSubject->y;
      pt(2) = pvtxSubject->z;

      normal(0) = pvtxSubject->nx;
      normal(1) = pvtxSubject->ny;
      normal(2) = pvtxSubject->nz;

      pvtxTemplate->curv = std::abs( dot( normal, img-pt ) );
    }
  } // next ui, pvtxTemplate, pvtxSubject

  MRISwriteCurvature( this->m_mrisTemplate,
                      const_cast<char*>( this->strOutput.c_str() )
                    );

  MRISfree( &this->m_mrisTemplate );
  MRISfree( &this->m_mrisSubject );
}

void
PointListProbeFilter::Execute()
{
  std::ifstream ifs( this->strInput.c_str() );
  if ( !ifs ) throw " Failed to open input file while applying PointListProbeFilter ";

  std::vector<Coords3d> outputImages;
  Coords3d pt, img; 

  //int counter = 0;

  if (tractPointList) 
    this->pmorph->invert();
  
  //unsigned int voxInvalid(0);
  while ( ifs )
    {
      ifs >> pt(0) >> pt(1) >> pt(2);
      
      img = this->pmorph->image(pt);
      /*  if ( !img.isValid() )
	  {
	  if (img.status()==cInvalid)
	  ++voxInvalid;
	  continue;
	  }*/
      
      outputImages.push_back( img );
      //counter ++;
    }

  ifs.close();
  //std::cout << "In counter :" <<  counter << "\n";

  std::ofstream ofs( this->strOutput.c_str() );
  if ( !ofs ) throw " Failed to open output stream while applying PointListProbeFilter";

  //counter = 0;
  for ( std::vector<Coords3d>::const_iterator cit = outputImages.begin();
        cit != outputImages.end(); ++cit )
    {
      if ( cit->isValid() )
        ofs << (*cit)(0) << " " << (*cit)(1) << " " << (*cit)(2) << std::endl;
      else
	ofs << 10000 << " "<< 10000 << " " << 10000 << std::endl; // hack not to lose order
      //  counter ++;
    } // next cit
  ofs.close();
  //  std::cout << "Out counter :" <<  counter << "\n";
  }

