
// STL includes
#include <math.h>
#include <fstream>
#include <iostream>

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
#include "argparse.h"
 
#include "mri.h"


#include "applyMorph.help.xml.h"

// required by FreeSurfer
const char *Progname;

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

  std::shared_ptr<gmp::VolumeMorph> pmorph;

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

  void parse(int ac, const char** av);
};

//------------------------------------------------------

int
main(int argc, const char** argv)
{
  // cmd-line
  IoParams params;

  try
  {
    params.parse(argc, argv);
  }
  catch (const char* msg)
  {
    std::cerr << " Exception caught while parsing cmd-line" << std::endl << msg << std::endl;
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
  std::shared_ptr<gmp::VolumeMorph> pmorph(new gmp::VolumeMorph);
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

  typedef std::vector<std::shared_ptr<AbstractFilter> > FilterContainerType;
  FilterContainerType filterContainer;

  for ( std::vector<DataItem>::const_iterator cit = params.items.begin();
        cit != params.items.end(); ++cit )
  {
    std::shared_ptr<AbstractFilter> p;
    switch (cit->m_type)
    {
    case DataItem::surf :
      {
        std::shared_ptr<SurfaceFilter> pTmp(new SurfaceFilter);
        pTmp->strAttached = cit->strAttached;
        pTmp->mriTemplate = mriTemplate;
        p = pTmp;
        break;
      }
    case DataItem::volume :
      {
        p = std::shared_ptr<AbstractFilter>(new VolumeFilter);
      }
      break;
    case DataItem::sprobe :
      {
        std::shared_ptr<SurfaceProbeFilter> pTmp(new SurfaceProbeFilter);
        pTmp->strDestinationSurf = cit->strAttached;
        p = pTmp;
      }
      break;
    case DataItem::snormals :
      {
        std::shared_ptr<SurfaceNormalsProbeFilter> pTmp(new SurfaceNormalsProbeFilter);
        pTmp->strDestinationSurf = cit->strAttached;
        p = pTmp;
      }
      break;
    case DataItem::pointList:
      {
        std::shared_ptr<PointListProbeFilter> pTmp(new PointListProbeFilter);
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

  // apply morph to one point for debug if needed
  if ( !g_vDbgCoords.empty() )
  {
    tDblCoords pt, img;
    for (unsigned int ui=0; ui<3; ++ui)
      pt(ui) = g_vDbgCoords[ui];
    img = pmorph->image(pt);
    std::cout << " computing image for point " << pt << std::endl
    << "\t = " << img << std::endl
    << (img.isValid() ? "valid": "not valid") << std::endl;
  }
  printf("#VMPC# fapplyMorph VmPeak  %d\n",GetVmPeak());
  return 0;
}

//---------------------

void
IoParams::parse(int ac, const char** av)
{
  ArgumentParser parser;
  // required
  parser.addArgument("inputs", '+', String, true);
  parser.addArgument("--template", 1, String, true);
  parser.addArgument("--transform", 1, String, true);
  // optional
  parser.addArgument("--zlib_buffer", 1, Int);
  parser.addArgument("--dbg_coords", 3, Int);
  // help text
  parser.addHelp(applyMorph_help_xml, applyMorph_help_xml_len);
  parser.parse(ac, av);

  strTemplate = parser.retrieve<std::string>("template");
  strTransform = parser.retrieve<std::string>("transform");

  zlibBuffer = 5;
  if (parser.exists("zlib_buffer")) {
    zlibBuffer = parser.retrieve<int>("zlib_buffer");
  }

  if (parser.exists("dbg_coords")) {
    g_vDbgCoords = parser.retrieve<std::vector<int>>("dbg_coords");
  }

  typedef std::vector<std::string> StringVector;
  StringVector container = parser.retrieve<StringVector>("inputs");


  StringVector::const_iterator cit = container.begin();
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
    else if ( *cit == "point_list" ) // This needs to be in data coordinates (and results are in data coordinates as well!)
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

  // int counter = 0;

  //if (tractPointList) 
  this->pmorph->invert();

  int numLines = 0;
  std::string unused;
  while ( std::getline(ifs, unused) )
   ++numLines;
  std::cout << " The number of lines in the input file is " << numLines << std::endl;
  // to rewind the file
  ifs.clear();
  ifs.seekg(0);

  //unsigned int voxInvalid(0);
  // while ( ifs )
  while ( numLines > 0 )
    {
      ifs >> pt(0) >> pt(1) >> pt(2);
      
      img.validate(); // LZ
      img = this->pmorph->image(pt);
      std::cout << " computing image for point " << pt << std::endl
                << "\t = " << img << std::endl
                << (img.isValid() ? "valid" : "not valid") << std::endl;
      /*  if ( !img.isValid() )
	  {
	  if (img.status()==cInvalid)
	  ++voxInvalid;
	  continue;
	  }*/
      
      outputImages.push_back( img );
      numLines --;
      // counter ++;
    }

  ifs.close();
  // std::cout << "In counter :" <<  counter << "\n";

  std::ofstream ofs( this->strOutput.c_str() );
  if ( !ofs ) throw " Failed to open output stream while applying PointListProbeFilter";

  // counter = 0;
  for ( std::vector<Coords3d>::const_iterator cit = outputImages.begin();
        cit != outputImages.end(); ++cit )
    {

      std::cout << " Computed image point " << (*cit)(0) << " " << (*cit)(1) << " " << (*cit)(2) << std::endl;
      if ( cit->isValid() )
        ofs << (*cit)(0) << " " << (*cit)(1) << " " << (*cit)(2) << std::endl;
      else
	ofs << 10000 << " "<< 10000 << " " << 10000 << std::endl; // hack not to lose order
      // counter ++;
    } // next cit
  ofs.close();
  // std::cout << "Out counter :" <<  counter << "\n";
}

