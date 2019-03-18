/*

Gheorghe Postelnicu, 2006

This binary creates a morph, given transforms

*/

// STL includes
#include <iterator>
#include <iostream>
#include <string>

// OWN
#include "morph.h"
#include "transformUtils.h"

// FS
#include "argparse.h"
 
#include "gcamorph.h"
#include "mri.h"


#include "createMorph.help.xml.h"

const char *Progname;

//==================================================

struct MorphItem
{
  std::string type;
  std::string file;
};

std::vector<int> g_vDbgCoords;

struct IoParams
{
  std::string strOut;
  std::string strTemplate;
  std::string strSubject;
  typedef std::vector<MorphItem> MorphContainerType;
  MorphContainerType items;

  typedef std::vector<int> IntVectorType;
  IntVectorType vDbgCoords;

  void parse(int ac, const char** av);
};

//==================================================


// assume field was already allocated
MRI*
CopyGcamToDeltaField(GCA_MORPH* gcam, MRI* field);

int
main(int argc, const char** argv)
{
  IoParams params;
  try
  {
    params.parse(argc, argv);
  }
  catch (const std::exception& excp)
  {
    std::cerr << " Exception while parsing cmd-line\n"
    << excp.what() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cerr << " unhandled exception caught while parsing cmd line\n";
    exit(1);
  }


  std::shared_ptr<gmp::VolumeMorph> pmorph(new gmp::VolumeMorph);

  if ( !params.strTemplate.empty() )
  {
    MRI* mri = MRIread( const_cast<char*>(params.strTemplate.c_str()) );
    if (!mri)
      std::cerr << " Failed to open template mri "
                << params.strTemplate << std::endl;
    else
      pmorph->set_volGeom_fixed(mri);
    MRIfree(&mri);
  }
  if ( !params.strSubject.empty() )
  {
    MRI* mri = MRIread( const_cast<char*>(params.strSubject.c_str()) );
    if (!mri)
      std::cerr << " Failed to open subject mri "
                << params.strSubject << std::endl;
    else
      pmorph->set_volGeom_moving(mri);
    MRIfree(&mri);
  }

  IoParams::MorphContainerType::iterator it;
  for ( it = params.items.begin();
        it != params.items.end();
        ++it )
  {
    if ( it->type == "affine" )
    {
      std::shared_ptr<gmp::AffineTransform3d>
      paffine( new gmp::AffineTransform3d);
      float* pf = read_transform( it->file.c_str() );
      if ( !pf )
      {
        std::cerr << " failed to read transform\n";
        exit(1);
      }
      paffine->set_pars( pf);
      pmorph->m_transforms.push_back( paffine );

    }
    else if ( it->type == "volume" )
    {
      std::shared_ptr<gmp::DeltaTransform3d>
      pvol( new gmp::DeltaTransform3d);
      MRI* field = MRIread
                   ( const_cast<char*>( it->file.c_str() ) );
      if ( !field )
      {
        std::cerr << " failed to read field "  << std::endl;
        exit(1);
      }
      pvol->set_field(field);

      pmorph->m_transforms.push_back( pvol );
    }
    else if ( it->type == "gcam" )
    {
      std::shared_ptr<gmp::DeltaTransform3d>
      pvol( new gmp::DeltaTransform3d);
      GCA_MORPH* gcam = GCAMread( const_cast<char*>( it->file.c_str() ) );
      if ( !gcam )
      {
        std::cerr << " failed to read GCAM " << std::endl;
        exit(1);
      }
      // create bogus MRI to hold the geometry
      MRI* mri_template = MRIread( const_cast<char*>
                                   ( params.strTemplate.c_str() ) );
      MRI* mriBuf = MRIalloc( gcam->image.width,
                              gcam->image.height,
                              gcam->image.depth,
                              MRI_UCHAR);
      useVolGeomToMRI( &gcam->image, mriBuf );
      GCAMrasToVox(gcam, mriBuf );
      MRIfree( &mriBuf );

      std::cout << " atlas geometry = "
      << gcam->atlas.width << " , " << gcam->atlas.height
      << " , " << gcam->atlas.depth << std::endl
      << " image geometry = "
      << gcam->image.width << " , " << gcam->image.height
      << " , " << gcam->image.depth << std::endl;


      MRI* mri = MRIallocSequence( mri_template->width,
                                   mri_template->height,
                                   mri_template->depth,
                                   MRI_FLOAT,
                                   3);
      MRIcopyHeader(mri_template, mri);
      g_vDbgCoords = params.vDbgCoords;
      try
      {
        //MRI* mask =
        CopyGcamToDeltaField(gcam, mri);

        pvol->set_field(mri);
        //pvol->set_mask(mask);
        pmorph->m_transforms.push_back( pvol );
      }
      catch (const std::string& e)
      {
        std::cerr << " Exception caught while processing GCAM node\n"
        << e << std::endl;
        exit(1);
      }
      MRIfree(&mri_template);
    }
    else if ( it->type == "morph" )
    {
      std::shared_ptr<gmp::VolumeMorph> tmpMorph(new gmp::VolumeMorph);
      try
      {
        tmpMorph->load( it->file.c_str() );
      }
      catch (const char* msg)
      {
        std::cerr << " Exception caught while loading morph in file "
        << it->file << std::endl;
        exit(1);
      }
      for ( gmp::VolumeMorph::TransformContainerType::iterator transformIter
            = tmpMorph->m_transforms.begin();
            transformIter != tmpMorph->m_transforms.end();
            ++transformIter )
        pmorph->m_transforms.push_back( *transformIter );

    }
    else if ( it->type == "mesh" )
    {
      std::shared_ptr<gmp::FemTransform3d> pfem(new gmp::FemTransform3d);
      std::shared_ptr<CMesh3d> pmesh(new CMesh3d);
      pmesh->load( it->file.c_str() );
      pfem->m_sharedMesh = pmesh;
      pmorph->m_transforms.push_back(pfem);
    }
    else
    {
      std::cerr << " unhandled transform type "
      << it->type << std::endl;
    }
  } // next it

  // finally write morph file
  try
  {
    pmorph->save( params.strOut.c_str() );
  }
  catch (const char* msg)
  {
    std::cerr << " Exception caught while saving morph\n"
    << msg << std::endl;
    exit(1);
  }

  printf("#VMPC# fcreateMorph VmPeak  %d\n",GetVmPeak());
  return 0;
}


void
IoParams::parse(int ac, const char** av)
{
  ArgumentParser parser;
  // required
  parser.addArgument("--in", '+', String, true);
  parser.addArgument("--out", 1, String, true);
  // optional
  parser.addArgument("--template", 1, String);
  parser.addArgument("--subject", 1, String);
  parser.addArgument("--dbg", 3, Int);
  // help text
  parser.addHelp(createMorph_help_xml, createMorph_help_xml_len);
  parser.parse(ac, av);

  typedef std::vector<std::string> StringVector;
  StringVector vItems = parser.retrieve<StringVector>("in");
  strOut = parser.retrieve<std::string>("out");

  strTemplate = parser.retrieve<std::string>("template");
  strSubject = parser.retrieve<std::string>("subject");

  if (parser.exists("dbg")) {
    vDbgCoords = parser.retrieve<std::vector<int>>("dbg");
  }

  if (vItems.empty()) {
    std::cout << " nothing to do - no transforms specified\n";
    exit(0);
  }

  if (vItems.size() % 2) {
    std::cerr << " odd number of transform tokens\n";
    exit(1);
  }

  for ( StringVector::iterator it = vItems.begin(); it != vItems.end(); ) {
    MorphItem item;
    item.type = *it++;
    item.file = *it++;
    this->items.push_back(item);
  }
}

inline int myNint(double x)
{
  return static_cast<int>(std::floor( x +.5));
}

// returns the mask if necessary
MRI*
CopyGcamToDeltaField(GCA_MORPH* gcam, MRI* field)
{
  if ( !field ) throw std::string(" Field not set in CopyGcamToDeltaField\n");

  // allocate the mask
  MRI* mask = MRIalloc( field->width,
                        field->height,
                        field->depth,
                        MRI_UCHAR );
  unsigned int maskCounter = 0;
  unsigned char ucOne(1);//, ucZero(0);
  // fill the mask with true values to start with

  for (int z=0; z<field->depth; ++z)
    for (int y=0; y<field->height; ++y)
      for (int x=0; x<field->width; ++x)
        MRIvox(mask, x,y,z) = ucOne;


  GMN* pnode = NULL;

  MATRIX* v2r_atlas = VGgetVoxelToRasXform( &gcam->atlas, NULL, 0);
  MATRIX* r2v_image = extract_r_to_i( field );
  MATRIX* transform = MatrixMultiply( r2v_image, v2r_atlas, NULL);

  std::cout << " vox 2 ras atlas = \n";
  MatrixPrint( stdout, v2r_atlas );
  std::cout << " ras 2 vox image = \n";
  MatrixPrint( stdout, r2v_image );
  std::cout << " product = \n";
  MatrixPrint( stdout, transform);
  std::cout << std::endl;

  // using the above matrix, go through the nodes of the GCAM
  // and associate the delta to the transformed node
  //
  // !!!! The strong underlying assumption is that the IMAGE and ATLAS associated
  // with the GCAM share the same RAS
  //
  VECTOR *vx = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(vx, 4) = 1.0;
  VECTOR_ELT(vx, 1) = 0.0;
  VECTOR_ELT(vx, 2) = 0.0;
  VECTOR_ELT(vx, 3) = 0.0;

  VECTOR *vfx= VectorAlloc(4, MATRIX_REAL);

  MatrixMultiply(transform, vx, vfx);
  std::cout << " transform at origin\n";
  MatrixPrint(stdout, vfx);

  int shift_x, shift_y, shift_z;

  shift_x = myNint( VECTOR_ELT(vfx, 1) );
  shift_y = myNint( VECTOR_ELT(vfx, 2) );
  shift_z = myNint( VECTOR_ELT(vfx, 3) );

  unsigned int invalidCount = 0;
  int img_x, img_y, img_z;
  for ( int z=0, maxZ=gcam->depth; z<maxZ; ++z)
    for ( int y=0, maxY=gcam->height; y<maxY; ++y)
      for ( int x=0, maxX=gcam->width; x<maxX; ++x)
      {
#if 0
        // indexing is 1-based - NR
        VECTOR_ELT(vx, 1) = x;
        VECTOR_ELT(vx, 2) = y;
        VECTOR_ELT(vx, 3) = z;

        MatrixMultiply(transform, vx, vfx);

        img_x = myNint( VECTOR_ELT(vfx, 1) );
        img_y = myNint( VECTOR_ELT(vfx, 2) );
        img_z = myNint( VECTOR_ELT(vfx, 3) );
#endif

        img_x = x + shift_x;
        img_y = y + shift_y;
        img_z = z + shift_z;

        if ( img_x < 0 || img_x > field->width-1 ||
             img_y < 0 || img_y > field->height-1 ||
             img_z < 0 || img_z > field->depth-1 )
        {
          maskCounter++;
          continue;
        }
        pnode = &gcam->nodes[x][y][z];

        //if (pnode->invalid) continue;
        if ( pnode->invalid == GCAM_POSITION_INVALID )
        {
          ++invalidCount;
          continue;
        }

#if 0
        if ( img_x == g_vDbgCoords[0] &&
             img_y == g_vDbgCoords[1] &&
             img_z == g_vDbgCoords[2] )
        {
          std::cout << " node " << img_x
          << " , " << img_y
          << " , " << img_z
          << " comes from "
          << x << " , " << y << " , " << z
          << " -> " << pnode->x
          << " , " << pnode->y
          << " , " << pnode->z << std::endl;
        }
#endif

        MRIsetVoxVal(mask, img_x, img_y, img_z, 0, ucOne );

        MRIsetVoxVal( field, img_x, img_y, img_z, 0,
                      pnode->x - img_x );
        MRIsetVoxVal( field, img_x, img_y, img_z, 1,
                      pnode->y - img_y );
        MRIsetVoxVal( field, img_x, img_y, img_z, 2,
                      pnode->z - img_z );

      }
  std::cout << " invalid voxel count = " << invalidCount << std::endl;
  if ( !g_vDbgCoords.empty() )
  {
    int x,y,z;

    pnode = &gcam->nodes[gcam->width/2][gcam->height/2][gcam->depth/2];

    for (unsigned int ui=0; ui<3; ++ui)
      VECTOR_ELT(vx, ui+1) = g_vDbgCoords[ui];

    MATRIX* minv = MatrixInverse(transform, NULL);
    MatrixMultiply(minv, vx, vfx);

    std::cout << " debugging at coords = \n";
    std::copy( g_vDbgCoords.begin(), g_vDbgCoords.end(),
               std::ostream_iterator<int>( std::cout, " ") );
    std::cout << std::endl;

    x = myNint( VECTOR_ELT( vfx, 1) );
    y = myNint( VECTOR_ELT( vfx, 2) );
    z = myNint( VECTOR_ELT( vfx, 3) );

    std::cout << " linear transf to get gcam xyz = "
    << x << " , " << y << " , " << z << std::endl;
    if ( x < 0 || x > gcam->width -1 ||
         y < 0 || y > gcam->height-1 ||
         z < 0 || z > gcam->depth -1 )
      std::cout << " out of bounds\n";
    else
    {
      pnode = &gcam->nodes[x][y][z];
      std::cout << " value of gcam = "
      << pnode->x << " , "
      << pnode->y << " , "
      << pnode->z << std::endl;
    }
  }


//   for( int z(0), maxZ(field->depth); z<maxZ; ++z)
//     for( int y(0), maxY(field->height); y<maxY; ++y)
//       for( int x(0), maxX(field->width); x<maxX; ++x)
//  {

//    if ( !GCAMsampleMorph(gcam, x,y,z,
//     &pxd, &pyd, &pzd) )
//      {
//        MRIvox(mask, x,y,z) = ucZero;
//        maskCounter++;
//        continue;
//      }
//    MRIsetVoxVal( field, x,y,z, 0,
//    pxd - x );
//    MRIsetVoxVal( field, x,y,z, 1,
//    pyd - y );
//    MRIsetVoxVal( field, x,y,z, 2,
//    pzd - z );

//  } // next x,y,z

  if ( !maskCounter )
  {
    MRIfree(&mask);
    return NULL;
  }
  return mask;
}
