
#include <stdexcept>

#include <itkLinearInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include "morph.h"

#include "tag_fio.h"

template<class T>
T mySqr(T x)
{
  return x*x;
}

namespace gmp
{

AffineTransform3d::AffineTransform3d()
    : Transform<3>(), m_pdata(NULL)
{ }

AffineTransform3d::AffineTransform3d(float* pdata)
    : Transform<3>(), m_pdata(NULL)
{
  /*
  if ( sizeof(pdata) != 12 * sizeof(float) )
    throw "AffineTransform3d constructor - incorrect transform size";
  */
  m_pdata = new float[12];
  std::copy( pdata, pdata+12, m_pdata);
}

AffineTransform3d::tCoords
AffineTransform3d::doOwnImg(const tCoords& pt) const
{
  tCoords ret;

  for (int j=0; j<3; ++j)
    ret(j) = m_pdata[j] * pt(0)
             + m_pdata[j+3] * pt(1)
             + m_pdata[j+6] * pt(2)
             + m_pdata[j+9];

  return ret;
}

void
AffineTransform3d::doInput(std::istream& is)
{
  // read 12 float values
  if (!m_pdata)
    m_pdata = new float[12];

  for (unsigned int ui=0; ui<12; ++ui)
    m_pdata[ui] = TRead<float>(is);
}

void
AffineTransform3d::doOutput(std::ostream& os) const
{
  if (!m_pdata)
    throw "AffineTransform3d save - NULL data";

  for (unsigned int ui=0; ui<12; ++ui)
    TWrite(os, m_pdata[ui]);
}

void
AffineTransform3d::invert()
{
  if ( !m_pdata )
    throw " AffineTransform3d invert - NULL data";

  //std::cout << "AffineTransform3d: invert" << std::endl;

  float* t = m_pdata;

  float fdet = t[0] * ( t[4]*t[8] - t[5]*t[7] )
               - t[3] * ( t[1]*t[8] - t[2]*t[7] )
               + t[6] * ( t[1]*t[5] - t[2]*t[4] );

  if ( std::abs(fdet) < 1.0e-5 )
  {
    std::cerr << " inv_transform -> null det \n";
    exit(1);
  }

  float *inv_t = new float[12];

  inv_t[0] = ( t[4]*t[8] - t[5]*t[7] ) / fdet;
  inv_t[1] = ( t[2]*t[7] - t[1]*t[8] ) / fdet;
  inv_t[2] = ( t[1]*t[5] - t[2]*t[4] ) / fdet;

  inv_t[3] = ( t[5]*t[6] - t[3]*t[8] ) / fdet;
  inv_t[4] = ( t[0]*t[8] - t[2]*t[6] ) / fdet;
  inv_t[5] = ( t[2]*t[3] - t[0]*t[5] ) / fdet;

  inv_t[6] = ( t[3]*t[7] - t[4]*t[6] ) / fdet;
  inv_t[7] = ( t[1]*t[6] - t[0]*t[7] ) / fdet;
  inv_t[8] = ( t[0]*t[4] - t[1]*t[3] ) / fdet;

  inv_t[9] = - t[9]* inv_t[0] - t[10]* inv_t[3] - t[11]* inv_t[6];
  inv_t[10]= - t[9]* inv_t[1] - t[10]* inv_t[4] - t[11]* inv_t[7];
  inv_t[11]= - t[9]* inv_t[2] - t[10]* inv_t[5] - t[11]* inv_t[8];

  delete[] m_pdata;
  m_pdata = inv_t;

}

//-=-----------------------------------------------

DenseDisplacementField::DenseDisplacementField() : m_fieldInterpolator(nullptr), m_maskInterpolator(nullptr) {}

// the field is not kept per se
//
// it is only there through the chosen interpolator
void
DenseDisplacementField::set_field(FieldPointer field)
{
  if ( !m_fieldInterpolator )
    // by default, create a linear interpolator
    m_fieldInterpolator = itk::VectorLinearInterpolateImageFunction<FieldType>::New();
  else
    m_fieldInterpolator = dynamic_cast<FieldInterpolatorType*>
                          (&*m_fieldInterpolator->CreateAnother());

  m_fieldInterpolator->SetInputImage( field );
}

// the mask is not kept per se either
//
// it is only available through the chosen interpolator
void
DenseDisplacementField::set_mask(MaskPointer mask)
{
  if ( !m_maskInterpolator )
    m_maskInterpolator = itk::LinearInterpolateImageFunction<MaskType>::New();
  else
    m_maskInterpolator = dynamic_cast<MaskInterpolatorType*>
                         (&*m_maskInterpolator->CreateAnother());

  m_maskInterpolator->SetInputImage( mask );
}

/*

The format for this file is TAG-based - this is not as good as XML,
but at least allows writing backward-compatible code.

*/
void DenseDisplacementField::doInput(std::istream& is)
{
  using namespace ftags;

  typedef itk::ImageRegionIterator<FieldType> FieldIterator;
  typedef itk::ImageRegionIterator<MaskType>  MaskIterator;
  TagReader tagReader(is);
  FieldType::RegionType region;

  FieldPointer field = FieldType::New();

  while ( tagReader.Read() )
  {
    switch ( tagReader.m_tag )
    {
    case tagSize:
    {
      FieldType::SizeType size;
      size.Fill(0);
      std::istringstream ss(tagReader.m_data);
      for (unsigned int ui=0; ui<3; ++ui)
        size[ui] = TRead<int>(ss);
      region.SetSize(size);
    }
    break;
    case  tagStart:
    {
      FieldType::IndexType start;
      start.Fill(0);
      std::istringstream ss(tagReader.m_data);
      for (unsigned int ui=0; ui<3; ++ui)
        start[ui] = TRead<int>(ss);
      region.SetIndex(start);
    }
    break;
    case tagSpacing:
    {
      double spacing[3];
      std::istringstream ss(tagReader.m_data);
      for (unsigned int ui=0; ui<3; ++ui)
        spacing[ui] = TRead<double>(ss);
      field->SetSpacing(spacing);
    }
    break;
    case tagOrigin:
    {
      double origin[3];
      std::istringstream ss(tagReader.m_data);
      for (unsigned int ui=0; ui<3; ++ui)
        origin[ui] = TRead<float>(ss);
      field->SetOrigin(origin);
    }
    break;
    case tagData:
    {
      // when getting here, allocate the image
      // is there a better timing to do this?
      field->SetRegions( region );
      field->Allocate();

      std::istringstream ss(tagReader.m_data);
      FieldType::PixelType v;

      FieldIterator it(field, field->GetRequestedRegion());
      for ( it.GoToBegin();
            !it.IsAtEnd(); ++it )
      {
        for (unsigned int ui=0;ui<3; ++ui)
          v[ui] = TRead<double>(ss);
        it.Set( v );
      } // next it
    }
    break;
    case tagMask:
    {
      // if getting here (not always the case),
      // set the region and allocate
      MaskPointer mask = MaskType::New();
      mask->SetRegions(region);
      mask->SetSpacing( field->GetSpacing() );
      mask->SetOrigin( field->GetOrigin() );

      std::istringstream ss(tagReader.m_data);

      MaskIterator it(mask, mask->GetRequestedRegion());
      for ( it.GoToBegin();
            !it.IsAtEnd();
            ++it )
      {
        it.Set( TRead<bool>(ss) );
      } // next it
    }
    break;

    }
  }
}

void DenseDisplacementField::doOutput(std::ostream& os) const
{
  if ( !m_fieldInterpolator ) throw " NULL interpolator ";
  FieldConstPointer field = m_fieldInterpolator->GetInputImage();

  using namespace ftags;

  std::string strTag;

  // write size
  strTag = CreateTag(tagSize,
                     this->PrepareTagSize(field)
                    );
  os.write( strTag.c_str(), strTag.size() );

  // write start
  strTag = CreateTag( tagStart,
                      this->PrepareTagStart(field)
                    );
  os.write( strTag.c_str(), strTag.size() );

  // write spacing
  strTag = CreateTag( tagSpacing,
                      this->PrepareTagSpacing(field)
                    );
  os.write( strTag.c_str(), strTag.size() );

  // write origin
  strTag = CreateTag( tagOrigin,
                      this->PrepareTagOrigin(field)
                    );
  os.write( strTag.c_str(), strTag.size() );

  // write data
  strTag = CreateTag( tagData,
                      this->PrepareTagData(field)
                    );
  os.write( strTag.c_str(), strTag.size() );

  // if present, write mask information
  if ( m_maskInterpolator )
  {
    strTag = CreateTag( tagMask,
                        this->PrepareTagMask(m_maskInterpolator->GetInputImage())
                      );
    os.write( strTag.c_str(), strTag.size() );
  }
}

std::string
DenseDisplacementField::PrepareTagSize(FieldConstPointer field) const
{
  std::ostringstream oss;
  for (unsigned int ui=0; ui<3; ++ui)
    TWrite(oss, double(field->GetLargestPossibleRegion().GetSize()[ui] ) );

  return oss.str();
}

std::string
DenseDisplacementField::PrepareTagStart(FieldConstPointer field) const
{
  std::ostringstream oss;
  for (unsigned int ui=0; ui<3; ++ui)
    TWrite(oss, int(field->GetLargestPossibleRegion().GetIndex()[ui]) );

  return oss.str();
}

std::string
DenseDisplacementField::PrepareTagSpacing(FieldConstPointer field) const
{
  std::ostringstream oss;
  for (unsigned int ui=0; ui<3; ++ui)
    TWrite(oss, double(field->GetSpacing()[ui]) );

  return oss.str();
}

std::string
DenseDisplacementField::PrepareTagOrigin(FieldConstPointer field) const
{
  std::ostringstream oss;
  for (unsigned int ui=0; ui<3; ++ui)
    TWrite(oss, double(field->GetOrigin()[ui]) );

  return oss.str();
}

std::string
DenseDisplacementField::PrepareTagData(FieldConstPointer field) const
{
  typedef itk::ImageRegionConstIterator<FieldType> FieldConstIterator;

  std::ostringstream oss;
  FieldConstIterator it(field, field->GetLargestPossibleRegion());

  for (it.GoToBegin();
       !it.IsAtEnd();
       ++it )
    for (unsigned int ui=0;ui<3; ++ui)
      TWrite(oss, it.Get()[ui]);

  return oss.str();
}

// if here, the mask data exists
std::string
DenseDisplacementField::PrepareTagMask(MaskConstPointer mask) const
{
  typedef itk::ImageRegionConstIterator<MaskType> MaskConstIterator;

  std::ostringstream oss;
  MaskConstIterator  it(mask, mask->GetLargestPossibleRegion());

  for (it.GoToBegin();
       !it.IsAtEnd();
       ++it)
    TWrite(oss, it.Get());

  return oss.str();
}

DenseDisplacementField::tCoords
DenseDisplacementField::doOwnImg(const tCoords& pt) const
{
  PointType point;
  for (unsigned int ui=0; ui<3; ++ui)
    point[ui] = pt(ui);

  if ( m_maskInterpolator )
  {
    if ( std::floor(m_maskInterpolator->Evaluate( point )+.5) != 1 )
      return tCoords();
  }

  tCoords retVal;

  OutputType ret = m_fieldInterpolator->Evaluate(  point );

  for (unsigned int ui=0; ui<3; ++ui)
    retVal(ui) = ret[ui];

  return retVal;
}

//--------------------------------------------------

DeltaTransform3d::DeltaTransform3d()
    : m_interpolation(SAMPLE_TRILINEAR), m_field(NULL),
    m_mask(NULL)
{}

DeltaTransform3d::~DeltaTransform3d()
{
  if ( m_field )  MRIfree(&m_field);
  if ( m_mask )   MRIfree(&m_mask);
}

DeltaTransform3d::tCoords
DeltaTransform3d::doOwnImg(const tCoords& pt) const
{

  double val;
  tCoords retVal;

  // sample the mask first - do this on a nearest neighbor method
  if ( m_mask )
  {
    MRIsampleVolumeFrameType( m_mask, pt(0), pt(1), pt(2), 0, m_interpolation, &val );
    if ( !val )
    {
      retVal.invalidate();
      return retVal;
    }
  }


  for (unsigned int ui=0; ui<3; ++ui)
  {
    MRIsampleVolumeFrameType( m_field, pt(0), pt(1), pt(2), ui,
                              m_interpolation, &val);
    retVal(ui) = val;
  }

  retVal += pt;
  return retVal;
}

void
DeltaTransform3d::doInput(std::istream& is)
{
  if ( m_field )
  {
    MRIfree(&m_field);
    m_field = NULL;
  }
  if ( m_mask )
  {
    MRIfree(&m_mask);
    m_mask = NULL;
  }


  int width, height, depth;
  width = TRead<int>(is);
  height = TRead<int>(is);
  depth = TRead<int>(is);

  // alloc volume
  m_field = MRIallocSequence( width,
                              height,
                              depth,
                              MRI_FLOAT,
                              3);

  // populate buffer
  for (int un=0; un<m_field->nframes; ++un)
    for (int uk=0; uk<m_field->depth; ++uk)
      for (int uj=0; uj<m_field->height; ++uj)
        for (int ui=0; ui<m_field->width; ++ui)
          MRIsetVoxVal( m_field, ui,uj,uk, un,
                        TRead<float>(is) );

  bool hasMask = TRead<bool>(is);
  if ( hasMask )
  {
    unsigned int validCount = 0;
    m_mask = MRIalloc( width, height, depth, MRI_UCHAR);
    for (int uk=0; uk<m_mask->depth; ++uk)
      for (int uj=0; uj<m_mask->height; ++uj)
        for (int ui=0; ui<m_mask->width; ++ui)
        {
          MRIvox(m_mask,ui,uj,uk) = TRead<unsigned char>(is);
          if ( MRIvox(m_mask,ui,uj,uk) ) validCount++;
        }
  }
}

void
DeltaTransform3d::doOutput(std::ostream& os) const
{
  if ( !m_field )
    throw "DeltaTransform3d save - NULL data";

  // write size of the buffer
  TWrite(os, m_field->width);
  TWrite(os, m_field->height);
  TWrite(os, m_field->depth);

  // write buffer sequentially, frame by frame
  for (int un=0; un<m_field->nframes; ++un)
    for (int uk=0; uk<m_field->depth; ++uk)
      for (int uj=0; uj<m_field->height; ++uj)
        for (int ui=0; ui<m_field->width; ++ui)
        {
          TWrite(os,
                 MRIgetVoxVal(m_field, ui,uj,uk, un) );
        }

  TWrite(os, (m_mask!=NULL));
  if ( m_mask )
  {
    for (int uk=0; uk<m_mask->depth; ++uk)
      for (int uj=0; uj<m_mask->height; ++uj)
        for (int ui=0; ui<m_mask->width; ++ui)
          TWrite(os, MRIvox(m_mask, ui,uj,uk));

  }

}

void
DeltaTransform3d::invert()
{ 
  double val;
  //  std::cout << "DeltaTransform3d: invert" << std::endl;
  int i, j, k;
  for (i = 0; i < m_field->width; i++)
    for (j = 0; j < m_field->height; j++)
      for (k = 0; k < m_field->depth; k++)
	{
	  MRIsampleVolume(m_field, i, j, k, &val);
	  MRIsetVoxVal(m_field, i, j, k, 0, -val); 
	}
}


//--------------------------------------------------

FemTransform3d::FemTransform3d()
    : Transform<3>()
{
  m_bdbg=false;
  m_signalTopology=false;
}

FemTransform3d::tCoords
FemTransform3d::doOwnImg(const tCoords& pt) const
{
  if (!m_sharedMesh)
    throw std::logic_error("FemTransform3d img -> NULL mesh");

  return m_sharedMesh->dir_img(pt, m_signalTopology);

}

void
FemTransform3d::doInput(std::istream& is)
{
  m_sharedMesh = std::shared_ptr<CMesh3d>(new CMesh3d);

  m_sharedMesh->load(is);
}

void
FemTransform3d::doOutput(std::ostream& os) const
{
  if (!m_sharedMesh)
    throw "FemTransform3d doOutput - NULL mesh";

  m_sharedMesh->save(os);
}

void
FemTransform3d::doOwnInit()
{
  m_sharedMesh->build_index_src();
}

void
FemTransform3d::invert()
{
  if ( m_pInitial )
    throw " FemTransform3d invert - non null initial";

  //std::cout << "FemTransform3d: invert" << std::endl;

  // ugly, but how else?
  std::shared_ptr<CMesh3d>
  pmesh( new CMesh3d(*dynamic_cast<CMesh3d*>(&*m_sharedMesh) ));
  pmesh->invert();
  pmesh->build_index_src();

  m_sharedMesh = pmesh;

}

FemTransform3d::TransformType*
FemTransform3d::convert_to_delta() const
{
  MRI* field = NULL;
  MRI* mask = NULL;
  int width(256), height(256), depth(256);

  // alloc volume
  field = MRIallocSequence( width,
                            height,
                            depth,
                            MRI_FLOAT,
                            3);
  mask = MRIalloc( width,
                   height,
                   depth,
                   MRI_UCHAR);

  tCoords pt,img;

  // populate the buffer
  for (int z=0; z<depth; ++z)
    for (int y=0; y<height; ++y)
      for (int x=0; x<width; ++x)
      {
        pt.validate();
        pt(0) = x;
        pt(1) = y;
        pt(2) = z;
        img = this->img(pt);
        if ( !img.isValid() )
        {
          MRIvox(mask,x,y,z) = 0;
        }
        else
        {
          MRIvox(mask,x,y,z) = 1;
          img -= pt;
          for (unsigned int a=0; a<3; ++a)
            MRIsetVoxVal(field, x,y,z, a, img(a));
        }
      } // next x,y,z

  // initialize the DeltaTransform
  DeltaTransform3d* pdelta = new DeltaTransform3d;
  pdelta->set_field(field);
  pdelta->set_mask(mask);

  return pdelta;
  // no need to delete pdelta, since it will be handled by the smart pointer
}

/****************

****************/


VolumeMorph::VolumeMorph()
{
  m_template = NULL;
  mriCache   = NULL;
  m_interpolationType = SAMPLE_TRILINEAR;

  initVolGeom(&m_vgFixed);
  initVolGeom(&m_vgMoving);
}

VolumeMorph::~VolumeMorph()
{
  if ( mriCache )
    MRIfree(&mriCache);
}

MRI*
VolumeMorph::apply_transforms(MRI* input,
                              bool cacheField,
                              const VG* vgOutput) const
{
  if ( !m_template )
    throw "VolumeMorph apply_transforms - NULL template";
  if ( !input )
    throw "VolumeMorph apply_transforms - NULL input";

  //---------------------------------------
  // allocate output
  VG vg;
  if ( vgOutput )
  {
    vg = *vgOutput;
  }
  else
  {
    // copy geometry from template
    getVolGeom( m_template, &vg );
  }

  MRI* mriOut;
  int nframes = 1;
  if ( input->nframes > 1 ) // for multi-frame volumes
  {
    std::cout << "Multi-frame input \n";
    nframes = input->nframes;
    mriOut = MRIallocSequence( vg.width,
                               vg.height,
                               vg.depth,
                               input->type,
                               input->nframes
                             );
  }
  else
    mriOut = MRIalloc( vg.width,
                       vg.height,
                       vg.depth,
                       input->type
                     );

  // use volume geometry to MRI
  useVolGeomToMRI( &vg, mriOut );

  // init output values
  MRIvalueFill(mriOut, 0.0f);

  // initialize matrices for RAS 2 matrix stuff
  MATRIX* mat_template = NULL;
  MATRIX* mat_subject  = NULL;
  // setup the matrix for the fixed side
  {
    MATRIX* vox2ras_crt = vg_i_to_r(&vg);
    MATRIX* ras2vox_morph = vg_r_to_i(&m_vgFixed);

    if ( !vox2ras_crt || !ras2vox_morph )
      throw " VolumeMorph apply_transforms - NULL matrix ";

    mat_template = MatrixMultiply( ras2vox_morph,
                                   vox2ras_crt,
                                   NULL
                                 );
    MatrixFree(&vox2ras_crt);
    MatrixFree(&ras2vox_morph);
  }

  // setup the matrix for the moving side
  {
    MATRIX* vox2ras_morph = vg_i_to_r(&m_vgMoving);
    MATRIX* ras2vox_crt   = extract_r_to_i(input);

    if ( !vox2ras_morph || !ras2vox_crt )
      throw " VolumeMorph apply_transforms - NULL matrix ";

    mat_subject = MatrixMultiply( ras2vox_crt,
                                  vox2ras_morph,
                                  NULL
                                );
    MatrixFree(&vox2ras_morph);
    MatrixFree(&ras2vox_crt);
  }

  //MatrixPrint( stdout, mat_template );
  //MatrixPrint( stdout, mat_subject );

  if ( cacheField )
  {
    if ( mriCache )
      MRIfree(&mriCache);
    mriCache = MRIallocSequence( m_template->width,
                                 m_template->height,
                                 m_template->depth,
                                 MRI_FLOAT, 4 ); // 4 frames - one for each direction + 1 to indicate a valid voxel
  }

  try
  {
    unsigned int voxInvalid(0), voxValid(0);
    tCoords pt, img;
    TransformContainerType::const_iterator cit;
    double val;
    float valvect[nframes];
    
    if ( cacheField )
      for (int z=0; z<mriOut->depth; ++z)
        for (int y=0; y<mriOut->height; ++y)
          for (int x=0; x<mriOut->width; ++x)
            MRIsetVoxVal(mriCache, x,y,z, 3, 0);

    VECTOR *vFixed, *vMoving, *vTmp;
    vFixed  = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(vFixed, 4) = 1.0;
    vMoving = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(vMoving,4) = 1.0;
    vTmp    = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(vTmp,4) = 1.0;

    for (int z=0; z<mriOut->depth; ++z)
    {
      if ( !(z%10) ) std::cout << " z = " << z << std::endl;
      for (int y=0; y<mriOut->height; ++y)
        for (int x=0; x<mriOut->width; ++x)
        {
          //std::cout << " in the FOR loop " << std::endl;
          //-------------------------
          // do RAS conversion
          VECTOR_ELT(vTmp, 1) = x;
          VECTOR_ELT(vTmp, 2) = y;
          VECTOR_ELT(vTmp, 3) = z;

          vFixed = MatrixMultiply( mat_template, vTmp, vFixed );
          //std::cout << " vFixed computed " << std::endl;

          pt.validate();
          pt(0) = V3_X( vFixed );
          pt(1) = V3_Y( vFixed );
          pt(2) = V3_Z( vFixed );
          //-------------------------

          img = this->image(pt);
          //std::cout << " img computed " << std::endl;

          if ( !img.isValid() )
          {
            if (img.status()==cInvalid)
              ++voxInvalid;
            continue;
          }

          //--------------------------
          // convert RAS on the
          //     moving side
          //

          V3_X(vTmp) = img(0);
          V3_Y(vTmp) = img(1);
          V3_Z(vTmp) = img(2);

          vMoving = MatrixMultiply( mat_subject, vTmp, vMoving );
          //std::cout << " vMoving computed " << std::endl;

          img(0) = V3_X( vMoving);
          img(1) = V3_Y( vMoving);
          img(2) = V3_Z( vMoving);
          //std::cout << " img assigned " << img(0) <<"," << img(1) <<"," <<img(2) << std::endl;

          //--------------------------
          // do nothing if out of bounds
          if ( img(0)<0 || img(0)>input->width-1 ||
               img(1)<0 || img(1)>input->height-1 ||
               img(2)<0 || img(2)>input->depth-1 ) continue;

          ++voxValid;
	  if (nframes == 1)
	    {
	      MRIsampleVolumeType( input, img(0), img(1), img(2), &val, m_interpolationType);
	      MRIsetVoxVal( mriOut, x,y,z, 0, val);
	    }
	  else
	    {
	      MRIsampleSeqVolumeType(  input, img(0), img(1), img(2), valvect, 0, nframes-1, m_interpolationType);
	      for (int ii = 0; ii < nframes; ii++)
		MRIsetVoxVal( mriOut, x,y,z,ii, valvect[ii]);
	    }

          if ( cacheField )
          {
            tCoords bufPt(img);
            bufPt -= pt;
            for ( unsigned int dir=0; dir<3; ++dir)
              MRIsetVoxVal( mriCache, x,y,z, dir, bufPt(dir) );
            MRIsetVoxVal( mriCache, x,y,z, 3, 1);
          }
        } // next x,y,z
    }
    std::cout << " Invalid voxels = " << voxInvalid << std::endl
    << " Valid = " << voxValid << std::endl;
  }
  catch (const gmpErr& excp)
  {
    std::cerr << " Exception caught -> " << excp.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << " Unhandled exception!!!\n";
    exit(1);
  }

  return mriOut;
}

// TODO: eliminate template reliance
MRI*
VolumeMorph::convert_transforms() const
//MRI* input,
//    bool cacheField,
//         const VG* vgOutput) const
{

  bool  cacheField = true;
  if ( !m_template )
    throw "VolumeMorph convert_transforms - NULL template";
  //if ( !input )
  //throw "VolumeMorph convert_transforms - NULL input";

  //---------------------------------------
  // allocate output
  VG vg;
  // copy geometry from template
  getVolGeom( m_template, &vg );

  //MRI* mriOut = MRIalloc( vg.width,
  //     vg.height,
  //     vg.depth,
  //     input->type
  //     );

  // use volume geometry to MRI -- can I do it with a sequence??
  //useVolGeomToMRI( &vg, mriOut );

  // init output values
  //MRIvalueFill(mriOut, 0.0f);

  // initialize matrices for RAS 2 matrix stuff
  MATRIX* mat_template = NULL;
  MATRIX* mat_subject  = NULL;
  // setup the matrix for the fixed side
  {
    //    MATRIX* vox2ras_crt = vg_i_to_r(&vg);
    MATRIX* vox2ras_crt = vg_i_to_r(&m_vgFixed);
    MATRIX* ras2vox_morph = vg_r_to_i(&m_vgFixed);

    if ( !vox2ras_crt || !ras2vox_morph )
      throw " VolumeMorph convert_transforms (1) - NULL matrix ";

    mat_template = MatrixMultiply( ras2vox_morph,
                                   vox2ras_crt,
                                   NULL
                                 );
    MatrixFree(&vox2ras_crt);
    MatrixFree(&ras2vox_morph);
  }

  // setup the matrix for the moving side
  {
    MATRIX* vox2ras_morph = vg_i_to_r(&m_vgMoving);
    //TODO: is this just the inverse of the above?
    MATRIX* ras2vox_crt = MatrixInverse(vox2ras_morph, NULL);

    if ( !vox2ras_morph || !ras2vox_crt )
      throw " VolumeMorph convert_transforms (2) - NULL matrix ";

    mat_subject = MatrixMultiply( ras2vox_crt,
                                  vox2ras_morph,
                                  NULL
                                );
    MatrixFree(&vox2ras_morph);
    MatrixFree(&ras2vox_crt);
  }


  if (cacheField)
  {
    if ( mriCache )
      MRIfree(&mriCache);
    mriCache = MRIallocSequence( m_vgFixed.width,  //m_template->width,
                                 m_vgFixed.height, //m_template->height,
                                 m_vgFixed.depth,  //m_template->depth,
                                 MRI_FLOAT, 4 ); // 4 frames - one for each direction + 1 to indicate a valid voxel
  }

  MRI* mriWarpAsVolume = MRIallocSequence( m_vgFixed.width,  //m_template->width,
                         m_vgFixed.height, //m_template->height,
                         m_vgFixed.depth,  //m_template->depth,
                         MRI_FLOAT, 3 );

  try
  {
    unsigned int voxInvalid(0), voxValid(0);
    tCoords pt, img;
    TransformContainerType::const_iterator cit;

    if ( cacheField )
      for (int z=0; z<m_vgFixed.depth; ++z)
        for (int y=0; y<m_vgFixed.height; ++y)
          for (int x=0; x<m_vgFixed.width; ++x)
            MRIsetVoxVal(mriCache, x,y,z, 3, 0);

    VECTOR *vFixed, *vMoving, *vTmp;
    vFixed  = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(vFixed, 4) = 1.0;
    vMoving = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(vMoving,4) = 1.0;
    vTmp    = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(vTmp,4) = 1.0;

    for (int z=0; z<m_vgFixed.depth; ++z)
    {
      if ( !(z%10) ) std::cout << " z = " << z << std::endl;
      for (int y=0; y<m_vgFixed.height; ++y)
        for (int x=0; x<m_vgFixed.width; ++x)
        {
          //-------------------------
          // do RAS conversion
          VECTOR_ELT(vTmp, 1) = x;
          VECTOR_ELT(vTmp, 2) = y;
          VECTOR_ELT(vTmp, 3) = z;

          vFixed = MatrixMultiply( mat_template, vTmp, vFixed );

          pt.validate();
          pt(0) = V3_X( vFixed );
          pt(1) = V3_Y( vFixed );
          pt(2) = V3_Z( vFixed );
          //-------------------------

          img = this->image(pt);
          if ( !img.isValid() )
          {
            if (img.status()==cInvalid)
              ++voxInvalid;
            continue;
          }

          //--------------------------
          // convert RAS on the
          //     moving side
          //

          V3_X(vTmp) = img(0);
          V3_Y(vTmp) = img(1);
          V3_Z(vTmp) = img(2);

          vMoving = MatrixMultiply( mat_subject, vTmp, vMoving );

          img(0) = V3_X( vMoving);
          img(1) = V3_Y( vMoving);
          img(2) = V3_Z( vMoving);

          //--------------------------
          // do nothing if out of bounds
          /*if ( img(0)<0 || img(0)>input->width-1 ||
          img(1)<0 || img(1)>input->height-1 ||
          img(2)<0 || img(2)>input->depth-1 ) continue;*/

          if ( img(0)<0 || img(0)>m_vgMoving.width-1 ||
               img(1)<0 || img(1)>m_vgMoving.height-1 ||
               img(2)<0 || img(2)>m_vgMoving.depth-1 ) continue;

          ++voxValid;
          //MRIsampleVolumeType( input, img(0), img(1), img(2), &val, m_interpolationType);
          //MRIsetVoxVal( mriOut, x,y,z,0, val);

          tCoords bufPt(img);
          bufPt -= pt;
          if ( cacheField )
          {
            for ( unsigned int dir=0; dir<3; ++dir)
              MRIsetVoxVal( mriCache, x,y,z, dir, bufPt(dir) );
            MRIsetVoxVal( mriCache, x,y,z, 3, 1);
          }
          for ( unsigned int dir=0; dir<3; ++dir)
            MRIsetVoxVal( mriWarpAsVolume, x,y,z, dir, bufPt(dir) );
        } // next x,y,z
    }
    std::cout << " Invalid voxels = " << voxInvalid << std::endl
    << " Valid = " << voxValid << std::endl;
  }
  catch (const gmpErr& excp)
  {
    std::cerr << " Exception caught -> " << excp.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << " Unhandled exception!!!\n";
    exit(1);
  }

  return mriCache;
  //return mriWarpAsVolume;
}


MRIS*
VolumeMorph::apply_transforms(MRIS* input) const
{
  MRI_SURFACE* const mris = MRISclone( input );

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  int const nvertices = input->nvertices;
  for (int ui = 0; ui < nvertices; ++ui) 
  {
    VERTEX* pvtxIn = &input->vertices[ui];

    tDblCoords pt;
    pt.validate();
    pt(0) = pvtxIn->x;
    pt(1) = pvtxIn->y;
    pt(2) = pvtxIn->z;

    pt = image( pt );
    if ( !pt.isValid() ) continue; // better leave it as it was if it's not working

    MRISsetXYZ(mris, ui, pt(0), pt(1), pt(2));
  }

  return mris;
}

tDblCoords
VolumeMorph::image(const tCoords& _pt) const
{
  TransformContainerType::const_iterator cit;
  tCoords pt(_pt), ret;
  bool bDone(false);

  for ( cit = m_transforms.begin();
        cit != m_transforms.end() && !bDone;
        ++cit)
  {
    ret = (*cit)->img(pt);
    pt = ret;
    if (!pt.isValid()) bDone = true;
  }
  return ret;

}

#if 0
void
VolumeMorph::save(const char* fname)
{
  std::ostringstream os(std::ios::binary);

  // write no. of transforms
  TWrite(os, (unsigned int)m_transforms.size() );

  // write transforms
  for ( TransformContainerType::iterator it = m_transforms.begin();
        it != m_transforms.end(); ++it )
    saveTransform( os, *it);

  // compress result
  ZlibStringCompressor compressor;
  std::string compressed = compressor.compress( os.str(),
                           Z_BEST_COMPRESSION);
  std::cout << " compressor buffer size = " << compressor.getBufferSize() << std::endl;

  // write compressed buffer to file
  std::ofstream ofs(fname, std::ios::binary);
  if ( !ofs )
    throw "VolumeMorph save - failed to open output stream";

  ofs.write(compressed.c_str(),
            sizeof(char) * compressed.size() );

}

void
VolumeMorph::save(const char* fname)
{
  // serialize transforms to help handle the compressor problem
  this->serialize();

  std::ostringstream os(std::ios::binary);

  // write no. of transforms
  TWrite(os, (unsigned int)m_transforms.size() );

  // write transforms
  for ( TransformContainerType::iterator it = m_transforms.begin();
        it != m_transforms.end(); ++it )
  {
    std::ostringstream osit(std::ios::binary);
    saveTransform( osit, *it);

    ZlibStringCompressor compressor;
    std::cout << " transform ostream size = " << osit.str().size() << std::endl;
    std::string compressed = compressor.compress( osit.str(),
                             Z_BEST_COMPRESSION);
    unsigned long bufferSize = sizeof(char) * compressed.size() ;

    TWrite( os, bufferSize );
    std::cout << " transform size = " << (int)bufferSize << std::endl;
    os.write( compressed.c_str(),
              bufferSize );
  } // next transform it

  // write compressed buffer to file
  std::ofstream ofs(fname, std::ios::binary);
  if ( !ofs )
    throw "VolumeMorph save - failed to open output stream";

  ofs.write(os.str().c_str(),
            sizeof(char) * os.str().size() );

}
#else

void
VolumeMorph::save(const char* fname)
{
  std::string strTag;

  // serialize transforms to help handle the compressor problem
  this->serialize();

  std::ostringstream os(std::ios::binary);

  // vol geom fixed
  strTag = ftags::CreateTag( tagVgFixed,
                             this->PrepareTagVolGeom(m_vgFixed)
                           );
  os.write( strTag.c_str(), strTag.size() );

  // vol geom moving
  strTag = ftags::CreateTag( tagVgMoving,
                             this->PrepareTagVolGeom(m_vgMoving)
                           );
  os.write( strTag.c_str(), strTag.size() );

  // write transforms
  ZlibStringCompressor compressor; // may have a mem leak, so move outside of loop
  for ( TransformContainerType::iterator it = m_transforms.begin();
        it != m_transforms.end(); ++it )
  {
    std::ostringstream osit(std::ios::binary);
    saveTransform(osit, *it);

    std::string strBuf = compressor.compress( osit.str(),
                         Z_BEST_COMPRESSION );
    std::cout << " writing transform size = " << strBuf.size() << std::endl;

    strTag =  ftags::CreateTag( tagTransform,
                                strBuf
                              );
    os.write( strTag.c_str(),
              strTag.size()
            );
  } // next transform

  std::cout << " writing morph to file " << fname << std::endl;
  // write compressed buffer to file
  std::ofstream ofs(fname, std::ios::binary);
  if ( !ofs )
    throw "VolumeMorph save - failed to open output stream";

  ofs.write(os.str().c_str(),
            sizeof(char) * os.str().size() );

}

#endif

void
VolumeMorph::load(const char* fname,
                  unsigned int bufferMultiplier,
                  bool clearExisting)
{
  // read the transform
  // for backwards compatibility, read old if extension

  // 1. get the extension
  std::string strFname(fname);
  std::string::size_type pos = strFname.find_last_of(".");
  if ( pos == std::string::npos ) throw " Couldn't find extension in filename. ";

  std::string strExtension = strFname.substr(pos+1);

  std::cout << " extension = " << strExtension << std::endl;
  if ( strExtension == "nm3d" )
    this->load_old(fname, bufferMultiplier, clearExisting);
  else
    this->load_new(fname, bufferMultiplier, clearExisting);

}

#if 0
void
VolumeMorph::load_old(const char* fname,
                      unsigned int bufferMultiplier,
                      bool clearExisting)
{
  std::cout << " VolumeMorph::load_old\n";
  std::ifstream ifs(fname, std::ios::binary);
  if ( !ifs )
    throw "VolumeMorph load - failed to open input stream";

  // determine the size of the file
  ifs.seekg(0, std::ios::end);
  unsigned int size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);

  char* dataBuffer = new char[size];
  ifs.read(dataBuffer, size);

  const std::string strCompressed(dataBuffer, size);

  ZlibStringCompressor compressor;
  compressor.m_bufferAllocationMultiplier = bufferMultiplier;
  const std::string strInflated = compressor.inflate( strCompressed );

  std::istringstream is(strInflated, std::ios::binary);

  // read number of transforms
  unsigned int noTransforms = TRead<unsigned int>(is);

  if ( noTransforms>0 && clearExisting ) m_transforms.clear();

  for (unsigned int ui=0;
       ui < noTransforms;
       ++ui )
  {
    TransformPointer t = loadTransform(is);
    m_transforms.push_back(t);
  }
}

#else
/*

modified algorithm compared to past.
1. read number of transforms
2. for each transform - read size of compressed buffer
3. read transform buffer
4. uncompress it
5. populate transform

*/
void
VolumeMorph::load_old(const char* fname,
                      unsigned int bufferMultiplier,
                      bool clearExisting)
{
  std::ifstream ifs(fname, std::ios::binary);
  if ( !ifs ) throw "VolumeMorph load - failed to open input stream";

  unsigned int noTransforms = TRead<unsigned int>(ifs);
  if ( noTransforms>0 && clearExisting ) m_transforms.clear();

  std::cout << " got here\n";
  for (unsigned int ui=0;
       ui < noTransforms;
       ++ui)
  {
    // read buffer size
    unsigned long bufferSize = TRead<unsigned long>(ifs);
    std::cout << " buffer size = " << (int)bufferSize << std::endl;
    // allocate transform buffer
    char* dataBuffer = new char[bufferSize];
    ifs.read(dataBuffer, bufferSize);
    const std::string strCompressed(dataBuffer,bufferSize);

    // init compressor
    ZlibStringCompressor compressor;
    compressor.m_bufferAllocationMultiplier = bufferMultiplier;
    const std::string strInflated = compressor.inflate(strCompressed);

    // read the transform
    std::istringstream is(strInflated, std::ios::binary);
    TransformPointer t = loadTransform(is);
    m_transforms.push_back(t);
  } // next ui

}
#endif

void
VolumeMorph::load_new(const char* fname,
                      unsigned int bufferMultiplier,
                      bool clearExisting)
{
  ZlibStringCompressor compressor; // memory leak in compressor?
  std::ifstream ifs(fname, std::ios::binary);
  if ( !ifs ) throw "VolumeMorph load - failed to open input stream";

  ftags::TagReader tagReader(ifs);

  if ( clearExisting ) m_transforms.clear();

  while ( tagReader.Read() )
  {
    switch ( tagReader.m_tag )
    {
    case tagVgFixed:
      this->ReadTagVolGeom
      (
        std::string(tagReader.m_data, tagReader.m_len),
        m_vgFixed
      );

      break;
    case tagVgMoving:
      this->ReadTagVolGeom
      (
        std::string(tagReader.m_data, tagReader.m_len),
        m_vgMoving
      );

      break;
    case tagTransform:
    {
      //ZlibStringCompressor compressor;// might be a memleak, so made function scope
      compressor.m_bufferAllocationMultiplier = bufferMultiplier;
      std::cout << " data size = " << tagReader.m_len << std::endl;

      const std::string strCompressed(tagReader.m_data, tagReader.m_len);
      const std::string strInflated = compressor.inflate(strCompressed);
      std::istringstream is(strInflated);
      TransformPointer t = loadTransform(is);
      m_transforms.push_back(t);
    }
    break;
    default:
      ;
    }
  } // tagReader
}

std::string
VolumeMorph::PrepareTagVolGeom(const VOL_GEOM& vg)
{
  std::ostringstream oss;
  TWrite(oss, vg.valid);

  TWrite(oss, int(vg.width) );
  TWrite(oss, int(vg.height) );
  TWrite(oss, int(vg.depth) );

  TWrite(oss, float(vg.xsize) );
  TWrite(oss, float(vg.ysize) );
  TWrite(oss, float(vg.zsize) );

  TWrite(oss, float(vg.x_r) );
  TWrite(oss, float(vg.x_a) );
  TWrite(oss, float(vg.x_s) );

  TWrite(oss, float(vg.y_r) );
  TWrite(oss, float(vg.y_a) );
  TWrite(oss, float(vg.y_s) );

  TWrite(oss, float(vg.z_r) );
  TWrite(oss, float(vg.z_a) );
  TWrite(oss, float(vg.z_s) );

  TWrite(oss, float(vg.c_r) );
  TWrite(oss, float(vg.c_a) );
  TWrite(oss, float(vg.c_s) );

  return oss.str();
}

void
VolumeMorph::ReadTagVolGeom(const std::string& strData,
                            VOL_GEOM& vg)
{
  std::istringstream iss(strData);

  vg.valid = TRead<int>(iss);

  vg.width = TRead<int>(iss);
  vg.height = TRead<int>(iss);
  vg.depth = TRead<int>(iss);

  vg.xsize = TRead<float>(iss);
  vg.ysize = TRead<float>(iss);
  vg.zsize = TRead<float>(iss);

  vg.x_r = TRead<float>(iss);
  vg.x_a = TRead<float>(iss);
  vg.x_s = TRead<float>(iss);

  vg.y_r = TRead<float>(iss);
  vg.y_a = TRead<float>(iss);
  vg.y_s = TRead<float>(iss);

  vg.z_r = TRead<float>(iss);
  vg.z_a = TRead<float>(iss);
  vg.z_s = TRead<float>(iss);

  vg.c_r = TRead<float>(iss);
  vg.c_a = TRead<float>(iss);
  vg.c_s = TRead<float>(iss);
}

GCA_MORPH*
VolumeMorph::exportGcam(MRI* mriMoving,
                        bool useTemplateBoundingBox,
                        int thresh,
                        int padding) const
{
  double dval, dmean(0), dnum(0);
  
  std::cout << "exportGcam: padding is = " << padding <<   "\n";

  if ( !m_template )
  {
    throw "VolumeMorph exportGcam -> NULL template";
  }
  if ( !mriCache )
  {
    std::cerr << "VolumeMorph::exportGcam -> no mriCache - this will take some time\n";
    VOL_GEOM vgLike;
    initVolGeom(&vgLike);
    getVolGeom(this->m_template, &vgLike);
    MRI* mriTmp = this->apply_transforms(mriMoving,true, &vgLike); // start by caching the data
    if ( !mriTmp ) return NULL;
    MRIfree(&mriTmp);
  }
  else
    std::cout << "exportGcam: mriCache exists\n";

  //allocate the Gcam
  MRI_REGION box;
  if ( useTemplateBoundingBox )
  {
    std::cout << "Using template bounding box\n";
    MRIboundingBox( m_template, thresh, &box);
    if ( box.x - padding < 0 ||
         box.y - padding < 0 ||
         box.z - padding < 0 ||
         box.x + box.dx + padding >= m_template->width ||
         box.y + box.dy + padding >= m_template->height ||
         box.z + box.dz + padding >= m_template->depth )
      padding = 0;
    box.dx += 2* padding;
    box.dy += 2* padding;
    box.dz += 2* padding;

    box.x -= padding;
    box.y -= padding;
    box.z -= padding;
    std::cout << " padding = " << padding << std::endl;
  }
  else
  {
    std::cout << "Not using bounding box\n";
    box.x = 0;
    box.y = 0;
    box.z = 0;
    box.dx = m_template->width;
    box.dy = m_template->height;
    box.dz = m_template->depth;
  }
  std::cout << " box = " << std::endl
  << box.x << " , " << box.y << " , " << box.z << std::endl
  << box.dx << " , " << box.dy << " , " << box.dz << std::endl;

  GCA_MORPH* gcam = GCAMalloc( box.dx,
                               box.dy,
                               box.dz );
  gcam->type = GCAM_VOX;
  gcam->spacing = 1;

  // assign atlas geometry
  //
  // the atlas must have the same size as the GCAM
  // so crop it first
  // to keep it consistent with the RAS,
  // update the center of the volume and keep the same i_to_r
  MRI* croppedTemplate = MRIextractRegion( this->m_template,
                         NULL,
                         &box);
  if (mriMoving)  GCAMinitVolGeom(gcam, mriMoving, croppedTemplate);
  MRIfree(&croppedTemplate);
  gcam->ninputs = 1;

  // populate the morph
  //  const unsigned int depth( (unsigned int)gcam->depth ),
  // height( (unsigned int)gcam->height ),
  //width( (unsigned int)gcam->width );

  unsigned int gcamInvalidVoxels = 0;

  for (int z(0), zbox(box.z); z<box.dz; ++z, ++zbox)
    for (int y(0), ybox(box.y); y<box.dy; ++y, ++ybox)
      for (int x(0), xbox(box.x); x<box.dx; ++x, ++xbox)
      {
        GMN* pnode = &gcam->nodes[x][y][z];

        if ( MRIgetVoxVal(mriCache, xbox,ybox,zbox, 3) < 1 )
        {
          pnode->invalid = GCAM_POSITION_INVALID;
          pnode->x = pnode->origx = 0;
          pnode->y = pnode->origy = 0;
          pnode->z = pnode->origz = 0;
          pnode->label = 0;
          ++gcamInvalidVoxels;
          continue;
        }

        pnode->invalid = GCAM_VALID;
        pnode->x = pnode->origx = MRIgetVoxVal( mriCache, xbox,ybox,zbox, 0) +x;
        pnode->y = pnode->origy = MRIgetVoxVal( mriCache, xbox,ybox,zbox, 1) +y;
        pnode->z = pnode->origz = MRIgetVoxVal( mriCache, xbox,ybox,zbox, 2) +z;

        pnode->orig_area = gcam->spacing*gcam->spacing*gcam->spacing;

        // allocate GC1D values
        pnode->gc = alloc_gcs(1, GCA_NO_MRF, 1); // 0 labels, 0 flags, 1 input

        dval = MRIgetVoxVal( m_template, xbox,ybox,zbox,0);

        pnode->gc->means[0] = dval;
        if ( std::abs(dval) > 1.0e-10 )
        {
          pnode->label = 128; // something real here - taken from GCAMcreateIntensityImage
          dmean += dval;
          dnum += 1.0;
        }
      }

  if ( !dmean )
  {
    std::cout << " NO VALID LABEL\n";
    return gcam;
  }

  dmean /= dnum;
  dmean = mySqr(.05 * dmean);
  for (int z(0); z<box.dz; ++z)
    for (int y(0); y<box.dy; ++y)
      for (int x(0); x<box.dx; ++x)
      {
        GMN* pnode = &gcam->nodes[x][y][z];
        if ( !pnode || pnode->invalid == GCAM_POSITION_INVALID ) continue;
        pnode->gc->covars[0] = dmean;
      }

  std::cout << " gcam export invalid voxels = " << gcamInvalidVoxels << std::endl;
  return gcam;
}

void
VolumeMorph::invert()
{
  TransformContainerType tmpContainer;

  //std::cout << "VolumeMorph: invert" << std::endl;

  //int counter = 0;

  for ( TransformContainerType::iterator it = m_transforms.begin();
        it != m_transforms.end(); ++it )
  {
    if ( (*it)->initial() )
      throw " VolumeMorph invert - transform has initial";

    (*it)->invert();
    tmpContainer.push_front( *it );

    //counter ++;
  } // next it
  m_transforms = tmpContainer;

  //std::cout << "in morph:invert ==> counter = " <<  counter-1 << std::endl;
}

void
VolumeMorph::serialize()
{
  TransformContainerType::iterator it = m_transforms.begin();
  TransformPointer pt;

  while ( it != m_transforms.end() )
  {
    if ( (pt = (*it)->initial()) )
    {
      (*it)->setInitial( TransformPointer() );
      it = m_transforms.insert(it, pt);
    }
    else
    {
      ++it;
    }
  }
}


std::shared_ptr<Transform<3> >
loadTransform(std::istream& is, unsigned int zlibBufferMultiplier)
{
  // read the string preceding the data
  unsigned int uilen = TRead<unsigned int>(is);
  std::cout << " uilen = " << uilen << std::endl;

  char *buffer = new char[uilen];
  is.read(buffer, uilen);

  std::string strDescription(buffer, uilen);
  delete[] buffer;

  std::string strName;
  std::string::size_type sepPos = strDescription.find("|");
  if ( sepPos != std::string::npos )
  {
    strName = strDescription.substr(sepPos+1);
    strDescription = strDescription.substr(0, sepPos);
  }

  std::shared_ptr<Transform<3> > bp;

  if ( strDescription == "affine" )
    bp = std::shared_ptr<Transform<3> >(new AffineTransform3d);
  else if ( strDescription == "fem" )
    bp = std::shared_ptr<Transform<3> >(new FemTransform3d);
  else if ( strDescription == "id" )
    bp = std::shared_ptr<Transform<3> >(new IdentityTransform3d);
  else if ( strDescription == "field" )
    bp = std::shared_ptr<Transform<3> >(new DeltaTransform3d);
  else
    throw "loadTransform - unknown transform type";

  std::cout << " loading transform id string = " << strName << std::endl;
  bp->load(is, zlibBufferMultiplier);

  return bp;

}

void
saveTransform(std::ostream& os,
              std::shared_ptr<Transform<3> > ptransform)
{
  std::cout << " saveTransform code\n";


  std::string strType;

  if ( dynamic_cast<AffineTransform3d*>( &*ptransform) )
    strType = "affine";
  else if ( dynamic_cast<FemTransform3d*>( &*ptransform) )
    strType = "fem";
  else if ( dynamic_cast<IdentityTransform3d*>( &*ptransform) )
    strType = "id";
  else if (dynamic_cast<DeltaTransform3d*>( &*ptransform))
    strType = "field";
  else
    throw "saveTransform - unsupported transform type";

  strType += "|";
  strType += ptransform->m_strName;
  TWrite(os, (unsigned int)strType.size() );
  os.write(strType.c_str(), strType.size()*sizeof(char));

  ptransform->save(os);
}


}

