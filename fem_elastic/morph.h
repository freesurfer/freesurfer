/**

Gheorghe Postelnicu,
June 2006

A volume morph is obtained through a transform. 
The role of a transform is to provide a mapping 
of a location of origin to a destination point. 
This can be done in a number of ways.  Will create
an abstract transform class and then transform
composition can be implemented in a transparent
polymorphic way.

**/

#ifndef _h_morph_h_
#define _h_morph_h_

// STL  includes
#include <fstream>
#include <iostream>
#include <list>

// ITK
#include <itkImage.h>
#include <itkVectorInterpolateImageFunction.h>
#include <itkInterpolateImageFunction.h>
#include <itkVector.h>

// OWN includes
#include "coords.h"
#include "mesh.h"
#include "gmpError.h"

#include "fem_3d.h"



#include "mri.h"
#define class xclass
#include "gcamorph.h"
#undef class
;


/*

The initial transform is to be applied to the left (arrow notation).

A full chain of transformations can thus be specified

*/

namespace gmp
{

template<int n> class Transform;

std::shared_ptr<Transform<3> > loadTransform(std::istream& is, unsigned int zlibBufferMultiplier=5);
void saveTransform(std::ostream& os, std::shared_ptr<Transform<3> > ptransform);

template <int n>
class Transform
{
public:
  typedef TCoords<double,n> tCoords;

  Transform() : m_strName(), m_pInitial()
  { }

  virtual ~Transform()
  {}

  void setInitial(std::shared_ptr<Transform> transform)
  {
    m_pInitial = transform;
  }
  std::shared_ptr<Transform> initial() const
  {
    return m_pInitial;
  }

  tCoords img(const tCoords& pt) const
  {
    tCoords ptBuf = this->processInitial(pt);
    if ( !ptBuf.isValid() ) return ptBuf;

    return this->doOwnImg(ptBuf);
  }
  virtual void invert() = 0;

  void performInit()
  {
    if ( m_pInitial ) m_pInitial->performInit();
    this->doOwnInit();
  }

  void load(std::istream& is, unsigned int zlibBufferMultiplier=5)
  {
    if (TRead<bool>(is))
      m_pInitial = loadTransform(is,zlibBufferMultiplier);
    doInput(is);
  }

  void save(std::ostream& os) const
  {
    TWrite(os, (m_pInitial!=NULL) );
    if ( m_pInitial )
      saveTransform(os, m_pInitial);
    doOutput(os);
  }

  void print(std::string indent=std::string()) const
  {
    std::cout << indent << " name = " << m_strName << std::endl;
    if (m_pInitial)
    {
      std::cout << "   Initial transform \n";
      m_pInitial->print("\t");
    }

  }
  std::string m_strName; // useful for debugging and other tracking purposes

protected:
  virtual tCoords processInitial(const tCoords& pt) const
  {
    if (m_pInitial) return m_pInitial->img(pt);
    return pt;
  }
  std::shared_ptr<Transform<n> >  m_pInitial;

  virtual void doInput(std::istream& is)=0;
  virtual void doOutput(std::ostream& os) const=0;
  virtual tCoords doOwnImg(const tCoords& pt) const=0;
  virtual void doOwnInit()
{};
};

class IdentityTransform3d : public Transform<3>
{
public:
  typedef gmp::Transform<3> Superclass;
  typedef Superclass::tCoords tCoords;

  void invert()
  {}
protected:
  void doInput(std::istream& is)
  {}
  void doOutput(std::ostream& os) const
    {}

  virtual tCoords doOwnImg(const tCoords& pt) const
  {
    return pt;
  }
};

class AffineTransform3d : public gmp::Transform<3>
{
public:
  typedef gmp::Transform<3> Superclass;
  typedef Superclass::tCoords tCoords;

  AffineTransform3d();
  AffineTransform3d(float* pdata);
  void set_pars(float* pdata)
  {
    if (m_pdata) delete[] m_pdata;
    m_pdata = pdata;
  }
  float* pars() const
  {
    return m_pdata;
  }

  void invert();
protected:
  void doInput(std::istream& is);
  void doOutput(std::ostream& os) const;

  virtual tCoords doOwnImg(const tCoords& pt) const;

private:
  float* m_pdata;
};

class DeltaTransform3d : public gmp::Transform<3>
{
public:
  typedef gmp::Transform<3> Superclass;
  typedef Superclass::tCoords tCoords;

  DeltaTransform3d();
  ~DeltaTransform3d();

  int m_interpolation;

  void invert();

  // assumes ownership
  void set_field(MRI* field)
  {
    if (m_field) MRIfree(&m_field);
    m_field = field;
  }
  void set_mask(MRI* mask)
  {
    if (m_mask) MRIfree(&m_mask);
    m_mask = mask;
  }

  const MRI* field() const
  {
    return m_field;
  }

protected:
  void doInput(std::istream& is);
  void doOutput(std::ostream& os) const;

  virtual tCoords doOwnImg(const tCoords& pt) const;

private:
  MRI* m_field;
  MRI* m_mask;
};

class DenseDisplacementField : public gmp::Transform<3>
{
public:
  typedef gmp::Transform<3> Superclass;
  typedef Superclass::tCoords tCoords;
  typedef itk::Vector<float, 3> FieldPixelType;
  typedef itk::Image<FieldPixelType,3> FieldType;
  typedef FieldType::Pointer  FieldPointer;
  typedef FieldType::ConstPointer FieldConstPointer;
  typedef itk::Image<bool,3> MaskType;
  typedef MaskType::Pointer  MaskPointer;
  typedef MaskType::ConstPointer MaskConstPointer;

  DenseDisplacementField();

  void set_field(FieldPointer field);
  void set_mask(MaskPointer mask);

protected:

  enum IoTags
  {
    tagSize = 1,
    tagStart,
    tagSpacing,
    tagOrigin,
    tagData,
    tagMask
  };

  void doInput(std::istream& is);
  void doOutput(std::ostream& os) const;

  virtual tCoords doOwnImg(const tCoords& pt) const;

private:
  typedef itk::VectorInterpolateImageFunction<FieldType> FieldInterpolatorType;
  typedef FieldInterpolatorType::Pointer FieldInterpolatorPointer;
  typedef itk::InterpolateImageFunction<MaskType> MaskInterpolatorType;
  typedef MaskInterpolatorType::Pointer MaskInterpolatorPointer;

  typedef FieldInterpolatorType::PointType PointType;
  typedef FieldInterpolatorType::OutputType OutputType;

  FieldInterpolatorPointer m_fieldInterpolator;
  MaskInterpolatorPointer  m_maskInterpolator;

  std::string PrepareTagSize(FieldConstPointer) const;
  std::string PrepareTagStart(FieldConstPointer) const;
  std::string PrepareTagSpacing(FieldConstPointer) const;
  std::string PrepareTagOrigin(FieldConstPointer) const;
  std::string PrepareTagData(FieldConstPointer) const;
  std::string PrepareTagMask(MaskConstPointer) const;
};

class FemTransform3d : public gmp::Transform<3>
{
public:
  typedef gmp::Transform<3> Superclass;
  typedef Superclass::tCoords tCoords;
  typedef gmp::Transform<3> TransformType;
  typedef std::shared_ptr<TransformType> TransformPointer;

  virtual ~FemTransform3d()
  {}
  FemTransform3d();

  void set_mesh(const std::shared_ptr<TMesh3d> cpmesh)
  {
    m_sharedMesh = cpmesh;
  }

  std::shared_ptr<TMesh3d> m_sharedMesh;

  void invert();

  // the following is to be applied when a mesh is loaded
  void initOctree();

  bool m_bdbg;
  bool m_signalTopology;

  TransformType* convert_to_delta() const;
protected:
  void doInput(std::istream& is);
  void doOutput(std::ostream& os) const;
  void doOwnInit();

  virtual tCoords doOwnImg(const tCoords& pt) const;
};


/*************************************************

Applies the transforms in the INDEX space

************************************************/

class VolumeMorph
{
public:
  VolumeMorph();
  ~VolumeMorph();

  typedef gmp::Transform<3> TransformType;
  typedef std::shared_ptr<TransformType> TransformPointer;
  typedef std::list<TransformPointer> TransformContainerType;
  typedef TransformType::tCoords tCoords;

  MRI*  m_template;

  int   m_interpolationType; // interpolation method
  // should have one of the values -
  // SAMPLE_NEAREST SAMPLE_TRILINEAR SAMPLE_CUBIC

  TransformContainerType m_transforms;

  MRI* convert_transforms() const;

  // if true, the following option will cache a volume with
  // the VF of the images
  MRI* apply_transforms(MRI*,
                        bool cacheField=false,
                        const VG* vgOutput=NULL) const;

  //MRI*,
  //    bool cacheField=false,
  //    const VG* vgOutput=NULL) const;

  MRIS* apply_transforms(MRIS*) const;

  const MRI* cache() const
  {
    return mriCache;
  }

  void save(const char* fname);

  // second param is there because of the ZLib and some huge meshes....
  void load(const char* fname, unsigned int bufferMultiplier = 5,
            bool clearExisting=true);

  // needs the mriCache
  //
  // since this is meant to be used by mri_nl_align
  // will need to use the bounding box to corect for this
  GCA_MORPH* exportGcam(MRI* moving=NULL, bool useTemplateBoundingBox=false, int thresh=0, int padding=1) const;

  // extract all initial transforms and place them in the regular list
  void serialize();

  // inverts the flow: this has 2 implications:
  // 1. each of the transforms is inverted
  // 2. the order of the transforms in the chain is inverted
  void invert();

  // apply the morph to a point
  tCoords image(const tCoords& pt) const;

  //----------
  // vol geom

  void set_volGeom_fixed(const VOL_GEOM& vg)
  {
    m_vgFixed = vg;
  }
  void set_volGeom_fixed(const MRI* mri)
  {
    getVolGeom(mri, &m_vgFixed);
  }
  void set_volGeom_moving(const MRI* mri)
  {
    getVolGeom(mri, &m_vgMoving);
  }
  void set_volGeom_moving(const VOL_GEOM& vg)
  {
    m_vgMoving = vg;
  }

  const VOL_GEOM& vgFixed() const
  {
    return m_vgFixed;
  }
  const VOL_GEOM& vgMoving() const
  {
    return m_vgMoving;
  }

protected:
  // purposely not implemented
  VolumeMorph(const VolumeMorph&);
  VolumeMorph& operator=(const VolumeMorph&);

  enum IoTags
  {
    tagVgFixed = 1,
    tagVgMoving,
    tagTransform
  };

private:
  VOL_GEOM m_vgFixed;   // aka, atlas, template
  VOL_GEOM m_vgMoving;  // aka, subject

  mutable MRI* mriCache;

  void load_old(const char* fname, unsigned int bufferMultiplier = 5,
                bool clearExisting = true);
  void load_new(const char* fname, unsigned int bufferMultiplier = 5,
                bool clearExisting = true);

  std::string PrepareTagVolGeom(const VOL_GEOM& vg);
  void ReadTagVolGeom(const std::string& strData, VOL_GEOM& vg);
};

}

typedef std::shared_ptr<gmp::Transform<3> > Transform3SPointer;
typedef std::shared_ptr<gmp::FemTransform3d> FemTransform3SPointer;
typedef std::shared_ptr<gmp::DeltaTransform3d> DeltaTransform3SPointer;
typedef std::shared_ptr<gmp::VolumeMorph> VolumeMorphSPointer;

#endif
