#ifndef __kvlAtlasMeshRasterizor_txx
#define __kvlAtlasMeshRasterizor_txx

#include "GL/gl.h"


#if 1
  // Better precision
  #define TRIANGLE_WALK_DOUBLE 1
  #define CEILF(x)   ceilf(x)
  #define FLOORF(x)  floorf(x)
  #define FABSF(x)   fabsf(x)
#endif




/**
 * From imports.h
 */ 
 
 
/***
 *** IROUND: return (as an integer) float rounded to nearest integer
 ***/
#if defined(USE_SPARC_ASM) && defined(__GNUC__) && defined(__sparc__)
static INLINE int iround(float f)
{
   int r;
   __asm__ ("fstoi %1, %0" : "=f" (r) : "f" (f));
   return r;
}
#define IROUND(x)  iround(x)
#elif defined(USE_X86_ASM) && defined(__GNUC__) && defined(__i386__) && \
			(!defined(__BEOS__) || (__GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 95)))
static INLINE int iround(float f)
{
   int r;
   __asm__ ("fistpl %0" : "=m" (r) : "t" (f) : "st");
   return r;
}
#define IROUND(x)  iround(x)
#elif defined(USE_X86_ASM) && defined(__MSC__) && defined(__WIN32__)
static INLINE int iround(float f)
{
   int r;
   _asm {
	 fld f
	 fistp r
	}
   return r;
}
#define IROUND(x)  iround(x)
#elif defined(__WATCOMC__) && defined(__386__)
long iround(float f);
#pragma aux iround =                    \
	"push   eax"                        \
	"fistp  dword ptr [esp]"            \
	"pop    eax"                        \
	parm [8087]                         \
	value [eax]                         \
	modify exact [eax];
#define IROUND(x)  iround(x)
#else
#define IROUND(f)  ((int) (((f) >= 0.0F) ? ((f) + 0.5F) : ((f) - 0.5F)))
#endif



/***
 *** IS_INF_OR_NAN: test if float is infinite or NaN
 ***/
#ifdef USE_IEEE
static INLINE int IS_INF_OR_NAN( float x )
{
   fi_type tmp;
   tmp.f = x;
   return !(int)((unsigned int)((tmp.i & 0x7fffffff)-0x7f800000) >> 31);
}
#elif defined(isfinite)
#define IS_INF_OR_NAN(x)        (!isfinite(x))
#elif defined(finite)
#define IS_INF_OR_NAN(x)        (!finite(x))
#elif defined(__VMS)
#define IS_INF_OR_NAN(x)        (!finite(x))
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define IS_INF_OR_NAN(x)        (!isfinite(x))
#else
#define IS_INF_OR_NAN(x)        (!finite(x))
#endif




/** In order to trick Mesa's Triangle Rasterizer Template into compiling, we need to 
 *  make sure the following statements can be made regarding "GLcontext* ctx": 
 * - GLcontext* ctx
 * - GLfloat bf = SWRAST_CONTEXT(ctx)->_BackfaceSign;
 * - span.facing = ctx->_Facing;
 * - if (ctx->Light.ShadeModel == GL_SMOOTH)
 * - ASSERT (ctx->Light.ShadeModel == GL_FLAT);
 */
struct ContextLight
{
  unsigned long  ShadeModel;
  
  ContextLight()
    {
    ShadeModel = GL_SMOOTH;   
    }

};


struct GLcontext
{
  GLfloat  _BackfaceSign;
  GLuint   _Facing; 
  ContextLight  Light; 

  GLcontext()
    {
    _BackfaceSign = -1.0;
    _Facing = 0;
    }
  
};

#define SWRAST_CONTEXT(ctx) ctx
#define ASSERT




/**
  * From mtypes.h
  */
typedef GLfloat GLchan;
#define CHAN_MAX 1.0
#define CHAN_MAXF 1.0F
#define CHAN_TYPE GL_FLOAT
typedef int GLfixed;
#define FIXED_FRAC_BITS 11
#define FIXED_SHIFT     FIXED_FRAC_BITS
#define FIXED_ONE       (1 << FIXED_SHIFT)
#define FIXED_SCALE     ((float) FIXED_ONE)
#define FloatToFixed(X) (IROUND((X) * FIXED_SCALE))
#define FixedToFloat(X) ((X) * (1.0F / FIXED_SCALE))
#define FIXED_EPSILON   1
#define FIXED_FRAC_MASK (FIXED_ONE - 1)
#define FIXED_INT_MASK  (~FIXED_FRAC_MASK)
#define FixedCeil(X)    (((X) + FIXED_ONE - FIXED_EPSILON) & FIXED_INT_MASK)
#define FixedToInt(X)   ((X) >> FIXED_SHIFT)
#define SignedFloatToFixed(X)   FloatToFixed(X)
#define FixedFloor(X)   ((X) & FIXED_INT_MASK)


/**
 * From config.h
 */
#define SUB_PIXEL_BITS 4
#define RCOMP 0
#define GCOMP 1
#define BCOMP 2
#define ACOMP 3



/**
 * From swrast.h, keeping only the necessary components
 */
struct SWvertex {
   /** win[0], win[1] are the screen-coords of SWvertex. win[2] is the
    * z-coord. what is win[3]? */
   GLfloat win[4];
   GLchan color[4];
};



/**
 * From s_context.h
 */  
struct sw_span {
   GLint x, y;

   /** Only need to process pixels between start <= i < end */
   /** At this time, start is always zero. */
   GLuint start, end;

   /** This flag indicates that mask[] array is effectively filled with ones */
   GLboolean writeAll;

   /** either GL_POLYGON, GL_LINE, GL_POLYGON, GL_BITMAP */
   GLenum primitive;

   /** 0 = front-facing span, 1 = back-facing span (for two-sided stencil) */
   GLuint facing;

   /**
    * This bitmask (of  \link SpanFlags SPAN_* flags\endlink) indicates
    * which of the x/xStep variables are relevant.
    */
   GLuint interpMask;

   /* For horizontal spans, step is the partial derivative wrt X.
    * For lines, step is the delta from one fragment to the next.
    */
#if CHAN_TYPE == GL_FLOAT
   GLfloat red, redStep;
   GLfloat green, greenStep;
   GLfloat blue, blueStep;
   GLfloat alpha, alphaStep;
   GLfloat specRed, specRedStep;
   GLfloat specGreen, specGreenStep;
   GLfloat specBlue, specBlueStep;
#else /* CHAN_TYPE == GL_UNSIGNED_BYTE or GL_UNSIGNED_SHORT */
   GLfixed red, redStep;
   GLfixed green, greenStep;
   GLfixed blue, blueStep;
   GLfixed alpha, alphaStep;
   GLfixed specRed, specRedStep;
   GLfixed specGreen, specGreenStep;
   GLfixed specBlue, specBlueStep;
#endif
   GLfixed index, indexStep;
   GLfixed z, zStep;

   /* partial derivatives wrt X and Y. */
   GLfloat dzdx, dzdy;
   GLfloat w, dwdx, dwdy;
   GLfloat drdx, drdy;
   GLfloat dgdx, dgdy;
   GLfloat dbdx, dbdy;
   GLfloat dadx, dady;
   GLfloat dsrdx, dsrdy;
   GLfloat dsgdx, dsgdy;
   GLfloat dsbdx, dsbdy;
   GLfloat dfogdx, dfogdy;

   /**
    * This bitmask (of \link SpanFlags SPAN_* flags\endlink) indicates
    * which of the fragment arrays in the span_arrays struct are relevant.
    */
   GLuint arrayMask;

};


#define INIT_SPAN(S, PRIMITIVE, END, INTERP_MASK, ARRAY_MASK)	\
do {								\
   (S).primitive = (PRIMITIVE);					\
   (S).interpMask = (INTERP_MASK);				\
   (S).arrayMask = (ARRAY_MASK);				\
   (S).start = 0;						\
   (S).end = (END);						\
   (S).facing = 0;						\
} while (0)

#define SPAN_RGBA         0x001
#define SPAN_FLAT         0x400  /**< flat shading? */


#if CHAN_TYPE == GL_FLOAT
#define ChanToFixed(X)  (X)
#define FixedToChan(X)  (X)
#else
#define ChanToFixed(X)  IntToFixed(X)
#define FixedToChan(X)  FixedToInt(X)
#endif




#if CHAN_TYPE == GL_FLOAT
  typedef GLfloat  ChanType;
#else
  typedef GLfixed  ChanType;
#endif




/*
 * ColorTemp is used for intermediate color values.
 */
#if CHAN_TYPE == GL_FLOAT
#define ColorTemp GLfloat
#else
#define ColorTemp GLint  /* same as GLfixed */
#endif


/*
 * Walk triangle edges with GLfixed or GLdouble
 */
#if TRIANGLE_WALK_DOUBLE
#define GLinterp        GLdouble
#define InterpToInt(X)  ((GLint) (X))
#define INTERP_ONE      1.0
#else
#define GLinterp        GLfixed
#define InterpToInt(X)  FixedToInt(X)
#define INTERP_ONE      FIXED_ONE
#endif





namespace kvl
{


//
//
//
template < class TFragmentProcessor >
AtlasMeshRasterizor< TFragmentProcessor >
::AtlasMeshRasterizor()
{
  m_LabelImage = 0;

}




//
//
//
template < class TFragmentProcessor >
void
AtlasMeshRasterizor< TFragmentProcessor >
::Rasterize( const AtlasMesh* mesh )
{

  //
  if ( !m_LabelImage )
    {
    itkExceptionMacro( "No label image available!" );
    }

  //
  m_FragmentProcessor.SetMesh( mesh );

  // Rasterize all triangles
  typedef AtlasMesh::CellsContainer  CellsContainerType;
  typedef AtlasMesh::CellType  CellType;
  typedef AtlasMesh::PointType  PointType;

  CellsContainerType::ConstIterator  cellIt = mesh->GetCells()->Begin();
  CellsContainerType::ConstIterator  cellEnd = mesh->GetCells()->End();
  while ( cellIt != cellEnd )
    {
    CellType*  cell = cellIt.Value();

    if( cell->GetType() == CellType::TETRAHEDRON_CELL )
      {
      // OK, found tetrahedron. Warn the FragmentProcessor
      if ( !m_FragmentProcessor.StartNewTetrahedron( cellIt.Index() ) )
        {
        // Something is wrong with this tetrahedron; abort rasterization.
        return;
        }
      
      // Retrieve position of vertices
      CellType::PointIdIterator  pit = cell->PointIdsBegin();
      PointType  p0;
      PointType  p1;
      PointType  p2;
      PointType  p3;
      mesh->GetPoint( *pit, &p0 );
      ++pit;
      mesh->GetPoint( *pit, &p1 );
      ++pit;
      mesh->GetPoint( *pit, &p2 );
      ++pit;
      mesh->GetPoint( *pit, &p3 );

      // Rasterize
      this->RasterizeTetrahedron( p0, p1, p2, p3 );
      }

    ++cellIt;
    }
  



}



struct sortableItem
  {
  sortableItem( const itk::Vector< float, 3 >& point, const itk::Vector< float, 4 >& pis )
    {
    m_Point = point;
    m_Pis = pis;
    } 

  itk::Vector< float, 3 >  m_Point;
  itk::Vector< float, 4 >  m_Pis;
  bool operator<(const sortableItem& it) const
    {
    return m_Point[ 2 ] > it.m_Point[ 2 ]; // Needed for std::sort()
    }
  };


//
//
//
template < class TFragmentProcessor >
void
AtlasMeshRasterizor< TFragmentProcessor >
::RasterizeTetrahedron( const AtlasMesh::PointType& p0, const AtlasMesh::PointType& p1,  const AtlasMesh::PointType& p2, const AtlasMesh::PointType& p3 )
{

  // The notation (A,B,C,O) for tetrahedra and (a,b,c) and (b,c,d) for triangles is taken from figure 6 in
  //     "Voxelization of solids using simplical coverings", A.J. Rueda, R.J. Segura, F.R. Feito, J.R. de Miras, C. Ogayar
  // 

  // std::cout << "Rasterizing tetrahedron with" << std::endl;
  // std::cout << "   p0: " << p0 << std::endl;
  // std::cout << "   p1: " << p1 << std::endl;
  // std::cout << "   p2: " << p2 << std::endl;
  // std::cout << "   p3: " << p3 << std::endl;

  // Sort the vertices according to their z-coordinates: A >= B >= C >= O 
  itk::Vector< float, 4 >  pisIn0( 0.0f );  pisIn0[ 0 ] = 1.0f;
  itk::Vector< float, 4 >  pisIn1( 0.0f );  pisIn1[ 1 ] = 1.0f;
  itk::Vector< float, 4 >  pisIn2( 0.0f );  pisIn2[ 2 ] = 1.0f;
  itk::Vector< float, 4 >  pisIn3( 0.0f );  pisIn3[ 3 ] = 1.0f;


  std::vector< sortableItem >  sortableItems;
  sortableItems.push_back( sortableItem( p0.GetDataPointer(), pisIn0 ) );
  sortableItems.push_back( sortableItem( p1.GetDataPointer(), pisIn1 ) );
  sortableItems.push_back( sortableItem( p2.GetDataPointer(), pisIn2 ) );
  sortableItems.push_back( sortableItem( p3.GetDataPointer(), pisIn3 ) );
  std::sort( sortableItems.begin(), sortableItems.end() );

  const itk::Vector< float, 3 >  pA = sortableItems[ 0 ].m_Point;
  const itk::Vector< float, 3 >  pB = sortableItems[ 1 ].m_Point;
  const itk::Vector< float, 3 >  pC = sortableItems[ 2 ].m_Point;
  const itk::Vector< float, 3 >  pO = sortableItems[ 3 ].m_Point;

  const itk::Vector< float, 4 >  pisInA = sortableItems[ 0 ].m_Pis;
  const itk::Vector< float, 4 >  pisInB = sortableItems[ 1 ].m_Pis;
  const itk::Vector< float, 4 >  pisInC = sortableItems[ 2 ].m_Pis;
  const itk::Vector< float, 4 >  pisInO = sortableItems[ 3 ].m_Pis;

  // std::cout << "Rasterizing tetrahedron with" << std::endl;
  // std::cout << "   pA: " << pA << "   pisInA: " <<  pisInA << std::endl;
  // std::cout << "   pB: " << pB << "   pisInB: " <<  pisInB << std::endl;
  // std::cout << "   pC: " << pC << "   pisInC: " <<  pisInC << std::endl;
  // std::cout << "   pO: " << pO << "   pisInO: " <<  pisInO << std::endl;



  // Loop over all zs, calculate the intersection locations and pis of the vertices of the resulting triangle, and rasterize
  const int  minZ = static_cast< int >( ceil( pO[ 2 ] ) );
  int  maxZ = static_cast< int >( floor( pA[ 2 ] ) );
  if ( static_cast< float >( maxZ ) == pA[ 2 ] )
    {
    // The upper limit is exactly on a Z-plane. To avoid that the voxels in this plane are rasterized mulitple times, 
    // our convention is that these voxels don't belong to the current tetrahedron 
    maxZ--;
    }
  //std::cout << "\n\n   Drawing tetrahedron from " << minZ << " to " << maxZ << std::endl;

  
  for ( int z = maxZ; z >= minZ; z-- )
    {
    // Determine the position and pis in a. a always lies on the edge AO no matter what
    const float  fractionalDistanceOnAO = ( z - pO[ 2 ] + 1e-10 ) / ( pA[ 2 ] - pO[ 2 ] + 1e-10 );
    const itk::Vector< float, 3 >  pa = pA * fractionalDistanceOnAO + ( 1 - fractionalDistanceOnAO ) * pO;
    const itk::Vector< float, 4 >  pisIna = fractionalDistanceOnAO * pisInA + ( 1 - fractionalDistanceOnAO ) * pisInO;
  
    // Determine the position of b. If z > B, b lies on AB; otherwise on BO
    itk::Vector< float, 3 >  pb;
    itk::Vector< float, 4 >  pisInb;
    if ( z > pB[ 2 ] )
      {
      const float  fractionalDistanceOnAB = ( z - pB[ 2 ] + 1e-10 ) / ( pA[ 2 ] - pB[ 2 ] + 1e-10 );
      pb = fractionalDistanceOnAB * pA + ( 1 - fractionalDistanceOnAB ) * pB;
      pisInb = fractionalDistanceOnAB * pisInA + ( 1 - fractionalDistanceOnAB ) * pisInB;
      }
    else
      {
      const float  fractionalDistanceOnBO = ( z - pO[ 2 ] + 1e-10 ) / ( pB[ 2 ] - pO[ 2 ] + 1e-10 );
      pb = fractionalDistanceOnBO * pB + ( 1 - fractionalDistanceOnBO ) * pO;
      pisInb = fractionalDistanceOnBO * pisInB + ( 1 - fractionalDistanceOnBO ) * pisInO;
      }

    // Determine the position of c. If z > C, c lies on AC; otherwise on CO
    itk::Vector< float, 3 >  pc;
    itk::Vector< float, 4 >  pisInc;
    if ( z > pC[ 2 ] )
      {
      const float  fractionalDistanceOnAC = ( z - pC[ 2 ] + 1e-10 ) / ( pA[ 2 ] - pC[ 2 ] + 1e-10 );
      pc = fractionalDistanceOnAC * pA + ( 1 - fractionalDistanceOnAC ) * pC;
      pisInc = fractionalDistanceOnAC * pisInA + ( 1 - fractionalDistanceOnAC ) * pisInC;
      }
    else
      {
      const float  fractionalDistanceOnCO = ( z - pO[ 2 ] + 1e-10 ) / ( pC[ 2 ] - pO[ 2 ] + 1e-10 );
      pc = fractionalDistanceOnCO * pC + ( 1 - fractionalDistanceOnCO ) * pO;
      pisInc = fractionalDistanceOnCO * pisInC + ( 1 - fractionalDistanceOnCO ) * pisInO;
      }

    // Rasterize the triangle (a,b,c)
    // std::cout << "\n\n      Rasterizing triangle in plane " << z << std::endl;
    // std::cout << "             pa: " << pa << "   pisIna: " <<  pisIna << std::endl;
    // std::cout << "             pb: " << pb << "   pisInb: " <<  pisInb << std::endl;
    // std::cout << "             pc: " << pc << "   pisInc: " <<  pisInc << std::endl;
    this->RasterizeTriangle( pa.GetDataPointer(), pb.GetDataPointer(), pc.GetDataPointer(),
                             pisIna.GetDataPointer(), pisInb.GetDataPointer(), pisInc.GetDataPointer(), 
                             z );


    // If B > z > C, there is an extra triangle (b,c,d) to rasterize. 
    if ( ( z > pC[ 2 ] ) && ( z < pB[ 2 ] ) )
      {
      // Determine the position and pis in d. d always lies on the edge BC no matter what
      const float  fractionalDistanceOnBC = ( z - pC[ 2 ] + 1e-10 ) / ( pB[ 2 ] - pC[ 2 ] + 1e-10 );
      const itk::Vector< float, 3 >  pd = fractionalDistanceOnBC * pB + ( 1 - fractionalDistanceOnBC ) * pC;
      const itk::Vector< float, 4 >  pisInd = fractionalDistanceOnBC * pisInB + ( 1 - fractionalDistanceOnBC ) * pisInC;

      // Rasterize the triangle (b,c,d)
      //std::cout << "\n\n      Rasterizing extra triangle in plane " << z << std::endl;
      //std::cout << "             pb: " << pb << "   pisInb: " <<  pisInb << std::endl;
      //std::cout << "             pc: " << pc << "   pisInc: " <<  pisInc << std::endl;
      //std::cout << "             pd: " << pd << "   pisInd: " <<  pisInd << std::endl;
      this->RasterizeTriangle( pb.GetDataPointer(), pc.GetDataPointer(), pd.GetDataPointer(),
                              pisInb.GetDataPointer(), pisInc.GetDataPointer(), pisInd.GetDataPointer(), 
                              z );
      } 


    }


}


//
//
//
template < class TFragmentProcessor >
void
AtlasMeshRasterizor< TFragmentProcessor >
::RasterizeTriangle( const float* vertex0, const float* vertex1, const float* vertex2, 
                     const float* pisIn0, const float* pisIn1, const float* pisIn2, 
                     const int zLevel )
{


  /** TODO
   *    - really clean up code
   *    - make rgba values in vertices input rather than hardcoded to red for vertex0, green for vertex1, blue for vertex2
   *    - make AtlasMeshRasterizor templated over pixel type
   *    - make AtlasMeshRasterizor work without there being an image
   */


  /** 
   * Some initial glue code; to be removed later on
   */
#if 0 
  GLcontext*  ctx = new GLcontext();
#else
  GLcontext  contextInstance;
  GLcontext*  ctx = &contextInstance;
#endif

#if 0
  SWvertex*  v0 = new SWvertex();
#else
  SWvertex  v0Instance;
  SWvertex*  v0 = &v0Instance;
#endif
  const float  gedverX = 0.5 - 1e-7;
  const float  gedverY = 0.5 - 1e-7;

  v0->win[ 0 ] = vertex0[ 0 ] + gedverX;
  v0->win[ 1 ] = vertex0[ 1 ] + gedverY;
  v0->win[ 2 ] = 0;
  v0->win[ 3 ] = 0;
  v0->color[ 0 ] = pisIn0[ 0 ];
  v0->color[ 1 ] = pisIn0[ 1 ];
  v0->color[ 2 ] = pisIn0[ 2 ];
  v0->color[ 3 ] = pisIn0[ 3 ];
#if 0
  SWvertex*  v1 = new SWvertex();
#else
  SWvertex  v1Instance;
  SWvertex*  v1 = &v1Instance;
#endif
  v1->win[ 0 ] = vertex1[ 0 ] + gedverX;
  v1->win[ 1 ] = vertex1[ 1 ] + gedverY;
  v1->win[ 2 ] = 0;
  v1->win[ 3 ] = 0;
  v1->color[ 0 ] = pisIn1[ 0 ];
  v1->color[ 1 ] = pisIn1[ 1 ];
  v1->color[ 2 ] = pisIn1[ 2 ];
  v1->color[ 3 ] = pisIn1[ 3 ];
#if 0
  SWvertex*  v2 = new SWvertex(); 
#else
  SWvertex  v2Instance;
  SWvertex*  v2 = &v2Instance;
#endif
  v2->win[ 0 ] = vertex2[ 0 ] + gedverX;
  v2->win[ 1 ] = vertex2[ 1 ] + gedverY;
  v2->win[ 2 ] = 0;
  v2->win[ 3 ] = 0;
  v2->color[ 0 ] = pisIn2[ 0 ];
  v2->color[ 1 ] = pisIn2[ 1 ];
  v2->color[ 2 ] = pisIn2[ 2 ];
  v2->color[ 3 ] = pisIn2[ 3 ];



  

  typedef struct {
      const SWvertex *v0, *v1;   /* Y(v0) < Y(v1) */
#if TRIANGLE_WALK_DOUBLE
      GLdouble dx;	/* X(v1) - X(v0) */
      GLdouble dy;	/* Y(v1) - Y(v0) */
      GLdouble dxdy;	/* dx/dy */
      GLdouble adjy;	/* adjust from v[0]->fy to fsy, scaled */
      GLdouble fsx;	/* first sample point x coord */
      GLdouble fsy;
      GLdouble fx0;	/*X of lower endpoint */
#else
      GLfloat dx;	/* X(v1) - X(v0) */
      GLfloat dy;	/* Y(v1) - Y(v0) */
      GLfloat dxdy;	/* dx/dy */
      GLfixed fdxdy;	/* dx/dy in fixed-point */
      GLfloat adjy;	/* adjust from v[0]->fy to fsy, scaled */
      GLfixed fsx;	/* first sample point x coord */
      GLfixed fsy;
      GLfixed fx0;	/* fixed pt X of lower endpoint */
#endif
      GLint lines;	/* number of lines to be sampled on this edge */
   } EdgeT;

   
   EdgeT eMaj, eTop, eBot;
   GLfloat oneOverArea;
   const SWvertex *vMin, *vMid, *vMax;  /* Y(vMin)<=Y(vMid)<=Y(vMax) */
   GLfloat bf = SWRAST_CONTEXT(ctx)->_BackfaceSign;
#if !TRIANGLE_WALK_DOUBLE
   const GLint snapMask = ~((FIXED_ONE / (1 << SUB_PIXEL_BITS)) - 1); /* for x/y coord snapping */
#endif
   GLinterp vMin_fx, vMin_fy, vMid_fx, vMid_fy, vMax_fx, vMax_fy;

   struct sw_span span;

   INIT_SPAN(span, GL_POLYGON, 0, 0, 0);


   /*
   printf("%s()\n", __FUNCTION__);
   printf("  %g, %g, %g\n", v0->win[0], v0->win[1], v0->win[2]);
   printf("  %g, %g, %g\n", v1->win[0], v1->win[1], v1->win[2]);
   printf("  %g, %g, %g\n", v2->win[0], v2->win[1], v2->win[2]);
   */
   /*
   ASSERT(v0->win[2] >= 0.0);
   ASSERT(v1->win[2] >= 0.0);
   ASSERT(v2->win[2] >= 0.0);
   */
   /* Compute fixed point x,y coords w/ half-pixel offsets and snapping.
    * And find the order of the 3 vertices along the Y axis.
    */
   {
#if TRIANGLE_WALK_DOUBLE
      const GLdouble fy0 = v0->win[1] - 0.5;
      const GLdouble fy1 = v1->win[1] - 0.5;
      const GLdouble fy2 = v2->win[1] - 0.5;
#else
      const GLfixed fy0 = FloatToFixed(v0->win[1] - 0.5F) & snapMask;
      const GLfixed fy1 = FloatToFixed(v1->win[1] - 0.5F) & snapMask;
      const GLfixed fy2 = FloatToFixed(v2->win[1] - 0.5F) & snapMask;
#endif
      if (fy0 <= fy1) {
         if (fy1 <= fy2) {
            /* y0 <= y1 <= y2 */
            vMin = v0;   vMid = v1;   vMax = v2;
            vMin_fy = fy0;  vMid_fy = fy1;  vMax_fy = fy2;
         }
         else if (fy2 <= fy0) {
            /* y2 <= y0 <= y1 */
            vMin = v2;   vMid = v0;   vMax = v1;
            vMin_fy = fy2;  vMid_fy = fy0;  vMax_fy = fy1;
         }
         else {
            /* y0 <= y2 <= y1 */
            vMin = v0;   vMid = v2;   vMax = v1;
            vMin_fy = fy0;  vMid_fy = fy2;  vMax_fy = fy1;
            bf = -bf;
         }
      }
      else {
         if (fy0 <= fy2) {
            /* y1 <= y0 <= y2 */
            vMin = v1;   vMid = v0;   vMax = v2;
            vMin_fy = fy1;  vMid_fy = fy0;  vMax_fy = fy2;
            bf = -bf;
         }
         else if (fy2 <= fy1) {
            /* y2 <= y1 <= y0 */
            vMin = v2;   vMid = v1;   vMax = v0;
            vMin_fy = fy2;  vMid_fy = fy1;  vMax_fy = fy0;
            bf = -bf;
         }
         else {
            /* y1 <= y2 <= y0 */
            vMin = v1;   vMid = v2;   vMax = v0;
            vMin_fy = fy1;  vMid_fy = fy2;  vMax_fy = fy0;
         }
      }

      /* fixed point X coords */
#if TRIANGLE_WALK_DOUBLE
      vMin_fx = vMin->win[0] + 0.5;
      vMid_fx = vMid->win[0] + 0.5;
      vMax_fx = vMax->win[0] + 0.5;
#else
      vMin_fx = FloatToFixed(vMin->win[0] + 0.5F) & snapMask;
      vMid_fx = FloatToFixed(vMid->win[0] + 0.5F) & snapMask;
      vMax_fx = FloatToFixed(vMax->win[0] + 0.5F) & snapMask;
#endif
   }

   /* vertex/edge relationship */
   eMaj.v0 = vMin;   eMaj.v1 = vMax;   /*TODO: .v1's not needed */
   eTop.v0 = vMid;   eTop.v1 = vMax;
   eBot.v0 = vMin;   eBot.v1 = vMid;

   /* compute deltas for each edge:  vertex[upper] - vertex[lower] */
#if TRIANGLE_WALK_DOUBLE
   eMaj.dx = vMax_fx - vMin_fx;
   eMaj.dy = vMax_fy - vMin_fy;
   eTop.dx = vMax_fx - vMid_fx;
   eTop.dy = vMax_fy - vMid_fy;
   eBot.dx = vMid_fx - vMin_fx;
   eBot.dy = vMid_fy - vMin_fy;
#else
   eMaj.dx = FixedToFloat(vMax_fx - vMin_fx);
   eMaj.dy = FixedToFloat(vMax_fy - vMin_fy);
   eTop.dx = FixedToFloat(vMax_fx - vMid_fx);
   eTop.dy = FixedToFloat(vMax_fy - vMid_fy);
   eBot.dx = FixedToFloat(vMid_fx - vMin_fx);
   eBot.dy = FixedToFloat(vMid_fy - vMin_fy);
#endif

   /* compute area, oneOverArea and perform backface culling */
   {
#if TRIANGLE_WALK_DOUBLE
      const GLdouble area = eMaj.dx * eBot.dy - eBot.dx * eMaj.dy;
#else
      const GLfloat area = eMaj.dx * eBot.dy - eBot.dx * eMaj.dy;
#endif

#if 0
      /* Do backface culling */
      if (area * bf < 0.0)
         return;
#endif

      if (IS_INF_OR_NAN(area) || area == 0.0F)
         return;

      oneOverArea = 1.0F / area;
   }

   span.facing = ctx->_Facing; /* for 2-sided stencil test */

   /* Edge setup.  For a triangle strip these could be reused... */
   {
#if TRIANGLE_WALK_DOUBLE
      eMaj.fsy = CEILF(vMin_fy);
      eMaj.lines = (GLint) CEILF(vMax_fy - eMaj.fsy);
#else
      eMaj.fsy = FixedCeil(vMin_fy);
      eMaj.lines = FixedToInt(FixedCeil(vMax_fy - eMaj.fsy));
#endif
      if (eMaj.lines > 0) {
         eMaj.dxdy = eMaj.dx / eMaj.dy;
#if TRIANGLE_WALK_DOUBLE
         eMaj.adjy = (eMaj.fsy - vMin_fy) * FIXED_SCALE;  /* SCALED! */
         eMaj.fx0 = vMin_fx;
         eMaj.fsx = eMaj.fx0 + (eMaj.adjy * eMaj.dxdy) / (GLdouble) FIXED_SCALE;
#else
         eMaj.fdxdy = SignedFloatToFixed(eMaj.dxdy);
         eMaj.adjy = (GLfloat) (eMaj.fsy - vMin_fy);  /* SCALED! */
         eMaj.fx0 = vMin_fx;
         eMaj.fsx = eMaj.fx0 + (GLfixed) (eMaj.adjy * eMaj.dxdy);
#endif
      }
      else {
         return;  /*CULLED*/
      }

#if TRIANGLE_WALK_DOUBLE
      eTop.fsy = CEILF(vMid_fy);
      eTop.lines = (GLint) CEILF(vMax_fy - eTop.fsy);
#else
      eTop.fsy = FixedCeil(vMid_fy);
      eTop.lines = FixedToInt(FixedCeil(vMax_fy - eTop.fsy));
#endif
      if (eTop.lines > 0) {
         eTop.dxdy = eTop.dx / eTop.dy;
#if TRIANGLE_WALK_DOUBLE
         eTop.adjy = (eTop.fsy - vMid_fy) * FIXED_SCALE; /* SCALED! */
         eTop.fx0 = vMid_fx;
         eTop.fsx = eTop.fx0 + (eTop.adjy * eTop.dxdy) / (GLdouble) FIXED_SCALE;
#else
         eTop.fdxdy = SignedFloatToFixed(eTop.dxdy);
         eTop.adjy = (GLfloat) (eTop.fsy - vMid_fy); /* SCALED! */
         eTop.fx0 = vMid_fx;
         eTop.fsx = eTop.fx0 + (GLfixed) (eTop.adjy * eTop.dxdy);
#endif
      }

#if TRIANGLE_WALK_DOUBLE
      eBot.fsy = CEILF(vMin_fy);
      eBot.lines = (GLint) CEILF(vMid_fy - eBot.fsy);
#else
      eBot.fsy = FixedCeil(vMin_fy);
      eBot.lines = FixedToInt(FixedCeil(vMid_fy - eBot.fsy));
#endif
      if (eBot.lines > 0) {
         eBot.dxdy = eBot.dx / eBot.dy;
#if TRIANGLE_WALK_DOUBLE
         eBot.adjy = (eBot.fsy - vMin_fy) * FIXED_SCALE;  /* SCALED! */
         eBot.fx0 = vMin_fx;
         eBot.fsx = eBot.fx0 + (eBot.adjy * eBot.dxdy) / (GLdouble) FIXED_SCALE;
#else
         eBot.fdxdy = SignedFloatToFixed(eBot.dxdy);
         eBot.adjy = (GLfloat) (eBot.fsy - vMin_fy);  /* SCALED! */
         eBot.fx0 = vMin_fx;
         eBot.fsx = eBot.fx0 + (GLfixed) (eBot.adjy * eBot.dxdy);
#endif
      }
   }

   /*
    * Conceptually, we view a triangle as two subtriangles
    * separated by a perfectly horizontal line.  The edge that is
    * intersected by this line is one with maximal absolute dy; we
    * call it a ``major'' edge.  The other two edges are the
    * ``top'' edge (for the upper subtriangle) and the ``bottom''
    * edge (for the lower subtriangle).  If either of these two
    * edges is horizontal or very close to horizontal, the
    * corresponding subtriangle might cover zero sample points;
    * we take care to handle such cases, for performance as well
    * as correctness.
    *
    * By stepping rasterization parameters along the major edge,
    * we can avoid recomputing them at the discontinuity where
    * the top and bottom edges meet.  However, this forces us to
    * be able to scan both left-to-right and right-to-left.
    * Also, we must determine whether the major edge is at the
    * left or right side of the triangle.  We do this by
    * computing the magnitude of the cross-product of the major
    * and top edges.  Since this magnitude depends on the sine of
    * the angle between the two edges, its sign tells us whether
    * we turn to the left or to the right when travelling along
    * the major edge to the top edge, and from this we infer
    * whether the major edge is on the left or the right.
    *
    * Serendipitously, this cross-product magnitude is also a
    * value we need to compute the iteration parameter
    * derivatives for the triangle, and it can be used to perform
    * backface culling because its sign tells us whether the
    * triangle is clockwise or counterclockwise.  In this code we
    * refer to it as ``area'' because it's also proportional to
    * the pixel area of the triangle.
    */

   {
      GLint scan_from_left_to_right;  /* true if scanning left-to-right */

      scan_from_left_to_right = (oneOverArea < 0.0F);


      /* compute d?/dx and d?/dy derivatives */
      span.interpMask |= SPAN_RGBA;
      if (ctx->Light.ShadeModel == GL_SMOOTH) {
         GLfloat eMaj_dr = (GLfloat) ((ColorTemp) vMax->color[RCOMP] - (ColorTemp) vMin->color[RCOMP]);
         GLfloat eBot_dr = (GLfloat) ((ColorTemp) vMid->color[RCOMP] - (ColorTemp) vMin->color[RCOMP]);
         GLfloat eMaj_dg = (GLfloat) ((ColorTemp) vMax->color[GCOMP] - (ColorTemp) vMin->color[GCOMP]);
         GLfloat eBot_dg = (GLfloat) ((ColorTemp) vMid->color[GCOMP] - (ColorTemp) vMin->color[GCOMP]);
         GLfloat eMaj_db = (GLfloat) ((ColorTemp) vMax->color[BCOMP] - (ColorTemp) vMin->color[BCOMP]);
         GLfloat eBot_db = (GLfloat) ((ColorTemp) vMid->color[BCOMP] - (ColorTemp) vMin->color[BCOMP]);
         GLfloat eMaj_da = (GLfloat) ((ColorTemp) vMax->color[ACOMP] - (ColorTemp) vMin->color[ACOMP]);
         GLfloat eBot_da = (GLfloat) ((ColorTemp) vMid->color[ACOMP] - (ColorTemp) vMin->color[ACOMP]);
         
         span.drdx = oneOverArea * (eMaj_dr * eBot.dy - eMaj.dy * eBot_dr);
         span.drdy = oneOverArea * (eMaj.dx * eBot_dr - eMaj_dr * eBot.dx);
         span.dgdx = oneOverArea * (eMaj_dg * eBot.dy - eMaj.dy * eBot_dg);
         span.dgdy = oneOverArea * (eMaj.dx * eBot_dg - eMaj_dg * eBot.dx);
         span.dbdx = oneOverArea * (eMaj_db * eBot.dy - eMaj.dy * eBot_db);
         span.dbdy = oneOverArea * (eMaj.dx * eBot_db - eMaj_db * eBot.dx);
         span.dadx = oneOverArea * (eMaj_da * eBot.dy - eMaj.dy * eBot_da);
         span.dady = oneOverArea * (eMaj.dx * eBot_da - eMaj_da * eBot.dx);
#  if CHAN_TYPE == GL_FLOAT
         span.redStep   = span.drdx;
         span.greenStep = span.dgdx;
         span.blueStep  = span.dbdx;
         span.alphaStep = span.dadx;
#  else
         span.redStep   = SignedFloatToFixed(span.drdx);
         span.greenStep = SignedFloatToFixed(span.dgdx);
         span.blueStep  = SignedFloatToFixed(span.dbdx);
         span.alphaStep = SignedFloatToFixed(span.dadx);
#  endif /* GL_FLOAT */
      }
      else {
         // ASSERT (ctx->Light.ShadeModel == GL_FLAT);
         span.interpMask |= SPAN_FLAT;
         span.drdx = span.drdy = 0.0F;
         span.dgdx = span.dgdy = 0.0F;
         span.dbdx = span.dbdy = 0.0F;
         span.dadx = span.dady = 0.0F;
#    if CHAN_TYPE == GL_FLOAT
	 span.redStep   = 0.0F;
	 span.greenStep = 0.0F;
	 span.blueStep  = 0.0F;
	 span.alphaStep = 0.0F;
#    else
	 span.redStep   = 0;
	 span.greenStep = 0;
	 span.blueStep  = 0;
	 span.alphaStep = 0;
#    endif /* GL_FLOAT */
      }

      
      /*
       * We always sample at pixel centers.  However, we avoid
       * explicit half-pixel offsets in this code by incorporating
       * the proper offset in each of x and y during the
       * transformation to window coordinates.
       *
       * We also apply the usual rasterization rules to prevent
       * cracks and overlaps.  A pixel is considered inside a
       * subtriangle if it meets all of four conditions: it is on or
       * to the right of the left edge, strictly to the left of the
       * right edge, on or below the top edge, and strictly above
       * the bottom edge.  (Some edges may be degenerate.)
       *
       * The following discussion assumes left-to-right scanning
       * (that is, the major edge is on the left); the right-to-left
       * case is a straightforward variation.
       *
       * We start by finding the half-integral y coordinate that is
       * at or below the top of the triangle.  This gives us the
       * first scan line that could possibly contain pixels that are
       * inside the triangle.
       *
       * Next we creep down the major edge until we reach that y,
       * and compute the corresponding x coordinate on the edge.
       * Then we find the half-integral x that lies on or just
       * inside the edge.  This is the first pixel that might lie in
       * the interior of the triangle.  (We won't know for sure
       * until we check the other edges.)
       *
       * As we rasterize the triangle, we'll step down the major
       * edge.  For each step in y, we'll move an integer number
       * of steps in x.  There are two possible x step sizes, which
       * we'll call the ``inner'' step (guaranteed to land on the
       * edge or inside it) and the ``outer'' step (guaranteed to
       * land on the edge or outside it).  The inner and outer steps
       * differ by one.  During rasterization we maintain an error
       * term that indicates our distance from the true edge, and
       * select either the inner step or the outer step, whichever
       * gets us to the first pixel that falls inside the triangle.
       *
       * All parameters (z, red, etc.) as well as the buffer
       * addresses for color and z have inner and outer step values,
       * so that we can increment them appropriately.  This method
       * eliminates the need to adjust parameters by creeping a
       * sub-pixel amount into the triangle at each scanline.
       */

      {
         GLint subTriangle;
         GLinterp fxLeftEdge = 0, fxRightEdge = 0;
         GLinterp fdxLeftEdge = 0, fdxRightEdge = 0;
         GLinterp fError = 0, fdError = 0;
         
         const unsigned char*  pRow = NULL;
         GLint dPRowOuter = 0, dPRowInner;  /* offset in bytes */
         
         ColorTemp rLeft = 0, fdrOuter = 0, fdrInner;
         ColorTemp gLeft = 0, fdgOuter = 0, fdgInner;
         ColorTemp bLeft = 0, fdbOuter = 0, fdbInner;
         ColorTemp aLeft = 0, fdaOuter = 0, fdaInner;

         for (subTriangle=0; subTriangle<=1; subTriangle++) {
            EdgeT *eLeft, *eRight;
            int setupLeft, setupRight;
            int lines;

            if (subTriangle==0) {
               /* bottom half */
               if (scan_from_left_to_right) {
                  eLeft = &eMaj;
                  eRight = &eBot;
                  lines = eRight->lines;
                  setupLeft = 1;
                  setupRight = 1;
               }
               else {
                  eLeft = &eBot;
                  eRight = &eMaj;
                  lines = eLeft->lines;
                  setupLeft = 1;
                  setupRight = 1;
               }
            }
            else {
               /* top half */
               if (scan_from_left_to_right) {
                  eLeft = &eMaj;
                  eRight = &eTop;
                  lines = eRight->lines;
                  setupLeft = 0;
                  setupRight = 1;
               }
               else {
                  eLeft = &eTop;
                  eRight = &eMaj;
                  lines = eLeft->lines;
                  setupLeft = 1;
                  setupRight = 0;
               }
               if (lines == 0)
                  return;
            }

            if (setupLeft && eLeft->lines > 0) {
               const SWvertex *vLower = eLeft->v0;
#if TRIANGLE_WALK_DOUBLE
               const GLdouble fsy = eLeft->fsy;
               const GLdouble fsx = eLeft->fsx;
               const GLdouble fx = CEILF(fsx);
               const GLdouble adjx = (fx - eLeft->fx0) * FIXED_SCALE;  /* SCALED! */
#else
               const GLfixed fsy = eLeft->fsy;
               const GLfixed fsx = eLeft->fsx;  /* no fractional part */
               const GLfixed fx = FixedCeil(fsx);  /* no fractional part */
               const GLfixed adjx = (GLinterp) (fx - eLeft->fx0); /* SCALED! */
#endif
               const GLinterp adjy = (GLinterp) eLeft->adjy;      /* SCALED! */
               GLint idxOuter;
#if TRIANGLE_WALK_DOUBLE
               GLdouble dxOuter;

               fError = fx - fsx - 1.0;
               fxLeftEdge = fsx;
               fdxLeftEdge = eLeft->dxdy;
               dxOuter = FLOORF(fdxLeftEdge);
               fdError = dxOuter - fdxLeftEdge + 1.0;
               idxOuter = (GLint) dxOuter;
               span.y = (GLint) fsy;
#else
               GLfloat dxOuter;
               GLfixed fdxOuter;

               fError = fx - fsx - FIXED_ONE;
               fxLeftEdge = fsx - FIXED_EPSILON;
               fdxLeftEdge = eLeft->fdxdy;
               fdxOuter = FixedFloor(fdxLeftEdge - FIXED_EPSILON);
               fdError = fdxOuter - fdxLeftEdge + FIXED_ONE;
               idxOuter = FixedToInt(fdxOuter);
               dxOuter = (GLfloat) idxOuter;
               span.y = FixedToInt(fsy);
#endif

               /* silence warnings on some compilers */
               (void) dxOuter;
               (void) adjx;
               (void) adjy;
               (void) vLower;
               
               
               
               LabelImageType::IndexType  index;
               index[ 0 ] = InterpToInt(fxLeftEdge);
               index[ 1 ] = span.y;
               index[ 2 ] = zLevel;
               pRow = m_LabelImage->GetBufferPointer() + m_LabelImage->ComputeOffset( index );
               //pRow = (PIXEL_TYPE *) PIXEL_ADDRESS( InterpToInt(fxLeftEdge), span.y );
#if 0
               dPRowOuter = -( (int) m_LabelImage->GetLargestPossibleRegion().GetSize( 0 ) ) + idxOuter;
#else
               dPRowOuter = ( (int) m_LabelImage->GetLargestPossibleRegion().GetSize( 0 ) ) + idxOuter;
#endif 
               //dPRowOuter = -((int)BYTES_PER_ROW) + idxOuter * sizeof(PIXEL_TYPE);
               /* negative because Y=0 at bottom and increases upward */
        
               
                      
               /*
                * Now we need the set of parameter (z, color, etc.) values at
                * the point (fx, fsy).  This gives us properly-sampled parameter
                * values that we can step from pixel to pixel.  Furthermore,
                * although we might have intermediate results that overflow
                * the normal parameter range when we step temporarily outside
                * the triangle, we shouldn't overflow or underflow for any
                * pixel that's actually inside the triangle.
                */

               if (ctx->Light.ShadeModel == GL_SMOOTH) {
#  if CHAN_TYPE == GL_FLOAT
                  rLeft = vLower->color[RCOMP] + (span.drdx * adjx + span.drdy * adjy) * (1.0F / FIXED_SCALE);
                  gLeft = vLower->color[GCOMP] + (span.dgdx * adjx + span.dgdy * adjy) * (1.0F / FIXED_SCALE);
                  bLeft = vLower->color[BCOMP] + (span.dbdx * adjx + span.dbdy * adjy) * (1.0F / FIXED_SCALE);
                  aLeft = vLower->color[ACOMP] + (span.dadx * adjx + span.dady * adjy) * (1.0F / FIXED_SCALE);
                  fdrOuter = span.drdy + dxOuter * span.drdx;
                  fdgOuter = span.dgdy + dxOuter * span.dgdx;
                  fdbOuter = span.dbdy + dxOuter * span.dbdx;
                  fdaOuter = span.dady + dxOuter * span.dadx;
#  else
                  rLeft = (GLint)(ChanToFixed(vLower->color[RCOMP]) + span.drdx * adjx + span.drdy * adjy) + FIXED_HALF;
                  gLeft = (GLint)(ChanToFixed(vLower->color[GCOMP]) + span.dgdx * adjx + span.dgdy * adjy) + FIXED_HALF;
                  bLeft = (GLint)(ChanToFixed(vLower->color[BCOMP]) + span.dbdx * adjx + span.dbdy * adjy) + FIXED_HALF;
                  aLeft = (GLint)(ChanToFixed(vLower->color[ACOMP]) + span.dadx * adjx + span.dady * adjy) + FIXED_HALF;
                  fdrOuter = SignedFloatToFixed(span.drdy + dxOuter * span.drdx);
                  fdgOuter = SignedFloatToFixed(span.dgdy + dxOuter * span.dgdx);
                  fdbOuter = SignedFloatToFixed(span.dbdy + dxOuter * span.dbdx);
                  fdaOuter = SignedFloatToFixed(span.dady + dxOuter * span.dadx);
#  endif
               }
               else {
                  // ASSERT (ctx->Light.ShadeModel == GL_FLAT);
#  if CHAN_TYPE == GL_FLOAT
                  rLeft = v2->color[RCOMP];
                  gLeft = v2->color[GCOMP];
                  bLeft = v2->color[BCOMP];
                  aLeft = v2->color[ACOMP];
                  fdrOuter = fdgOuter = fdbOuter = 0.0F;
                  fdaOuter = 0.0F;
#  else
                  rLeft = ChanToFixed(v2->color[RCOMP]);
                  gLeft = ChanToFixed(v2->color[GCOMP]);
                  bLeft = ChanToFixed(v2->color[BCOMP]);
                  aLeft = ChanToFixed(v2->color[ACOMP]);
                  fdrOuter = fdgOuter = fdbOuter = 0;
                  fdaOuter = 0;
#  endif
               }

            } /*if setupLeft*/


            if (setupRight && eRight->lines>0) {
#if TRIANGLE_WALK_DOUBLE
               fxRightEdge = eRight->fsx;
               fdxRightEdge = eRight->dxdy;
#else
               fxRightEdge = eRight->fsx - FIXED_EPSILON;
               fdxRightEdge = eRight->fdxdy;
#endif
            }

            if (lines==0) {
               continue;
            }


            /* Rasterize setup */
            dPRowInner = dPRowOuter + 1;
            
            fdrInner = fdrOuter + span.redStep;
            fdgInner = fdgOuter + span.greenStep;
            fdbInner = fdbOuter + span.blueStep;
            fdaInner = fdaOuter + span.alphaStep;

            while (lines > 0) {
               /* initialize the span interpolants to the leftmost value */
               /* ff = fixed-pt fragment */
               const GLint right = InterpToInt(fxRightEdge);
               span.x = InterpToInt(fxLeftEdge);
               if (right <= span.x)
                  span.end = 0;
               else
                  span.end = right - span.x;

               span.red = rLeft;
               span.green = gLeft;
               span.blue = bLeft;
               span.alpha = aLeft;

#if 0
               if (span.end > 1) {
                  /* Under rare circumstances, we might have to fudge the
                   * colors. XXX does this really happen anymore???
                   */
                  const GLint len = span.end - 1;
                  (void) len;
                  {
                     GLfixed ffrend = span.red + len * span.redStep;
                     GLfixed ffgend = span.green + len * span.greenStep;
                     GLfixed ffbend = span.blue + len * span.blueStep;
                     if (ffrend < 0) {
                        span.red -= ffrend;
                        if (span.red < 0)
                           span.red = 0;
                     }
                     if (ffgend < 0) {
                        span.green -= ffgend;
                        if (span.green < 0)
                           span.green = 0;
                     }
                     if (ffbend < 0) {
                        span.blue -= ffbend;
                        if (span.blue < 0)
                           span.blue = 0;
                     }
                  }
                  {
                     GLfixed ffaend = span.alpha + len * span.alphaStep;
                     if (ffaend < 0) {
                        span.alpha -= ffaend;
                        if (span.alpha < 0)
                           span.alpha = 0;
                     }
                  }
               } /* span.end > 1 */
#endif

               /* This is where we actually generate fragments */
               if (span.end > 0) {
               
                  m_FragmentProcessor.StartNewSpan( span.x, span.y, zLevel, pRow );
               
#if 0
                  std::cout << "Rasterizing span" << std::endl;
                  std::cout << "   (x, y, z): (" << span.x << ", " << span.y << ", " << zLevel << ")" << std::endl;
                  std::cout << "   length: " << span.end << std::endl;
                  std::cout << "   r: " << span.red << std::endl;
                  std::cout << "   g: " << span.green << std::endl;
                  std::cout << "   b: " << span.blue << std::endl;
                  std::cout << "   a: " << span.alpha << std::endl;
                    
                  std::cout << "   pRow: " << (void*) pRow << std::endl; 
                  std::cout << "   *pRow: " << static_cast< unsigned int >( *pRow ) << std::endl;
#endif  
                  
                  /** if (!clip_span(ctx, span)) { return; } */
                  
                  const GLuint n = span.end;
                  
                  ChanType  r = span.red;
                  ChanType  g = span.green;
                  ChanType  b = span.blue;
                  ChanType  a = span.alpha;
                  
                  const ChanType  dr = span.redStep;
                  const ChanType  dg = span.greenStep;
                  const ChanType  db = span.blueStep;
                  const ChanType  da = span.alphaStep;
                                    
                  for ( unsigned int i = 0; i < n; i++ )
                    {
                    m_FragmentProcessor( FixedToChan( r ), FixedToChan( g ), FixedToChan( b ), FixedToChan( a ) );
                    r += dr;
                    g += dg;
                    b += db;
                    a += da;
                    }
               
               }

               /*
                * Advance to the next scan line.  Compute the
                * new edge coordinates, and adjust the
                * pixel-center x coordinate so that it stays
                * on or inside the major edge.
                */
               span.y++;
               lines--;

               fxLeftEdge += fdxLeftEdge;
               fxRightEdge += fdxRightEdge;

               fError += fdError;
               if (fError >= 0) {
                  fError -= INTERP_ONE;

                  pRow = pRow + dPRowOuter;
                  
                  rLeft += fdrOuter;
                  gLeft += fdgOuter;
                  bLeft += fdbOuter;
                  aLeft += fdaOuter;
               }
               else {
                  pRow = pRow + dPRowInner;
                  
                  rLeft += fdrInner;
                  gLeft += fdgInner;
                  bLeft += fdbInner;
                  aLeft += fdaInner;
               }
            } /*while lines>0*/

         } /* for subTriangle */

      }
   }
}

  
  
  
  
  
} // end namespace kvl

#endif
