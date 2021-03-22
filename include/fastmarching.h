/**
 * @brief fast marching algorithm
 *
 */
/*
 * Original Author: Florent Segonne  
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

// include guard
#ifndef fastmarching_h
#define fastmarching_h

#include "mri.h"
#include "mrisurf.h"
#include "diag.h"

MRI *MRIextractDistanceMap(MRI *mri_src,
                           MRI *mri_dst,
                           int label,
                           float max_distance,
                           int mode,
                           MRI *mri_mask);

void MRISextractOutsideDistanceMap(MRIS *mris,
                                   MRI *mri_src,
                                   int label,
                                   int offset,
                                   float resolution,
                                   float max_distance);

#include <queue>
#include <functional>
#include <climits>
#include <list>



#include "mrisurf.h"
#include "error.h"


#define mapMRI_XYZ(mri,x,y,z) for(int z =0 ; z < mri->depth ; z++) \
                 for(int y = 0 ; y < mri->height ; y++) \
                 for(int x = 0 ; x < mri->width ; x++)

const float DELTA = 1.0f;
const float EPS = 1e-6f;
const float BIG = 1e6f;

class stCoord
{
public:
  int x,y,z;
  stCoord(int _x, int _y, int _z) : x(_x), y(_y), z(_z)
  {}
};

typedef std::list<stCoord> CoordList;
typedef CoordList::iterator CoordListIterator;
typedef CoordList::const_iterator const_CoordListIterator;

typedef enum
{
  eAlive=0, eTrial=1, eFar=2, eForbidden=3
}
eState;

template <int sign = +1>
class FastMarching
{

  // type definition for the min-heap
class HeapCompare : std::binary_function<stCoord, stCoord, bool>
  {
  protected:
    MRI *mri;
  public:
    HeapCompare(MRI *_mri) : mri(_mri)
    {}
    
    bool operator() (const stCoord &a, const stCoord &b) const
    {
      return (sign * 
              MRIFvox(mri,a.x,a.y,a.z) > sign * 
              MRIFvox(mri,b.x,b.y,b.z));
    }
    
  };

  typedef std::priority_queue<stCoord,std::vector<stCoord>,HeapCompare> CoordHeap;

protected:

  //members

public:

  MRI *mri;
  MRI *mri_mask ;
  CoordHeap trial_heap;
  CoordList alive_list;
  MRI *status;
  float limit;
  int width,height,depth;

public:

  //constructor, destructor
  FastMarching(MRI *_mri, MRI *_mri_mask): mri(_mri), mri_mask(_mri_mask), trial_heap( HeapCompare( mri ) )
  {
    status = MRIalloc( mri->width, mri->height, mri->depth, MRI_UCHAR );
    width = mri->width;
    height = mri->height;
    depth = mri->depth;
    Init();
    // set mask value in status volume to eForbidden here 
    if (mri_mask)
    {
      mapMRI_XYZ(status,x,y,z) 
        if ((int)MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
        {
          MRIvox(status,x,y,z)=eForbidden; // Por qua? Por qua pa?
          MRIFvox(mri, x, y, z) = limit ;
        }
    }
  }

  ~FastMarching()
  {
    if (mri_mask)
    {
      mapMRI_XYZ(status,x,y,z) 
        if ((int)MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
        {
          MRIFvox(mri, x, y, z) = limit ;
        }
    }
    MRIfree(&status);
  }

  // Iterators on the list of alive points
  const_CoordListIterator begin() const
  {
    return alive_list.begin();
  }
  const_CoordListIterator end() const
  {
    return alive_list.end();
  }

  void Init()
  {
    while (!trial_heap.empty()) trial_heap.pop();
    alive_list.clear();
    mapMRI_XYZ(status,x,y,z) MRIvox(status,x,y,z)=eFar;
  }
  
  void SetLimit(float _limit)
  {
    limit=_limit;
  }

  void AddTrialPoint(const int x, const int y, const int z)
  {
    if (MRIvox(status,x,y,z) == eFar)
    {
      MRIvox(status,x,y,z) = eTrial;
      MRIFvox(mri,x,y,z)=_UpdateValue(x,y,z);
      trial_heap.push(stCoord(x,y,z));
    }
  }
  
  /**
   * If the status of the point is eFar, then add it to the alive list
   */
  void AddAlivePoint(const int x, const int y, const int z)
  {
    if ( MRIvox(status,x,y,z) == eFar)
    {
      MRIvox(status,x,y,z) = eAlive;
      alive_list.push_back(stCoord(x,y,z));
    }
  }
  
  void InitTrialFromAlive()
  {

    for (const_CoordListIterator i=begin(); i!=end(); i++)
    {
      
      const int x = i->x;
      const int y = i->y;
      const int z = i->z;

      if (x > 0) 
        AddTrialPoint(x-1,y,z);

      if (x<width-1) 
        AddTrialPoint(x+1,y,z);

      if (y>0) 
        AddTrialPoint(x,y-1,z);

      if (y<height-1) 
        AddTrialPoint(x,y+1,z);

      if (z>0) 
        AddTrialPoint(x,y,z-1);

      if (z<depth-1) 
        AddTrialPoint(x,y,z+1);

    }
  }

  void Run(float _limit)
  {
    // I used google translate to translate the French...
    
    // On memorise l'etendue du fast marching
    // One stores the extent of the fast marching
    limit = _limit;

    // Tant qu'il y a des points trial
    // As long as there are points trial
    while (!trial_heap.empty())
    {
      
      // On enleve de la pile le point trial de valeur la plus faible
      // One removes pile the weakest point trial of the value 
      stCoord pt = trial_heap.top();

      // Si on a atteint la limite
      // If the limit were reached
      if ( sign * MRIFvox( mri, pt.x, pt.y, pt.z ) >= sign * limit )
      {
        // On vide la pile et on marque tous ses points comme far
        // One empties the pile and one marks all his points like far
        while (!trial_heap.empty())
        {
          pt = trial_heap.top();
          if (MRIvox(status, pt.x, pt.y, pt.z) == eForbidden)
            DiagBreak() ;
          MRIvox( status, pt.x, pt.y, pt.z ) = eFar;
          MRIFvox( mri, pt.x, pt.y, pt.z ) = limit;
          trial_heap.pop();
        }
        return;
      }
      // On enleve le point de la pile et on l'ajoute a la liste des points alive
      // One removes the point of the pile and one adds it to the alive list of the points
      if (MRIvox(status, pt.x, pt.y, pt.z) == eForbidden)
        DiagBreak() ;
      trial_heap.pop();
      MRIvox(status,pt.x,pt.y,pt.z) = eAlive;
      alive_list.push_back(pt);

      // On determine les voisins du point et on les met a jour
      // Given the neighbors of the point and one they are updated
      const int x = pt.x;
      const int y = pt.y;
      const int z = pt.z;

      if ( x > 0 ) 
        _UpdatePoint( x-1, y, z );        
      if ( x < width-1 ) 
        _UpdatePoint( x+1, y, z );

      if ( y > 0 ) 
        _UpdatePoint( x, y-1, z );
      if ( y < height-1 ) 
        _UpdatePoint( x, y+1, z );

      if ( z > 0 ) 
        _UpdatePoint( x, y, z-1 );
      if ( z < depth-1 ) 
        _UpdatePoint( x, y, z+1 );
    }
  }

private:
  void _UpdatePoint(int x, int y, int z)
  {
    const eState st = (eState)MRIvox(status,x,y,z);
    if (x == Gx && y == Gy && z == Gz)
    DiagBreak() ;
    if (st == eFar)
    {
      MRIvox(status,x,y,z)= eTrial;
      MRIFvox(mri,x,y,z) = _UpdateValue(x,y,z);
      trial_heap.push(stCoord(x,y,z));
    }
    else if (st == eTrial)
    {
      MRIFvox(mri,x,y,z) = _UpdateValue(x,y,z);
    }
  }

  float _GetValue(const int x, const int y, const int z) const
  {
    if (MRIvox(status,x,y,z) == eFar) return limit;
    return MRIFvox(mri,x,y,z);
  }

  float _UpdateValue(const int x, const int y, const int z) const
  {

    // On prend le minimum de chaque paire de voisins
    // One takes the minimum of each pair of neighbors
    float A = (x==0) ? sign*_GetValue(x+1,y,z) : (x==width-1) ? sign*_GetValue(x-1,y,z) : MIN( sign*_GetValue(x+1,y,z), sign*_GetValue(x-1,y,z) );
    float B = (y==0) ? sign*_GetValue(x,y+1,z) : (y==height-1) ? sign*_GetValue(x,y-1,z) : MIN( sign*_GetValue(x,y+1,z), sign*_GetValue(x,y-1,z) );
    float C = (z==0) ? sign*_GetValue(x,y,z+1) : (z==depth-1) ? sign*_GetValue(x,y,z-1) : MIN( sign*_GetValue(x,y,z+1), sign*_GetValue(x,y,z-1) );

    // On reordonne les valeurs pour avoir C>=B>=A
    // One reorders the values to have C>=B>=A
    if (x == Gx && y == Gy && z == Gz)
    DiagBreak() ;
    if (A>B)
    {
      const float tmp = A;
      A = B;
      B = tmp;
    }

    if (B>C)
    {
      const float tmp = B;
      B = C;
      C = tmp;
    }
    if (A>B)
    {
      const float tmp = A;
      A = B;
      B = tmp;
    }

    float returnValue = sign * ( A+1 );

    // On suppose sol>=C : premier trinome
    float a = 3;
    float b = -(A+B+C);
    float c = A*A + B*B + C*C - 1;
    float delta = b*b - a*c;
    float sol = ( -b + ::sqrt(delta) ) / a;    

    // On a bien sol>=C, on a gagne
    if( delta >= 0 && sol+EPS>=C ) {

      returnValue = sign*sol;

    } else {

      // On supppose B<=sol<C : deuxieme trinome
      a = 2;
      b = -(A+B);
      c = A*A + B*B - 1;
      delta = b*b - a*c;
      
      if( delta >= 0 ) {
        
        const float sol = ( -b + ::sqrt(delta) ) / a;
        // On a bien sol>=B, on a gagne
        if( sol + EPS >= B ) {
          returnValue = sign*sol;
        }
      }

    }
        
    // On suppose A<=sol<B
    return returnValue;
  }

public:

  void InitFromMRI( MRI *_mri, const int label)
  {
    
    // sign is the template parameter (+/-1)
    if ( sign > 0 )
    {
      // this is a macro that goes through all elements of _mri.  All pixels 
      // indices in _mri that aren't equal to label will be set the limit
      mapMRI_XYZ(_mri,x,y,z)
      {
        const int val = static_cast< int >( round( MRIgetVoxVal( _mri, x, y, z, 0 ) ) );      
        
        // if the pixel value isn't the label, then set mri to be at the limit
        if( val != label )
        {
          // mri was in the constructor
          MRIFvox(mri,x,y,z) = limit;
        } else {
          MRIvox( status, x, y, z ) = eForbidden;
        }
      }
    }
    else
    {
      // all pixels in mri that correspond to label _mri will be set to the
      // limit
      mapMRI_XYZ(_mri,x,y,z)
      {
        // this the the value in our label volume
        const int voxel = static_cast< int >( round( MRIgetVoxVal( _mri, x, y, z, 0 ) ) );      

        if ( voxel == label ) {
          MRIFvox( mri, x, y, z ) = limit;
        } else {
          MRIvox( status, x, y, z ) = eForbidden;
        }
      }
    }
    
    mapMRI_XYZ(_mri,x,y,z)
    {

      // these are indices and will be either equal to the corresponding
      // index or the index + 1 dependong on if we've iterated beyond the
      // dimensions of our image
      const int px = (x < width-1) ? x+1:x;
      const int py = (y < height-1) ? y+1:y;
      const int pz = (z < depth-1) ? z+1:z;
      
      // get the values of the volume that has the label
      const int val1 = static_cast< int >( round( MRIgetVoxVal( _mri, x, y, z, 0 ) ) );

      // this was the old way of getting the label.  It wasn't returning the
      // right results, so I changed it to the current
      int val2 = static_cast< int >( round( MRIgetVoxVal( _mri, px, y, z, 0 ) ) ); 

      bool add = false;

      // if the values aren't equal and one of them is the label, then
      // run _AddAlivePoint

      if (val1 != val2 && ( val1 == label || val2 == label ) )
      {
        add = true;
        _AddAlivePoint( px, y, z );
      }

      val2 = static_cast< int >( round( MRIgetVoxVal( _mri, x, py, z, 0 ) ) ); 

      if ( val1!=val2 && (val1==label || val2==label))
      {
        add = true;
        _AddAlivePoint(x,py,z);
      }

      val2 = static_cast< int >( round( MRIgetVoxVal( _mri, x, y, pz, 0 ) ) ); 

      if ( val1 != val2 && ( val1 == label || val2 == label ) )
      {
        add = true;
        _AddAlivePoint(x,y,pz);
      }

      if ( add )
        _AddAlivePoint(x,y,z);
        
    }

    InitTrialFromAlive();

  }

  /* volume floats */
#define xVOL(mri,x) (mri->xsize*(x-mri->xstart))
#define yVOL(mri,y) (mri->ysize*(y-mri->ystart))
#define zVOL(mri,z) (mri->zsize*(z-mri->zstart))
  /* volume integers */
#define iVOL(mri,x) ((int)(xVOL(mri,x)+0.5))
#define jVOL(mri,y) ((int)(yVOL(mri,y)+0.5))
#define kVOL(mri,z) ((int)(zVOL(mri,z)+0.5))

#define xSURF(mri,x) (mri->xstart+(float)x/mri->xsize)
#define ySURF(mri,y) (mri->ystart+(float)y/mri->ysize)
#define zSURF(mri,z) (mri->zstart+(float)z/mri->zsize)

  void InitForOutsideMatch(MRI *mri_distance, MRI* _mri,int label)
  {
    double xv,yv,zv;
    MRI *mri_seg=MRIalloc(width,height,depth,MRI_UCHAR);
    mapMRI_XYZ(mri,x,y,z)
    {
      //find coordinates for _mri: mri->surf->_mri
      MRIsurfaceRASToVoxel(_mri,xSURF(mri_distance,x),ySURF(mri_distance,y),zSURF(mri_distance,z),&xv,&yv,&zv);
      //find value
      int val=MRIvox(_mri,(int)xv,(int)yv,(int)zv);
      if (val!=label)
      {
        MRIvox(mri_seg,x,y,z) = 0;
        MRIFvox(mri,x,y,z) = limit;
        continue;
      }
      MRIvox(mri_seg,x,y,z)=1;
      MRIvox(status,x,y,z)=eForbidden;
    }

    mapMRI_XYZ(mri,x,y,z)
    {
      if (MRIFvox(mri_distance,x,y,z)<=0.0f)
      {
        MRIFvox(mri,x,y,z)=limit;
        MRIvox(status,x,y,z)=eForbidden;
      }
    }

    mapMRI_XYZ(mri,x,y,z)
    {
      int val1=MRIvox(mri_seg,x,y,z),val2;
      int px,py,pz;
      px = (x < width-1) ? x+1:x;
      py = (y < height-1) ? y+1:y;
      pz = (z < depth-1) ? z+1:z;
      bool add = false;
      val2=MRIvox(mri_seg,px,y,z);
      if (val1!=val2 )
      {
        add = true;
        _AddAlivePoint(px,y,z);
      }
      val2=MRIvox(mri_seg,x,py,z);
      if (val1!=val2 )
      {
        add = true;
        _AddAlivePoint(x,py,z);
      }
      val2=MRIvox(mri_seg,x,y,pz);
      if (val1!=val2 )
      {
        add = true;
        _AddAlivePoint(x,y,pz);
      }
      if (add) _AddAlivePoint(x,y,z);
    }
    MRIfree(&mri_seg);
    InitTrialFromAlive();

  }

  /**
   * If the state is not eForbidden, add it to the heap.
   */
  void _AddAlivePoint(int x, int y, int z)
  {

    const eState st = (eState)MRIvox(status,x,y,z);

    if (st != eForbidden) {
      MRIFvox(mri,x,y,z) = sign * 0.5;
      AddAlivePoint(x,y,z);
    }
    
  }

};

#endif // ifndef fastmarching_h
