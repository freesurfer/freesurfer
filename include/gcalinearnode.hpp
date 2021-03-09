/**
 * @brief Class to hold a volume of GCA nodes in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
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
 */

#ifndef GCA_LINEAR_NODE_HPP
#define GCA_LINEAR_NODE_HPP

#include <vector>

#include "timer.h"
#include "gca.h"


// Forward declaration
namespace GPU
{
namespace Classes
{
class GCAnodeGPU;
}
}

//! Catch-all namespace
namespace Freesurfer
{

// Forward declaration
class const_GCAnode;
class const_GCAnode_GC1D;

//! Class to hold a 3D volume of GCA nodes in linear memory (ninputs=1)
class GCAlinearNode
{
public:
  GCAlinearNode( void ) : hasGibbsNeighbourhood(true),
			  xDim(0), yDim(0), zDim(0),
			  n4D(0), n5D(0), n6D(0),
			  gc1dNeighbourDim(GIBBS_NEIGHBORHOOD),
			  offsets4D(), offsets5D(), offsets6D(),
			  nodeMaxLabels(),
			  nodeTotalTraining(),
			  nodeLabels(),
			  means(), variances(),
			  nJustPriors(), nTraining(), regularised(),
			  gc1dDirecLabelPriors(),
			  gc1dDirecLabels() {};
  
  // -------------------------------------------------
  bool hasGibbsNeighbourhood;

  // -------------------------------------------------
  // 3D data access

  //! Accessor for nodeMaxLabels
  inline int maxLabels( const int ix, const int iy, const int iz ) const
  {
    const size_t idx = this->index3D(ix,iy,iz);

    return( this->nodeMaxLabels.at(idx) );
  }

  //! Mutator for nodeMaxLabels
  inline int& maxLabels( const int ix, const int iy, const int iz )
  {
    const size_t idx = this->index3D(ix,iy,iz);

    return( this->nodeMaxLabels.at(idx) );
  }

  // --------------

  //! Accessor for nodeTotalTraining
  inline int totalTraining( const int ix,
                            const int iy,
                            const int iz ) const
  {
    const size_t idx = this->index3D(ix,iy,iz);

    return( this->nodeTotalTraining.at(idx) );
  }

  //! Mutator for nodeTotalTraining
  inline int& totalTraining( const int ix,
                             const int iy,
                             const int iz )
  {
    const size_t idx = this->index3D(ix,iy,iz);

    return( this->nodeTotalTraining.at(idx) );
  }

  // -------------------------------------------------


  //! Dynamic computation of gc1dCount
  inline int gc1dCount( const int ix,
                        const int iy,
                        const int iz ) const
  {
    /*!
    This is computed as the difference between
    consecutive entries of the offsets4D array
    */
    const size_t idx3D = this->index3D(ix,iy,iz);

    long long currOffset = this->offsets4D.at( idx3D );
    long long nextOffset = this->offsets4D.at( idx3D+1 );

    return( nextOffset - currOffset );
  }


  // -------------------------------------------------
  // 4D data access

  //! Accessor for nodeLabels
  inline unsigned short labelsAtNode( const int ix,
                                      const int iy,
                                      const int iz,
                                      const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nodeLabels.at(idx) );
  }

  //! Mutator for nodeLabels
  inline unsigned short& labelsAtNode( const int ix,
                                       const int iy,
                                       const int iz,
                                       const int iGC1D )
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nodeLabels.at(idx) );
  }

  // --

  //! Accessor for means
  inline float meansAtNodeGC1D( const int ix,
                                const int iy,
                                const int iz,
                                const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->means.at(idx) );
  }

  //! Mutator for means
  inline float& meansAtNodeGC1D( const int ix,
                                 const int iy,
                                 const int iz,
                                 const int iGC1D )
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->means.at(idx) );
  }

  // --

  //! Accessor for variances
  inline float variancesAtNodeGC1D( const int ix,
                                    const int iy,
                                    const int iz,
                                    const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->variances.at(idx) );
  }

  //! Mutator for variances
  inline float& variancesAtNodeGC1D( const int ix,
                                     const int iy,
                                     const int iz,
                                     const int iGC1D )
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->variances.at(idx) );
  }

  // --

  //! Accessor for nJustPriors
  inline short nJustPriorsAtNodeGC1D( const int ix,
                                      const int iy,
                                      const int iz,
                                      const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nJustPriors.at(idx) );
  }

  //! Mutator for nJustPriors
  inline short& nJustPriorsAtNodeGC1D( const int ix,
                                       const int iy,
                                       const int iz,
                                       const int iGC1D )
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nJustPriors.at(idx) );
  }

  // --

  //! Accessor for nTraining
  inline int nTrainingAtNodeGC1D( const int ix,
                                  const int iy,
                                  const int iz,
                                  const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nTraining.at(idx) );
  }

  //! Mutator for nTraining
  inline int& nTrainingAtNodeGC1D( const int ix,
                                   const int iy,
                                   const int iz,
                                   const int iGC1D )
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nTraining.at(idx) );
  }

  // --

  //! Accessor for regularised
  inline char regularisedAtNodeGC1D( const int ix,
                                     const int iy,
                                     const int iz,
                                     const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->regularised.at(idx) );
  }

  //! Mutator for regularised
  inline char& regularisedAtNodeGC1D( const int ix,
                                      const int iy,
                                      const int iz,
                                      const int iGC1D )
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->regularised.at(idx) );
  }

  // -------------------------------------------------


  //! Dynamic computation of nLabelsAtNodeGC1Ddirection
  inline int nLabelsAtNodeGC1Ddirection( const int ix,
                                         const int iy,
                                         const int iz,
                                         const int iGC1D,
                                         const int iDirec ) const
  {
    /*!
    This is computed as the difference between
    consecutive entries of the offsets6D array
    */
    const size_t idx5D = this->index5D(ix,iy,iz,iGC1D,iDirec);

    long long currOffset = this->offsets6D.at( idx5D );
    long long nextOffset = this->offsets6D.at( idx5D+1 );

    return( nextOffset - currOffset );
  }

  // -------------------------------------------------
  // 6D data access


  //! Accesor for gc1dDirecLabelPriors
  inline float labelPriorsAtNodeGC1Ddirection( const int ix,
      const int iy,
      const int iz,
      const int iGC1D,
      const int iDirec,
      const int iLabel ) const
  {
    const size_t idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

    return( this->gc1dDirecLabelPriors.at(idx) );
  }

  //! Mutator for gc1dDirecLabelPriors
  inline float& labelPriorsAtNodeGC1Ddirection( const int ix,
      const int iy,
      const int iz,
      const int iGC1D,
      const int iDirec,
      const int iLabel )
  {
    const size_t idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

    return( this->gc1dDirecLabelPriors.at(idx) );
  }

  // --

  //! Accessor for gc1dDirecLabels
  inline unsigned short labelsAtNodeGC1Ddirection( const int ix,
      const int iy,
      const int iz,
      const int iGC1D,
      const int iDirec,
      const int iLabel ) const
  {
    const size_t idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

    return( this->gc1dDirecLabels.at(idx) );
  }

  //! Mutator for gc1dDirecLabels
  inline unsigned short& labelsAtNodeGC1Ddirection( const int ix,
      const int iy,
      const int iz,
      const int iGC1D,
      const int iDirec,
      const int iLabel )
  {
    const size_t idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

    return( this->gc1dDirecLabels.at(idx) );
  }

  // -------------------------------------------------

  //! Method to extract node data from a GCA
  void Exhume( const GCA* const src );

  //! Method to place node data into a GCA
  void Inhume( GCA* dst ) const;

  //! Method to remove nodes from a GCA prior to inhumation of new data
  void ScorchNodes( GCA* targ ) const;


  //! Method to print out timers
  void PrintStats( std::ostream& os = std::cout ) const;

  // -------------------------------------------------
private:

  // Dimensions
  long long xDim;
  long long yDim;
  long long zDim;
  size_t n4D;
  size_t n5D;
  size_t n6D;
  const int gc1dNeighbourDim;


  //! Stores offsets of (variable length) 4th dimension
  std::vector<size_t> offsets4D;

  //! Stores offsets of 5th dimension (5th dim itself always of
  // length gc1dNeighbourDim)
  std::vector<size_t> offsets5D;

  //! Stores offsets of (variable length) 6th dimension
  std::vector<size_t> offsets6D;


  //! Stores the value of max_labels for each node (3D)
  std::vector<int> nodeMaxLabels;
  //! Stores value of total_training for each node (3D)
  std::vector<int> nodeTotalTraining;


  //! Stores a list of labels for each node (4D)
  std::vector<unsigned short> nodeLabels;

  //! Stores the mean for each GC1D of each node (4D)
  std::vector<float> means;
  //! Stores the variance for each GC1D of each node (4D)
  std::vector<float> variances;
  //! Stores n_just_priors for each GC1D of each node (4D)
  std::vector<short> nJustPriors;
  //! Stores ntraining for each GC1D of each node (4D)
  std::vector<int> nTraining;
  //! Stores regularized for each GC1D of each node (4D)
  std::vector<char> regularised;


  //! Stores label_priors for each direction of each GC1D of each node (6D)
  std::vector<float> gc1dDirecLabelPriors;
  //! Stores labels for each direction of each GC1D of each node (6D)
  std::vector<unsigned short> gc1dDirecLabels;


  // -------------------------------------------------

  //! Method to extract the max. dimensions from a GCA
  void ExtractDims( const GCA* const src );

  //! Method to allocate memory according to internals
  void Allocate( void );

  //! Create a const_GCAnode
  const_GCAnode GetConstNode( const int ix,
                              const int iy,
                              const int iz ) const;

  // -------------------------------------------------

  //! Exhumation time
  mutable long exhumeTime;
  //! Inhumation time
  mutable long inhumeTime;

  // -------------------------------------------------

  //! Index computation for 3D matrices
  inline size_t index3D( const int ix,
                         const int iy,
                         const int iz ) const
  {
    if ( (ix<0) || (ix>=this->xDim) ||
         (iy<0) || (iy>=this->yDim) ||
         (iz<0) || (iz>=this->zDim ) )
    {
      std::cerr << __FUNCTION__
                << ": Index out of range" << std::endl;
      abort();
    }

    size_t index;
    index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );

    return( index );
  }

  //! Index computation for 4D arrays
  inline size_t index4D( const int ix,
                         const int iy,
                         const int iz,
                         const int iGC1D ) const
  {
    if ( (ix<0) || (ix>=this->xDim) ||
         (iy<0) || (iy>=this->yDim) ||
         (iz<0) || (iz>=this->zDim ) )
    {
      std::cerr << __FUNCTION__
                << ": Index out of range" << std::endl;
      abort();
    }

    const size_t idx3D = this->index3D(ix,iy,iz);

    const long long currOffset = this->offsets4D.at( idx3D );
    const long long nextOffset = this->offsets4D.at( idx3D + 1 );

    if ( (iGC1D<0) || (iGC1D>=(nextOffset-currOffset) ) )
    {
      std::cerr << __FUNCTION__
                << ": iGC1D out of range" << std::endl;
      abort();
    }

    return( currOffset + iGC1D );
  }


  //! Index computation for 5D arrays
  inline size_t index5D( const int ix,
                         const int iy,
                         const int iz,
                         const int iGC1D,
                         const int iDirec ) const
  {
    /*!
    No actual 5D arrays exist in the class, since nlabels in GC1D
    can be dynamically computed.
    However, the Exhume method does have a 5D temporary
    to generate the 6D offset array
    */

    if ( (ix<0) || (ix>=this->xDim) ||
         (iy<0) || (iy>=this->yDim) ||
         (iz<0) || (iz>=this->zDim ) )
    {
      std::cerr << __FUNCTION__
                << ": Index out of range" << std::endl;
      abort();
    }

    const size_t idx3D = this->index3D(ix,iy,iz);

    const long long currGC1dOffset = this->offsets4D.at( idx3D );
    const long long nextGC1dOffset = this->offsets4D.at( idx3D + 1 );

    if ( (iGC1D<0) || (iGC1D>=(nextGC1dOffset-currGC1dOffset) ) )
    {
      std::cerr << __FUNCTION__
                << ": iGC1D out of range" << std::endl;
      abort();
    }

    const long long currOffset = this->offsets5D.at( currGC1dOffset+iGC1D );
    const long long nextOffset = this->offsets5D.at( currGC1dOffset+iGC1D+1 );

    if ( (iDirec<0) || (iDirec>=(nextOffset-currOffset) ) )
    {
      std::cerr << __FUNCTION__
                << ": iDirec out of range" << std::endl;
      abort();
    }

    return( currOffset + iDirec );
  }



  //! Index computation for 6D arrays
  inline size_t index6D( const int ix,
                         const int iy,
                         const int iz,
                         const int iGC1D,
                         const int iDirec,
                         const int iLabel ) const
  {

    if ( (ix<0) || (ix>=this->xDim) ||
         (iy<0) || (iy>=this->yDim) ||
         (iz<0) || (iz>=this->zDim ) )
    {
      std::cerr << __FUNCTION__
                << ": Index out of range" << std::endl;
      abort();
    }

    const size_t idx3D = this->index3D(ix,iy,iz);

    const long long currgc1DOffset = this->offsets4D.at( idx3D );
    const long long nextgc1DOffset = this->offsets4D.at( idx3D + 1 );

    if ( (iGC1D<0) || (iGC1D>=(nextgc1DOffset-currgc1DOffset) ) )
    {
      std::cerr << __FUNCTION__
                << ": iGC1D outof range" << std::endl;
      abort();
    }

    const long long currDirecOffset =
      this->offsets5D.at( currgc1DOffset+iGC1D );
    const long long nextDirecOffset =
      this->offsets5D.at( currgc1DOffset+iGC1D+1 );

    if ( (iDirec<0) || (iDirec>=(nextDirecOffset-currDirecOffset) ) )
    {
      std::cerr << __FUNCTION__
                << ": iDirec out of range" << std::endl;
      abort();
    }

    const long long currOffset =
      this->offsets6D.at( currDirecOffset+iDirec );
    const long long nextOffset =
      this->offsets6D.at( currDirecOffset+iDirec+1 );

    if ( (iLabel<0) || (iLabel>=(nextOffset-currOffset)) )
    {
      std::cerr << __FUNCTION__
                << ": iLabel out of range" << std::endl;
      abort();
    }

    return( currOffset + iLabel );
  }


  friend class const_GCAnode;
  friend class const_GCAnode_GC1D;
  friend class GPU::Classes::GCAnodeGPU;
};


//! Equivalent of a 'const GCA_NODE'
/*!
  This class provides quick addess to a particular 3D
  point of a GCAlinearNode.
  The indexing must match that in GCAlinearNode, or everything
  will fall apart horribly.
  A lot of the precomputations are performed in the constructor
*/

class const_GCAnode
{
public:
  //! Constructor from location
  const_GCAnode( const int _ix,
                 const int _iy,
                 const int _iz,
                 const GCAlinearNode& src ) : gcaln(src),
    idx3d(src.index3D(_ix,_iy,_iz)),
    offset4d(src.offsets4D.at(idx3d)),
    myGC1Dcount(src.offsets4D.at(idx3d+1)-offset4d),
    ix(_ix),
    iy(_iy),
    iz(_iz) {}

  // --------------------

  //! Accessor for max_labels
  inline int maxLabels( void ) const
  {
    return( this->gcaln.nodeMaxLabels.at(this->idx3d) );
  }

  //! Accessor for totalTraining
  inline int totalTraining( void ) const
  {
    return( this->gcaln.nodeTotalTraining.at(this->idx3d) );
  }

  //! Accessor for gc1dCount
  inline int gc1dCount( void ) const
  {
    return( this->myGC1Dcount );
  }

  //! Accessor for nodeLabels
  inline unsigned short labels( const int iGC1D ) const
  {
    if ( (iGC1D<0) || (iGC1D>=this->myGC1Dcount) )
    {
      std::cerr << __FUNCTION__
                << ": Out of range " << iGC1D
                << std::endl;
      abort();
    }

    return( this->gcaln.nodeLabels.at(this->offset4d+iGC1D) );
  }

  // ----------------------

  //! Create the GC1D
  const_GCAnode_GC1D GetConstGC1D( const int iGC1D ) const;

private:
  //! The GCAlinearNode we're part of
  const GCAlinearNode& gcaln;
  //! Precomputed linear index for 3D data
  const size_t idx3d;
  //! Precomputed linear start index for 4D data
  const size_t offset4d;
  //! Length of the 4th dimension
  const int myGC1Dcount;

  //! Original location
  const int ix, iy, iz;

};



//! Equivalent of a 'const GC1D' for a GCA_NODE
class const_GCAnode_GC1D
{
public:
  const_GCAnode_GC1D( const int _ix,
                      const int _iy,
                      const int _iz,
                      const int _iGC1D,
                      const GCAlinearNode& src ) : gcaln(src),
    offset4d(src.index4D(_ix,_iy,_iz,_iGC1D)) {}

  // -----------------------------
  // 4D data

  //! Accessor for the mean
  inline float mean( void ) const
  {
    return( this->gcaln.means.at(this->offset4d) );
  }

  //! Accessor for the variance
  inline float variance( void ) const
  {
    return( this->gcaln.variances.at(this->offset4d) );
  }

  //! Accessor for nJustPriors
  inline short nJustPriors( void ) const
  {
    return( this->gcaln.nJustPriors.at(this->offset4d) );
  }

  //! Accessor for nTraining
  inline int nTraining( void ) const
  {
    return( this->gcaln.nTraining.at(this->offset4d) );
  }

  //! Accessor for regularised
  inline int regularised( void ) const
  {
    return( this->gcaln.regularised.at(this->offset4d ) );
  }

  // -------------------------------
  // 5D data

  //! Accessor for nlabels
  inline short nLabels( const int iDirec ) const
  {
    /*!
    This is basically an alternative for
    GCAlinearNode::nLabelsAtNodeGC1Ddirection
    */
    if ( iDirec >= this->gcaln.gc1dNeighbourDim )
    {
      std::cerr << __FUNCTION__
                << ": iDirec out of range"
                << std::endl;
    }

    const size_t idx5D = this->gcaln.offsets5D.at( this->offset4d ) + iDirec;

    const long long currOffset = this->gcaln.offsets6D.at( idx5D );
    const long long nextOffset = this->gcaln.offsets6D.at( idx5D+1 );

    return( nextOffset - currOffset );
  }

  // -------------------------------
  // 6D data

  //! Accessor for labels
  inline unsigned short labels( const int iDirec, const int iLabel ) const
  {
    /*!
    This is basically an alternative for
    GCAlinearNode::labelsAtNodeGC1Ddirection
    */
    const size_t idx6D = this->index6D(iDirec,iLabel);

    return( this->gcaln.gc1dDirecLabels.at(idx6D) );
  }

  //! Accessor for labelPriors
  inline float labelPriors( const int iDirec, const int iLabel ) const
  {
    /*!
    This is an alternative for
    GCAlinearNode::labelPriorsAtNodeGC1Ddirection
    */
    const size_t idx6D = this->index6D(iDirec,iLabel);

    return( this->gcaln.gc1dDirecLabelPriors.at(idx6D) );
  }


private:
  //! The GCAlinearNode we're part of
  const GCAlinearNode& gcaln;
  //! Precomputed linear start index for 4D data
  const size_t offset4d;


  //! Index computation for 6D data
  inline size_t index6D( const int iDirec,
                         const int iLabel ) const
  {
    if ( (iDirec<0) || (iDirec>= this->gcaln.gc1dNeighbourDim ) )
    {
      std::cerr << __FUNCTION__
                << ": iDirec out of range"
                << std::endl;
    }

    const size_t idx5D = this->gcaln.offsets5D.at( this->offset4d ) + iDirec;

    const long long currOffset = this->gcaln.offsets6D.at( idx5D );
    const long long nextOffset = this->gcaln.offsets6D.at( idx5D+1 );

    if ( (iLabel<0) || (iLabel>=(nextOffset-currOffset)) )
    {
      std::cerr << __FUNCTION__
                << ": iLabel out of range" << std::endl;
      abort();
    }

    return( currOffset + iLabel );
  }

};

}

#endif
