/**
 * @file  fixedmap.hpp
 * @brief Hold a fixed size map for the GPU
 *
 * This is developed for VoxInLabelWithPartialVolume, but will be
 * more generally useful. It's a map (in the pattern of std::map),
 * but for a fixed map size
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.2 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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


#include "cudacheck.h"
#include "cudatypeutils.hpp"

#ifndef FIXED_MAP_CUDA_HPP
#define FIXED_MAP_CUDA_HPP


namespace GPU
{

namespace Classes
{

//! Class to hold a key-value pair
template<typename T,typename U>
class KVpair
{
public:
  typedef T key_type;
  typedef U value_type;

  key_type key;
  value_type value;

  __device__ __host__
  KVpair( const key_type _key,
          const value_type _value ) : key(_key),
    value(_value) {};

};

// -------------------------------------------------------

//! Device class for a map of fixed maximum size
template<typename T,typename U,unsigned int maxPairs>
class FixedMap
{
public:
  typedef T key_type;
  typedef U value_type;
  typedef unsigned int index_t;

  typedef KVpair<T,U> KVtype;

  //! Default constructor (can supply default value)
  __device__
  FixedMap( const key_type _empty,
            const value_type _val ) : currPairs(0),
    keys(),
    values(),
    defaultValue(_val),
    emptyKey(_empty)
  {
    for( index_t i=0; i<maxPairs; i++ )
    {
      this->values[i] = this->defaultValue;
      this->keys[i] = this->emptyKey;
    }
  }


  //! Method to add a key to the map
  __device__
  index_t AddKey( const key_type myKey )
  {
    // Check to space
    if( maxPairs == this->currPairs )
    {
      float *p = NULL;
      *p = 0;
    }

    // Locate the insertion location
    index_t insLoc = this->KeyIndexLB( myKey );

    if( insLoc == maxPairs )
    {

      if( this->currPairs == 0 )
      {
        insLoc = 0;
      }
      else if( this->keys[0] > myKey )
      {
        insLoc = 0;
      }
      else if( this->keys[this->currPairs-1] < myKey )
      {
        insLoc = this->currPairs;
      }
    }

    if( this->keys[insLoc] == myKey )
    {
      return( insLoc );
    }

    // Move other data
    for( index_t i=this->currPairs; i>insLoc; i-- )
    {
      this->keys[i] = this->keys[i-1];
      this->values[i] = this->values[i-1];
    }


    // Add the key
    this->keys[insLoc] = myKey;

    // Initialise the value
    this->values[insLoc] = this->defaultValue;

    // Increment number stored
    this->currPairs++;

    return( insLoc );
  }



  //! Method to delete a key from the map
  __device__
  void DeleteKey( const key_type myKey )
  {
    // Locate the key
    index_t delLoc = this->KeyIndex( myKey );

    if( delLoc == maxPairs )
    {
      float *p = NULL;
      *p = 10;
    }

    // Move elements down
    for( index_t i=delLoc; i<this->currPairs-1; i++ )
    {
      this->keys[i] = this->keys[i+1];
      this->values[i] = this->values[i+1];
    }

    this->currPairs--;

    // Put in last element
    this->keys[this->currPairs] = this->emptyKey;
    this->values[this->currPairs] = this->defaultValue;
  }




  //! Method to set a key's value
  __device__
  void SetValue( const key_type key, const U val )
  {
    index_t valLoc = this->KeyIndex( key );

    if( valLoc == maxPairs )
    {
      float *p = NULL;
      *p = 10;
    }

    this->values[valLoc] = val;
  }


  //! Access a value by key (with insertion)
  __device__
  value_type& operator[]( const key_type key )
  {
    index_t valLoc = this->KeyIndex( key );

    if( valLoc == maxPairs )
    {
      valLoc = this->AddKey( key );
    }

    return( this->values[valLoc] );
  }

  //! Access a value by key (const)
  __device__
  value_type& operator[]( const key_type key ) const
  {
    index_t valLoc = this->KeyIndex( key );

    if( valLoc == maxPairs )
    {
      float *p = NULL;
      *p = 10;
    }

    return( this->values[valLoc] );
  }

  //! Access a value by index
  __device__
  value_type& ValueByIndex( const index_t idx )
  {
    if( idx >= this->currPairs )
    {
      float *p = NULL;
      *p = 10;
    }

    return( this->values[idx] );
  }


  //! Method to get the number of keys
  __device__
  index_t size( void ) const
  {
    return( this->currPairs );
  }

  //! Method to clear a map
  __device__
  void clear( void )
  {
    for( index_t i=0; i<this->currPairs; i++ )
    {
      this->values[i] = this->defaultValue;
      this->keys[i] = this->emptyKey;
    }
    this->currPairs = 0;
  }


  //! Access by internal index
  __device__
  KVtype AccessByIndex( const index_t i ) const
  {
    if( i > this->currPairs )
    {
      float *p = NULL;
      *p = 10;
    }

    return( KVtype( this->keys[i], this->values[i] ) );
  }

private:
  //! The current number of key-value pairs stored
  index_t currPairs;
  //! Array of the keys
  key_type keys[maxPairs];
  //! Array of the values
  value_type values[maxPairs];

  //! Default value to use
  const value_type defaultValue;

  //! Value to indicate empty key
  const key_type emptyKey;

  //! Method to get the location of a key
  __device__
  index_t KeyIndex( const int key ) const
  {

    index_t keyLoc;

    if( this->currPairs == 0 )
    {
      keyLoc = maxPairs;
    }
    else if( this->keys[0] > key )
    {
      // Check bottom of range
      keyLoc = maxPairs;
    }
    else if( this->keys[this->currPairs-1] < key )
    {
      // Check top of range
      keyLoc = maxPairs;
    }
    else if( this->keys[0] == key )
    {
      // Check equality at start
      keyLoc = 0;
    }
    else if( this->keys[this->currPairs-1] == key )
    {
      // Check equality at end
      keyLoc = this->currPairs-1;
    }
    else
    {
      // Do bisection search
      index_t lb = 0;
      index_t ub = this->currPairs;

      keyLoc = maxPairs;

      while( ub != (lb+1) )
      {
        index_t mp = (lb+ub)/2;
        if( this->keys[mp] == key )
        {
          keyLoc = mp;
          break;
        }
        else if( this->keys[mp] < key )
        {
          lb = mp;
        }
        else
        {
          ub = mp;
        }
      }

    }

    return( keyLoc );
  }


  //! Get the first index greater than given key
  __device__
  index_t KeyIndexLB( const int key ) const
  {

    index_t keyLoc;

    if( this->currPairs == 0 )
    {
      keyLoc = maxPairs;
    }
    else if( this->keys[0] > key )
    {
      // Check bottom of range
      keyLoc = maxPairs;
    }
    else if( this->keys[this->currPairs-1] < key )
    {
      // Check top of range
      keyLoc = maxPairs;
    }
    else if( this->keys[0] == key )
    {
      // Check equality at start
      keyLoc = 0;
    }
    else if( this->keys[this->currPairs-1] == key )
    {
      // Check equality at end
      keyLoc = this->currPairs-1;
    }
    else
    {
      // Do bisection search
      index_t lb = 0;
      index_t ub = this->currPairs;

      keyLoc = lb;

      do
      {
        index_t mp = (lb+ub)/2;
        if( this->keys[mp] == key )
        {
          keyLoc = mp;
          break;
        }
        else if( this->keys[mp] < key )
        {
          lb = mp;
        }
        else
        {
          ub = mp;
        }
        keyLoc = ub;
      }
      while( ub != (lb+1) );

    }

    return( keyLoc );
  }


  //! Can't default construct
  FixedMap( void ) : currPairs(),
    defaultValue(),
    emptyKey()
  {
    float *p = NULL;
    *p = 10;
  }

};


}
}


#endif
