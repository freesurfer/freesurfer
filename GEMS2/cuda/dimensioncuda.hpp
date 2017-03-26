#pragma once

#include <limits>
#include <stdexcept>

namespace kvl {
  namespace cuda {
    template<unsigned char nDims, typename IndexType = size_t>
    class Dimension {
    public:
      Dimension() {
	static_assert(std::numeric_limits<IndexType>::is_integer,
		      "Must have integral IndexType");
	static_assert(!std::numeric_limits<IndexType>::is_signed,
		      "Must have unsigned IndexType");
      }

      bool operator==( const Dimension& d ) const {
	bool res = true;
	
	for( unsigned char i=0; i<nDims; i++ ) {
	  res &= (this->lengths[i] == d.lengths[i]);
	}
	
	return res;
      }
      
      bool operator!=( const Dimension& d ) const {
	return !(this->operator==(d));
      }
      
      IndexType& operator[]( const unsigned char idx ) {
	if( idx>=nDims ) {
	  throw std::range_error("Index out of range");
	}
	
	return this->lengths[idx];
      }

      const IndexType& operator[]( const unsigned char idx ) const {
	if( idx>=nDims ) {
	  throw std::range_error("Index out of range");
	}
	
	return this->lengths[idx];
      }

      // Variadic function to check if a point is in range
      template<typename NextType, typename... Values>
      bool PointInRange(const NextType iVal, const Values... vals) const {
	static_assert((1+sizeof...(Values))==nDims,
		      "Must call with nDims arguments");
	this->CheckForIndexType<NextType>();
	
	IndexType location[nDims];
	
	this->copyToArray(location, iVal, vals...);
	
	return this->PointInRangeFromArray(location);
      }
      
      bool PointInRangeFromArray(const IndexType location[nDims]) const {
	bool result = true;

	for( unsigned char i=0; i<nDims; i++ ) {
	  result = result && (location[i] < this->lengths[i]);
	}

	return result;
      }

      template<typename NextType, typename... Values>
      size_t GetLinearIndex(const NextType iVal, const Values... vals) const {
	static_assert((1+sizeof...(Values))==nDims,
		      "Must call with nDims arguments");
	this->CheckForIndexType<NextType>();
	
	if( !this->PointInRange( iVal, vals... ) ) {
	  throw std::range_error("At least one index out of range");
	}
	
	IndexType location[nDims];
	
	this->copyToArray(location, iVal, vals...);
	
	return this->GetLinearIndexFromArray(location);
      }

      size_t GetLinearIndexFromArray(const IndexType location[nDims]) const {
	size_t result;
      
	result = location[0];
	
	for( unsigned char i=1; i<nDims; i++ ) {
	  result = location[i] + (result * this->lengths[i]);
	}

	return result;
      }


    private:
      IndexType lengths[nDims];


      // Base case of copyToArray
      void copyToArray(IndexType*) const { }

      // Extracts variadic template arguments into a array of length nDims
      template<typename NextType, typename... Values>
      void copyToArray(IndexType* nextLoc, const NextType iVal, const Values... vals) const {
	this->CheckForIndexType<NextType>();
      
	*nextLoc = iVal;
	this->copyToArray(++nextLoc,vals...);
      }
    
      template<typename T>
      void CheckForIndexType() const {
	static_assert(std::is_same<IndexType,T>::value,
		      "Type passed must be IndexType");
      }
    };
  }
}
