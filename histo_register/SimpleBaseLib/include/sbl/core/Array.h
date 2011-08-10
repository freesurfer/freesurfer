#ifndef _SBL_ARRAY_H_
#define _SBL_ARRAY_H_
#include <sbl/core/Display.h>
namespace sbl {


// make sure NULL is defined
#ifndef NULL
#define NULL 0
#endif


// The Array class and PtrArray class each holds a vector pointers.
// The Array class deallocates its contents, whereas the PtrArray class does not.
// As such, the PtrArray can store const types (whereas the Array class cannot).
// The PtrArray class is faster and simpler.


// fix(clean): unify this code; lots of duplication between Array and PtrArray


/// The ArrayItem class is a reference-counted pointer wrapper used for making shallow copies of Array objects.
template<typename T> class ArrayItem {
public:

    /// store the pointer
    explicit ArrayItem( T *obj ) {
        assertDebug( obj != NULL );
        m_refCount = 1;
        m_object = obj;
    }

    /// the pointer object should not be deleted until there are no references to it
    ~ArrayItem() {
        assertDebug( m_refCount == 0 );
        assertDebug( m_object == NULL );
    }

    /// increment reference count
    inline void incRef() { m_refCount++; }

    /// decrement reference count
    inline void decRef() { 
        assertDebug( m_refCount );
        m_refCount--; 
        if (m_refCount == 0) {
            delete m_object;
            m_object = NULL;
        }
    }

    /// the number of references to this item
    inline int refCount() const { return m_refCount; }

    /// the item being referenced
    inline T &ref() { return *m_object; }
    inline const T &ref() const { return *m_object; }

private:

    // the number of references to this item
    int m_refCount;

    // the item being referenced
    T *m_object;

    // disable copy constructor and assignment operator
    ArrayItem( const ArrayItem &x );
    ArrayItem &operator=( const ArrayItem &x );
};


//-------------------------------------------
// TEMPLATE-BASED DYNAMIC ARRAY CLASS
//-------------------------------------------


/// The Array class holds a dynamically-resizable array of pointers to objects.  
/// It will deallocate the objects when it is deallocated.
template<typename T> class Array {
public:

    /// create empty array
    Array() { m_count = 0; m_allocCount = 0; m_set = NULL; }

    /// shallow copy constructor
    Array( const Array<T> &arr );

    /// dealloc array elements (if no other array referencing them)
    ~Array();

    /// assignment operator
    Array<T> &operator=( const Array<T> &arr );

    /// returns number of elements in array (some may be NULL if not set using set())
    inline int size() const { return m_count; } 
    inline int count() const { return m_count; }

    /// returns memory used by Array, but *NOT* by objects in array
    // fix(later): handle ArrayItem
    // fix(later): include memory used by objects (via obj.memUsed()
    inline int memUsed() const { return sizeof(Array<T>) + m_allocCount * sizeof(int *); }

    /// reference an item in the set
    inline T &ref( int index ) { return operator[]( index ); }
    inline const T &ref( int index ) const { return operator[]( index ); }
    inline T &operator[]( int index ) { 
        assertDebug( 0 <= index && index < m_count );
        assertDebug( m_set[ index ] );
        return m_set[ index ]->ref();
    }
    inline const T &operator[]( int index ) const { 
        assertDebug( 0 <= index && index < m_count );
        assertDebug( m_set[ index ] );
        return m_set[ index ]->ref();
    }

    /// get index of matching object (according to == operator), or -1 if not found
    int find( const T &obj ) const;

    /// set an item in the set (deallocates / replaces old item, if any)
    /// takes ownership of pointer (do not externally deallocate)
    void set( int index, T *obj );

    /// appends an item
    /// takes ownership of pointer (do not externally deallocate)
    void append( T *obj ) { extend( m_count + 1 ); set( m_count - 1, obj ); }
    void appendCopy( const T &obj ) { append( new T( obj )); }

    /// delete the specified item and slide all subsequent items toward start
    void remove( int index );

    /// remove all items and reset to empty array (as if destruct then construct)
    void reset();

private:

    /// extend the set (keeping old items, if any)
    void extend( int newLength );

    // returns number of elements in array (some may be NULL if not set using set())
    int m_count;

    // the allocated length of the array (>= m_count)
    int m_allocCount;

    // the elements (some may be NULL)
    ArrayItem<T> **m_set;
};


/// shallow copy constructor
template <typename T> Array<T>::Array( const Array<T> &arr ) {
    m_count = 0;
    m_allocCount = 0;
    m_set = NULL;
    extend( arr.m_count );
    for (int i = 0; i < m_count; i++) {
        m_set[ i ] = arr.m_set[ i ];
        m_set[ i ]->incRef();
    }
}


/// deallocate object set and all objects within the array
template <typename T> Array<T>::~Array() {
    if (m_set) {
        for (int i = 0; i < m_count; i++)
            if (m_set[ i ]) {
                m_set[ i ]->decRef();
                if (m_set[ i ]->refCount() == 0)
                    delete m_set[ i ];
            }
        delete [] m_set;
    }
}


/// assignement operator 
template <typename T> Array<T> &Array<T>::operator=( const Array<T> &arr ) {
    reset();
    extend( arr.m_count );
    for (int i = 0; i < m_count; i++) {
        m_set[ i ] = arr.m_set[ i ];
        m_set[ i ]->incRef();
    }
    return *this;
}


/// resize the array (keeping old items, if any)
template <typename T> void Array<T>::extend( int newLength ) {

    // at the end of this function: every element from 0 to m_allocCount - 1 should be valid or NULL

    // check params
    assertDebug( newLength >= m_count );

    // if old items    
    if (m_count) {

        // check whether to create new array
        if (newLength > m_allocCount) {

            // create new array
            m_allocCount *= 2;
            if (m_allocCount < newLength)
                m_allocCount = newLength;
            ArrayItem<T> **newSet = new ArrayItem<T>*[ m_allocCount ];
            assertDebug( newSet );

            // copy old items
            for (int i = 0; i < m_count; i++)
                newSet[ i ] = m_set[ i ];

            // clear allocated space
            for (int i = m_count; i < m_allocCount; i++)
                newSet[ i ] = NULL;

            // copy set
            delete [] m_set;
            m_set = newSet;
        }
        m_count = newLength;

    // length was zero, just alloc space
    } else {
        if (m_set) // could have been allocated zero-length array
            delete [] m_set;
        m_set = new ArrayItem<T>*[ newLength ];
        for (int i = 0; i < newLength; i++)
            m_set[ i ] = NULL;
        m_count = newLength;
        m_allocCount = newLength;
    }
}


/// get index of matching object (according to == operator), or -1 if not found
template <typename T> int Array<T>::find( const T &obj ) const {
    for (int i = 0; i < m_count; i++) 
        if (m_set[ i ]->ref() == obj)
            return i;
    return -1;
}


/// set an item in the set (deallocates / replaces old item, if any)
/// takes ownership of pointer (do not externally deallocate)
template <typename T> void Array<T>::set( int index, T *obj ) {
    assertDebug( 0 <= index && index < m_count );
    assertDebug( obj );
    if (m_set[ index ]) {
        m_set[ index ]->decRef();
        if (m_set[ index ]->refCount() == 0)
            delete m_set[ index ];
    }
    m_set[ index ] = new ArrayItem<T>( obj );
}


/// delete the specified item and slide all subsequent items toward start
template <typename T> void Array<T>::remove( int index ) {
    assertDebug( index >= 0 && index < m_count );

    // remove the item
    m_set[ index ]->decRef();
    if (m_set[ index ]->refCount() == 0)
        delete m_set[ index ];

    // slide all subsequent items toward the start
    for (int i = index; i < m_count - 1; i++) 
        m_set[ i ] = m_set[ i + 1 ];
    m_count--;
    m_set[ m_count ] = NULL;
}


/// remove all items and reset to empty array (as if destruct then construct)
template <typename T> void Array<T>::reset() {

    // deallocate
    if (m_set) {
        for (int i = 0; i < m_count; i++)
            if (m_set[ i ]) {
                m_set[ i ]->decRef();
                if (m_set[ i ]->refCount() == 0)
                    delete m_set[ i ];
            }
        delete [] m_set;
    }

    // re-initialize
    m_count = 0;
    m_allocCount = 0;
    m_set = NULL;
}


//-------------------------------------------
// TEMPLATE-BASED DYNAMIC ARRAY CLASS
//-------------------------------------------


/// The PtrArray class holds a dynamically-resizable array of pointers to objects.  
/// It will not deallocate the objects when it is deallocated.
template<typename T> class PtrArray {
public:

    /// create empty array
    PtrArray() { m_count = 0; m_allocCount = 0; m_set = NULL; }

    /// shallow copy constructor
    PtrArray( const PtrArray<T> &arr );
    
    // basic destructor
    ~PtrArray();

    /// returns number of elements in set (some may be NULL if not set using set())
    inline int size() const { return m_count; }
    inline int count() const { return m_count; }

    /// returns memory used by PtrArray, but *NOT* by objects in set
    inline int memUsed() const { return sizeof(Array<T>) + m_allocCount * sizeof(int *); }

    /// reference an item in the set
    inline T &ref( int index ) { return operator[]( index ); }
    inline const T &ref( int index ) const { return operator[]( index ); }
    inline T &operator[]( int index ) { 
        assertDebug( 0 <= index && index < m_count );
        assertDebug( m_set[ index ] );
        return *m_set[ index ];
    }
    inline const T &operator[]( int index ) const { 
        assertDebug( 0 <= index && index < m_count );
        assertDebug( m_set[ index ] );
        return *m_set[ index ];
    }

    /// get index of matching object (according to == operator), or -1 if not found
    int find( T &obj ) const;

    /// set an item in the set (replaces old item, if any); does not take ownership of pointer
    void set( int index, T *obj );

    /// appends an item; does takes ownership of pointer 
    void append( T *obj ) { extend( m_count + 1 ); set( m_count - 1, obj ); }

    /// remove the specified item and slide all subsequent items toward start
    void remove( int index );

    /// remove all items and reset to empty array (as if destruct then construct)
    void reset();

private:

    /// extend the set (keeping old items, if any)
    void extend( int newLength );

    // the data
    int m_count;
    int m_allocCount;
    T **m_set; 

    // disable assignment operator
    PtrArray<T> &operator=( const PtrArray<T> &arr );
};


/// copy constructor (shallow; does not copy objects, just pointers)
template <typename T> PtrArray<T>::PtrArray( const PtrArray<T> &arr ) {
    m_count = 0;
    m_allocCount = 0;
    m_set = NULL;
    extend( arr.m_count );
    for (int i = 0; i < m_count; i++) 
        m_set[ i ] = arr.m_set[ i ];
}


/// deallocate array
template <typename T> PtrArray<T>::~PtrArray() {
    if (m_set) 
        delete [] m_set;
}


/// resize the array (keeping old items, if any)
template <typename T> void PtrArray<T>::extend( int newLength ) {

    // at the end of this function: every element from 0 to m_allocCount - 1 should be valid or NULL

    // check params
    assertDebug( newLength >= m_count );

    // if old items    
    if (m_count) {

        // check whether to create new array
        if (newLength > m_allocCount) {

            // create new array
            m_allocCount *= 2;
            if (m_allocCount < newLength)
                m_allocCount = newLength;
            T **newSet = new T*[ m_allocCount ];
            assertDebug( newSet );

            // copy old items
            for (int i = 0; i < m_count; i++)
                newSet[ i ] = m_set[ i ];

            // clear allocated space
            for (int i = m_count; i < m_allocCount; i++)
                newSet[ i ] = NULL;

            // copy set
            delete [] m_set;
            m_set = newSet;
        }
        m_count = newLength;

    // length was zero, just alloc space
    } else {
        if (m_set) // could have been allocated zero-length array
            delete [] m_set;
        m_set = new T*[ newLength ];
        for (int i = 0; i < newLength; i++)
            m_set[ i ] = NULL;
        m_count = newLength;
        m_allocCount = newLength;
    }
}


/// get index of matching object (according to == operator), or -1 if not found
template <typename T> int PtrArray<T>::find( T &obj ) const {
    for (int i = 0; i < m_count; i++) 
        if (m_set[ i ] == obj)
            return i;
    return -1;
}


/// set an item in the set (deallocates / replaces old item, if any)
/// takes ownership of pointer (do not externally deallocate)
template <typename T> void PtrArray<T>::set( int index, T *obj ) {
    assertDebug( 0 <= index && index < m_count );
    assertDebug( obj );
    m_set[ index ] = obj;
}


/// delete the specified item and slide all subsequent items toward start
template <typename T> void PtrArray<T>::remove( int index ) {
    assertDebug( index >= 0 && index < m_count );

    // slide all subsequent items toward the start
    for (int i = index; i < m_count - 1; i++) 
        m_set[ i ] = m_set[ i + 1 ];
    m_count--;
    m_set[ m_count ] = NULL;
}


/// remove all items and reset to empty array (as if destruct then construct)
template <typename T> void PtrArray<T>::reset() {

    // deallocate
    if (m_set) 
        delete [] m_set;

    // re-initialize
    m_count = 0;
    m_allocCount = 0;
    m_set = NULL;
}


} // end namespace sbl
#endif // _SBL_ARRAY_H_

