#ifndef _SBL_POINTER_H_
#define _SBL_POINTER_H_
namespace sbl {


/*! \file Pointer.h
    \brief The Pointer module provides a std::auto_ptr clone (based directly
    on the std::auto_ptr code).  We define our own stand-alone version so as 
    to minimize header pollution.
*/


// forward declaration used by aptr_ref
template<class T> class aptr;


/// The aptr_ref class is used to assist copying aptr instances.  It is a direct clone of std::auto_ptr_ref; we define our own so as not to include all of the other definitions.
template<class T> struct aptr_ref {
    explicit aptr_ref( T *p ) : ptr( p ) { }
    T *ptr;
};


/// The aptr class represents a simple smart pointer.  It is a direct clone of std::auto_ptr; we define our own so as not to include all of the other definitions.
template<class T> class aptr {
public:

    /// basic constructor (take ownership)
    explicit aptr( T *p = 0 ) : m_ptr( p ) {}

    /// copy constructor (take ownership)
    aptr( aptr<T> &ap ) : m_ptr( ap.release() ) {}

    /// proxy copy constructor (take ownership)
    aptr( aptr_ref<T> apr ) { 
        T *p = apr.ptr;
        apr.ptr = 0;
        m_ptr = p;
    }

    /// pointer conversion (give up ownership)
    template<class T2> operator aptr<T2>() {
        return aptr<T2>( *this );
    }

    /// pointer proxy conversion (give up ownership)
    template<class T2> operator aptr_ref<T2>() {
        T2 *p = m_ptr;
        aptr_ref<T2> apr( p );
        m_ptr = 0;
        return apr;
    }

    /// pointer conversion assignment operator (take ownership)
    template<class T2> aptr<T> &operator=( aptr<T2> &p ) {
        reset( p.release() );
        return *this;
    }

    /// pointer conversion constructor (take ownership)
    template<class T2> aptr( aptr<T2> &p ) : m_ptr( p.release() ) {}

    /// assignment operator (take ownership)
    aptr<T> &operator=( aptr<T> &ap ) {
        reset( ap.release() );
        return *this;
    }

    /// proxy assignment operator (take ownership)
    aptr<T> &operator=( aptr_ref<T> apr ) {
        T *p = apr.ptr;
        apr.ptr = 0;
        reset( p );
        return *this;
    }

    /// deallocate our object
    ~aptr() { delete m_ptr; }

    /// access the object
    T &operator*() const { return *m_ptr; }

    /// access the object
    T *operator->() const {    return m_ptr; }

    /// access the object
    T *get() const { return m_ptr; }

    /// give up ownership
    T *release() { 
        T *p = m_ptr; 
        m_ptr = 0; 
        return p; 
    }

    /// point to new object
    void reset( T *p = 0 ) { 
        if (p != m_ptr) 
            delete m_ptr; 
        m_ptr = p; 
    }

private:

    // the object to which we are pointing
    T *m_ptr;
};


} // end namespace sbl
#endif // _SBL_POINTER_H_

