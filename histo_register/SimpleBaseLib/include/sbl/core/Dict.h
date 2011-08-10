#ifndef _SBL_DICT_H_
#define _SBL_DICT_H_
#include <sbl/core/Array.h>
#include <sbl/core/StringUtil.h>
namespace sbl {


// fix(clean): make key type a template parameter (so just one dictionary class with string keys, integer keys, etc.)


//-------------------------------------------
// INTEGER-KEY TEMPLATED DICTIONARY CLASS
//-------------------------------------------


/// A DictItem is a key-value pair in an integer-keyed dictionary.
template<typename T> class DictItem {
public:
    T *item;
    int key;
    DictItem<T> *next;
    inline DictItem( int newKey, T *newItem, DictItem<T> *newNext ) { key = newKey; item = newItem; next = newNext; }
};


/// The Dict class stores objects index by integer keys.
/// The Dict class uses a hash table for fast look-up.
/// (The class uses a power-of-2 number of bins, so certain key patterns could have many collisions, resulting in slow performance.)
template<typename T> class Dict {
public:

    // basic constructor/destructor
    Dict();
    ~Dict();

    /// get item with given key; returns NULL if not found;
    /// does not remove or relinquish control of the object
    T *find( int key ); 
    const T *find( int key ) const;

    /// access as array
    inline int refArrayKey( int index ) const { return m_keyArray[ index ]; }
    inline const T &refArray( int index ) const { return m_array[ index ]; }
    inline T &refArray( int index ) { return m_array[ index ]; }
    inline int count() const { return m_array.count(); }
    inline const Array<T> &arrayRef() const { return m_array; }

    /// add an item to the dictionary; takes ownership of pointer unless dealloc disabled;
    /// replaces item with same key if already exists
    void add( int key, T *obj );

    /// disables object deallocation; the dictionary no longer takes ownership of object pointers
    void disableObjectDealloc() { m_deallocObjs = false; }

    /// remove all items (as if call destructor then constructor)
    void reset();

private:

    // if true, dictionary owns the objects and should delete them on destruction
    bool m_deallocObjs;

    // the number of items in the dictionary
    int m_count;

    // each bin contains a list of items
    DictItem<T> **m_bin;

    // the number of bins
    int m_binCount;

    // a bit mask used to obtain a bin index from an item ID
    int m_binMask;

    // an array of items in the dictionary (for iterating over the members of a dictionary)
    PtrArray<T> m_array; // use non-owning array because dict will own the items
    Array<int> m_keyArray;

    // disable copy constructor and assignment operator
    Dict( const Dict &x );
    Dict &operator=( const Dict &x );
};


// create an empty dictionary
template <typename T> Dict<T>::Dict() {
    m_deallocObjs = true;
    m_binMask = 0xffff;
    m_binCount = m_binMask + 1;
    m_count = 0;
    m_bin = new DictItem<T>*[ m_binCount ];
    assertDebug( m_bin );
    for (int i = 0; i < m_binCount; i++) 
        m_bin[ i ] = NULL;
}


/// delete the items
template <typename T> Dict<T>::~Dict() {
    for (int i = 0; i < m_binCount; i++) {
        DictItem<T> *bin = m_bin[ i ];
        while (bin) {
            DictItem<T>    *delItem = bin;
            bin = bin->next;
            if (m_deallocObjs)
                delete delItem->item;
            delete delItem;
        }
    }
    delete [] m_bin;
}


/// get item with given key; returns NULL if not found;
/// does not remove or relinquish control of the object
template <typename T> T *Dict<T>::find( int key ) {
    DictItem<T> *bin = m_bin[ key & m_binMask ];
    while (bin) {
        if (bin->key == key)
            return bin->item;
        bin = bin->next;
    }
    return NULL;
}


/// get item with given key; returns NULL if not found;
/// does not remove or relinquish control of the object
template <typename T> const T *Dict<T>::find( int key ) const {
    DictItem<T> *bin = m_bin[ key & m_binMask ];
    while (bin) {
        if (bin->key == key)
            return bin->item;
        bin = bin->next;
    }
    return NULL;
}


/// add an item to the dictionary; takes ownership of pointer unless dealloc disabled;
/// replaces item with same key if already exists
template <typename T> void Dict<T>::add( int key, T *obj ) {
    assertDebug( obj );

    // look for existing item
    int binIndex = key & m_binMask;
    DictItem<T> *bin = m_bin[ binIndex ];
    while (bin) {
        if (bin->key == key) {
            assertAlways( true ); // fix(later): does not update array
            if (m_deallocObjs)
                delete bin->item;
            bin->item = obj; 
            return;
        }
        bin = bin->next;
    }
    
    // create new item
    m_bin[ binIndex ] = new DictItem<T>( key, obj, m_bin[ binIndex ] );
    assertDebug( m_bin[ binIndex ] );
    m_count++;

    // add to array
    m_array.append( obj );
    m_keyArray.appendCopy( key );
}


/// remove all items (as if call destructor then constructor)
template <typename T> void Dict<T>::reset() {

    // delete everything
    for (int i = 0; i < m_binCount; i++) {
        DictItem<T> *bin = m_bin[ i ];
        while (bin) {
            DictItem<T>    *delItem = bin;
            bin = bin->next;
            if (m_deallocObjs)
                delete delItem->item;
            delete delItem;
        }
    }
    delete [] m_bin;

    // re-init
    m_deallocObjs = true;
    m_binMask = 0xffff;
    m_binCount = m_binMask + 1;
    m_count = 0;
    m_bin = new DictItem<T>*[ m_binCount ];
    assertDebug( m_bin );
    for (int i = 0; i < m_binCount; i++) 
        m_bin[ i ] = NULL;

    // clear arrays
    m_array.reset();
    m_keyArray.reset();
}


//-------------------------------------------
// STRING-KEY TEMPLATED DICTIONARY CLASS
//-------------------------------------------


/// A StringDictItem is a key-value pair in an string-keyed dictionary.
template<typename T> class StringDictItem {
public:
    T *item;
    int hash;
    String key;
    StringDictItem<T> *next;
    inline StringDictItem( const String &newKey, int newHash, T *newItem, StringDictItem<T> *newNext ) 
        { key = newKey; hash = newHash; item = newItem; next = newNext; }
};


/// The Dict class stores objects index by integer keys.
/// The Dict class uses a hash table for fast look-up.
/// (The class uses a power-of-2 number of bins, so certain key patterns could have many collisions, resulting in slow performance.)
template<typename T> class StringDict {
public:

    // basic constructor/destructor
    StringDict();
    ~StringDict();

    /// get item with given key; returns NULL if not found;
    /// does not remove or relinquish control of the object
    T *find( const String &key ); 
    const T *find( const String &key ) const; 

    /// access as array
    inline const String &refArrayKey( int index ) const { return m_keyArray[ index ]; }
    inline const T &refArray( int index ) const { return m_array[ index ]; }
    inline T &refArray( int index ) { return m_array[ index ]; }
    inline int count() const { return m_array.count(); }

    /// add an item to the dictionary; takes ownership of pointer unless dealloc disabled;
    /// replaces item with same key if already exists
    void add( const String &key, T *obj );

    /// disables object deallocation; the dictionary no longer takes ownership of object pointers
    void disableObjectDealloc() { m_deallocObjs = false; }

    /// remove all items (as if call destructor then constructor)
    void reset();

private:

    // if true, dictionary owns the objects and should delete them on destruction
    bool m_deallocObjs;

    // the number of items in the dictionary
    int m_count;

    // each bin contains a list of items
    StringDictItem<T> **m_bin;

    // the number of bins
    int m_binCount;

    // a bit mask used to obtain a bin index from an item ID
    int m_binMask;

    // an array of items in the dictionary (for iterating over the members of a dictionary)
    PtrArray<T> m_array; // use non-owning array because dict will own the items
    Array<String> m_keyArray;

    // disable copy constructor and assignment operator
    StringDict( const StringDict &x );
    StringDict &operator=( const StringDict &x );
};


// basic constructor
template <typename T> StringDict<T>::StringDict() {
    m_deallocObjs = true;
    m_binMask = 0xffff;
    m_binCount = m_binMask + 1;
    m_count = 0;
    m_bin = new StringDictItem<T>*[ m_binCount ];
    assertDebug( m_bin );
    for (int i = 0; i < m_binCount; i++) 
        m_bin[ i ] = NULL;
}


// basic destructor
template <typename T> StringDict<T>::~StringDict() {
    for (int i = 0; i < m_binCount; i++) {
        StringDictItem<T> *bin = m_bin[ i ];
        while (bin) {
            StringDictItem<T> *delItem = bin;
            bin = bin->next;
            if (m_deallocObjs)
                delete delItem->item;
            delete delItem;
        }
    }
    delete [] m_bin;
}


/// get item with given key; returns NULL if not found;
/// does not remove or relinquish control of the object
template <typename T> T *StringDict<T>::find( const String &key ) {
    int hash = strHash( key );
    StringDictItem<T> *bin = m_bin[ hash & m_binMask ];
    while (bin) {
        if (bin->hash == hash && bin->key == key) // note: using string equality operator
            return bin->item;
        bin = bin->next;
    }
    return NULL;
}


/// get item with given key; returns NULL if not found;
/// does not remove or relinquish control of the object
template <typename T> const T *StringDict<T>::find( const String &key ) const {
    int hash = strHash( key );
    StringDictItem<T> *bin = m_bin[ hash & m_binMask ];
    while (bin) {
        if (bin->hash == hash && bin->key == key) // note: using string equality operator
            return bin->item;
        bin = bin->next;
    }
    return NULL;
}


/// add an item to the dictionary; takes ownership of pointer unless dealloc disabled;
/// replaces item with same key if already exists
template <typename T> void StringDict<T>::add( const String &key, T *obj ) {
    assertDebug( obj );
    int hash = strHash( key );

    // look for existing item
    int binIndex = hash & m_binMask;
    StringDictItem<T> *bin = m_bin[ binIndex ];
    while (bin) {
        if (bin->hash == hash && bin->key == key) { // note: using string equality operator
            assertAlways( true ); // fix(later): does not update array
            if (m_deallocObjs)
                delete bin->item;
            bin->item = obj;
            return;
        }
        bin = bin->next;
    }
    
    // create new item
    m_bin[ binIndex ] = new StringDictItem<T>( key, hash, obj, m_bin[ binIndex ] );
    assertDebug( m_bin[ binIndex ] );
    m_count++;

    // add to array
    m_array.append( obj );
    m_keyArray.appendCopy( key );
}


/// remove all items (as if call destructor then constructor)
template <typename T> void StringDict<T>::reset() {

    // delete everything
    for (int i = 0; i < m_binCount; i++) {
        StringDictItem<T> *bin = m_bin[ i ];
        while (bin) {
            StringDictItem<T> *delItem = bin;
            bin = bin->next;
            if (m_deallocObjs)
                delete delItem->item;
            delete delItem;
        }
    }
    delete [] m_bin;

    // re-init
    m_deallocObjs = true;
    m_binMask = 0xffff;
    m_binCount = m_binMask + 1;
    m_count = 0;
    m_bin = new StringDictItem<T>*[ m_binCount ];
    assertDebug( m_bin );
    for (int i = 0; i < m_binCount; i++) 
        m_bin[ i ] = NULL;

    // clear arrays
    m_array.reset();
    m_keyArray.reset();
}


} // end namespace sbl
#endif

