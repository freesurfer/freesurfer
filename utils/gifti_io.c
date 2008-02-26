#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "gifti_io.h"

/*! global history and version strings, for printing */
static char * gifti_history[] =
{
  "----------------------------------------------------------------------\n"
  "history (of gifti library changes):\n"
  "\n",
  "0.0  18 July, 2007\n"
  "     (Rick Reynolds of the National Institutes of Health, SSCC/DIRP/NIMH)\n"
  "     - initial version\n"
  "0.1  31 July, 2007\n",
  "     - changed dim0..dim5 to dims[]\n"
  "     - changed nvals to size_t\n"
  "     - added gifti_init_darray_from_attrs and some validation functions\n"
  "0.2  29 October, 2007\n",
  "     - renamed gifti.[ch] to gifti_io.[ch]\n"
  "     - main data structures all start with gii (or gifti_)\n"
  "     - added user indenting\n"
  "0.3  21 November, 2007\n",
  "     - added base64 encoding/decoding, via b64_en/decode_table\n"
  "     - added gifti_list_index2string, gifti_disp_hex_data, \n"
  "             gifti_check_swap, gifti_swap_Nbytes, etc.\n"
  "     - pop_darray: check for b64 errors and byte swapping\n"
  "     - dind is size_t\n",
  "0.4  29 November, 2007\n"
  "     - added more checks and fixed nvpair value allocation\n"
  "0.5  03 December, 2007: applied changes for GIFTI Format 1.0 (11/21)\n",
  "     - replaced Category with Intent\n"
  "     - replaced Location attribute with ExternalFileName/Offset\n"
  "     - added NumberOfDataArrays attribute to GIFTI element\n"
  "     - applied new index_order strings\n"
  "0.6  10 December, 2007:\n",
  "     - can read/write Base64Binary datasets (can set compress level)\n"
  "     - removed datatype lists (have gifti_type_list)\n"
  "     - added gifti_read_da_list(), with only partial ability\n"
  "     - added GIFTI numDA attribute\n"
  "     - change size_t to long long\n"
  "0.7  11 December, 2007:\n",
  "     - added GIFTI_B64_CHECK defines\n"
  "     - set b64_check default to SKIPNCOUNT\n"
  "     - added disp_gxml_data\n"
  "0.8  12 December, 2007:\n",
  "     - added sub-surface selection, via dalist in gifti_read_da_list()\n"
  "     - added gifti_copy_DataArray, and other structures\n"
  "0.9  28 December, 2007:\n",
  "     - made zlib optional, via -DHAVE_ZLIB in compile\n"
  "       (without zlib, the user will get warnings)\n"
  "     - now users only #include gifti_io.h, not gifti_xml, expat or zlib\n"
  "     - added more comments and made tables more readable\n"
  "     - added all user-variable access functions and reset_user_vars()\n",
  "     - added gifti_free_image_contents(), gifti_disp_raw_data(),\n"
  "             gifti_clear_float_zeros() and gifti_set_DA_atrs()\n"
  "     - changed gifti_gim_DA_size to long long\n"
  "     - added GIFTI_B64_CHECK_UNDEF as 0\n"
  "     - fixed 0-width indenting and accumulating base64 errors\n"
  "0.10 03 January, 2008:\n",
  "     - added top-level gifti_create_image() interface\n"
  "     - must now link libniftiio\n"
  "     - gifti_add_empty_darray() now takes num_to_add\n"
  "     - if data was expected but not read, free it\n"
  "         (can add via gifti_alloc_DA_data())\n"
  "     - many minor changes\n"
  "0.11 11 January, 2008:\n",
  "     - attribute/data setting functions are more flexible\n"
  "     - added gifti_disp_dtd_url, gifti_set_DA_meta, gifti_valid_int_list,\n"
  "       DA_data_exists, gifti_add_to_meta \n"
  "0.12 16 January, 2008:\n",
  "     - added gifti_copy_gifti_image() and gifti_convert_to_float()\n"
  "     - added gifti_valid_LabelTable(), gifticlib_version(),\n"
  "             gifti_copy_LabelTable(), gifti_updaet_nbyper() and\n"
  "             gifti_valid_gifti_image()\n"
  "     - added control over library updates to metadata\n"
  "     - expanded checks in gifti_valid_dims\n"
  "0.13 20 February, 2008:\n",
  "     - added gifti_get_meta_value() and gifti_image_has_data()\n"
  "0.14 25 February, 2008:\n",
  "     - consider data-less metadata as valid\n"
};

static char gifti_version[] = "gifti library version 0.14, 20 February, 2008";

/* ---------------------------------------------------------------------- */
/*! global lists of XML strings */

/*! this should match GIFTI_IND_ORD_* */
char * gifti_index_order_list[] = {"Undefined", "RowMajorOrder",
                                                "ColumnMajorOrder"};
/* {"Undefined", "HighestFirst", "LowestFirst"}; */

/* char * gifti_dataloc_list[] = {"Undefined", "Internal", "External"}; */

/*! gifti_type_list is an array of gifti_type_ele structs which list, for
    each type, the bytes per value, swapsize and corresponding name string
    (these type values are defined in nifti1.h) */
static gifti_type_ele gifti_type_list[] = {
    /* type                    nbyper  swapsize   name  */
    { DT_UNKNOWN,                 0,      0,      "Undefined"             },
    { NIFTI_TYPE_UINT8,           1,      0,      "NIFTI_TYPE_UINT8"      },
    { NIFTI_TYPE_INT16,           2,      2,      "NIFTI_TYPE_INT16"      },
    { NIFTI_TYPE_INT32,           4,      4,      "NIFTI_TYPE_INT32"      },
    { NIFTI_TYPE_FLOAT32,         4,      4,      "NIFTI_TYPE_FLOAT32"    },
    { NIFTI_TYPE_COMPLEX64,       8,      4,      "NIFTI_TYPE_COMPLEX64"  },
    { NIFTI_TYPE_FLOAT64,         8,      8,      "NIFTI_TYPE_FLOAT64"    },
    { NIFTI_TYPE_RGB24,           3,      0,      "NIFTI_TYPE_RGB24"      },
    { NIFTI_TYPE_INT8,            1,      0,      "NIFTI_TYPE_INT8"       },
    { NIFTI_TYPE_UINT16,          2,      2,      "NIFTI_TYPE_UINT16"     },
    { NIFTI_TYPE_UINT32,          4,      4,      "NIFTI_TYPE_UINT32"     },
    { NIFTI_TYPE_INT64,           8,      8,      "NIFTI_TYPE_INT64"      },
    { NIFTI_TYPE_UINT64,          8,      8,      "NIFTI_TYPE_UINT64"     },
    { NIFTI_TYPE_FLOAT128,        6,     16,      "NIFTI_TYPE_FLOAT128"   },
    { NIFTI_TYPE_COMPLEX128,     16,      8,      "NIFTI_TYPE_COMPLEX128" },
    { NIFTI_TYPE_COMPLEX256,     32,     16,      "NIFTI_TYPE_COMPLEX256" }
};

/*! this list provides a link between intent codes and their name strings */
typedef struct { int code; char * name; } gifti_intent_ele;
static gifti_intent_ele gifti_intent_list[] = {
    { NIFTI_INTENT_NONE,             "NIFTI_INTENT_NONE"        },
    { NIFTI_INTENT_CORREL,           "NIFTI_INTENT_CORREL"      },
    { NIFTI_INTENT_TTEST,            "NIFTI_INTENT_TTEST"       },
    { NIFTI_INTENT_FTEST,            "NIFTI_INTENT_FTEST"       },
    { NIFTI_INTENT_ZSCORE,           "NIFTI_INTENT_ZSCORE"      },
    { NIFTI_INTENT_CHISQ,            "NIFTI_INTENT_CHISQ"       },
    { NIFTI_INTENT_BETA,             "NIFTI_INTENT_BETA"        },
    { NIFTI_INTENT_BINOM,            "NIFTI_INTENT_BINOM"       },
    { NIFTI_INTENT_GAMMA,            "NIFTI_INTENT_GAMMA"       },
    { NIFTI_INTENT_POISSON,          "NIFTI_INTENT_POISSON"     },
    { NIFTI_INTENT_NORMAL,           "NIFTI_INTENT_NORMAL"      },
    { NIFTI_INTENT_FTEST_NONC,       "NIFTI_INTENT_FTEST_NONC"  },
    { NIFTI_INTENT_CHISQ_NONC,       "NIFTI_INTENT_CHISQ_NONC"  },
    { NIFTI_INTENT_LOGISTIC,         "NIFTI_INTENT_LOGISTIC"    },
    { NIFTI_INTENT_LAPLACE,          "NIFTI_INTENT_LAPLACE"     },
    { NIFTI_INTENT_UNIFORM,          "NIFTI_INTENT_UNIFORM"     },
    { NIFTI_INTENT_TTEST_NONC,       "NIFTI_INTENT_TTEST_NONC"  },
    { NIFTI_INTENT_WEIBULL,          "NIFTI_INTENT_WEIBULL"     },
    { NIFTI_INTENT_CHI,              "NIFTI_INTENT_CHI"         },
    { NIFTI_INTENT_INVGAUSS,         "NIFTI_INTENT_INVGAUSS"    },
    { NIFTI_INTENT_EXTVAL,           "NIFTI_INTENT_EXTVAL"      },
    { NIFTI_INTENT_PVAL,             "NIFTI_INTENT_PVAL"        },
    { NIFTI_INTENT_LOGPVAL,          "NIFTI_INTENT_LOGPVAL"     },
    { NIFTI_INTENT_LOG10PVAL,        "NIFTI_INTENT_LOG10PVAL"   },
    { NIFTI_INTENT_ESTIMATE,         "NIFTI_INTENT_ESTIMATE"    },
    { NIFTI_INTENT_LABEL,            "NIFTI_INTENT_LABEL"       },
    { NIFTI_INTENT_NEURONAME,        "NIFTI_INTENT_NEURONAME"   },
    { NIFTI_INTENT_GENMATRIX,        "NIFTI_INTENT_GENMATRIX"   },
    { NIFTI_INTENT_SYMMATRIX,        "NIFTI_INTENT_SYMMATRIX"   },
    { NIFTI_INTENT_DISPVECT,         "NIFTI_INTENT_DISPVECT"    },
    { NIFTI_INTENT_VECTOR,           "NIFTI_INTENT_VECTOR"      },
    { NIFTI_INTENT_POINTSET,         "NIFTI_INTENT_POINTSET"    },
    { NIFTI_INTENT_TRIANGLE,         "NIFTI_INTENT_TRIANGLE"    },
    { NIFTI_INTENT_QUATERNION,       "NIFTI_INTENT_QUATERNION"  },
    { NIFTI_INTENT_DIMLESS,          "NIFTI_INTENT_DIMLESS"     },
    { NIFTI_INTENT_TIME_SERIES,      "NIFTI_INTENT_TIME_SERIES" },
    { NIFTI_INTENT_NODE_INDEX,       "NIFTI_INTENT_NODE_INDEX"  },
    { NIFTI_INTENT_RGB_VECTOR,       "NIFTI_INTENT_RGB_VECTOR"  },
    { NIFTI_INTENT_RGBA_VECTOR,      "NIFTI_INTENT_RGBA_VECTOR" },
    { NIFTI_INTENT_SHAPE,            "NIFTI_INTENT_SHAPE"       }
};

/*! this should match GIFTI_ENCODING_* */
char * gifti_encoding_list[] = {
    "Undefined", "ASCII", "Base64Binary", "GZipBase64Binary",
    "ExternalFileBinary"
};

/*! this should match GIFTI_ENDIAN_* */
char * gifti_endian_list[] = {"Undefined", "BigEndian", "LittleEndian"};

/* ---------------------------------------------------------------------- */
/* local prototypes */
static int str2list_index(char *list[], int max, const char *str);
static int DA_data_exists(gifti_image * gim, const int * dalist, int len);
static int copy_data_as_float(void * dest, int dtype, void * src, int stype,
                              long long nvals);

/* ---------------------------------------------------------------------- */
/*! giftilib globals */
static gifti_globals G = { 1 };

/* ====================================================================== */

/* ---------------------------------------------------------------------- */
/*! user variable accessor functions - basically use gxml interface       */

int gifti_get_verb( void )          { return G.verb; }
int gifti_set_verb( int level )     { G.verb = level;  return 1; }
int gifti_get_indent( void )        { return gxml_get_indent(); }
int gifti_set_indent( int level )   { return gxml_set_indent(level); }
int gifti_get_b64_check( void )     { return gxml_get_b64_check(); }
int gifti_set_b64_check( int level ){ return gxml_set_b64_check(level); }
int gifti_get_update_ok( void )     { return gxml_get_update_ok(); }
int gifti_set_update_ok( int level ){ return gxml_set_update_ok(level); }
int gifti_get_zlevel( void )        { return gxml_get_zlevel(); }
int gifti_set_zlevel( int level )
{
    /* note that the default currently results in 6 */
    if( level != GZ_DEFAULT_COMPRESSION && (level < 0 || level > 9 ) ) {
        fprintf(stderr,"** invalid zlevel, must be %d (default) or {0..9}\n",
                GZ_DEFAULT_COMPRESSION);
        return 1;
    }
    return gxml_set_zlevel(level);
}

int gifti_get_xml_buf_size(void)        { return gxml_get_buf_size(); }
int gifti_set_xml_buf_size(int buf_size){ return gxml_set_buf_size(buf_size); }

/*! reset user variables to their defaults(via set to -1) */
int gifti_reset_user_vars(void)
{
    gxml_set_verb(-1);
    gxml_set_dstore(-1);
    gxml_set_indent(-1);
    gxml_set_buf_size(-1);
    gxml_set_b64_check(-1);
    gxml_set_update_ok(-1);
    gxml_set_zlevel(-1);

    return 0;
}
/* end user variable accessor functions                                   */
/* ---------------------------------------------------------------------- */

/* ====================================================================== */

/* ---------------------------------------------------------------------- */
/* begin general library functions                                        */

/*----------------------------------------------------------------------
 *! apply the attr=value GIFTI attribute to the gifti_image
 *
 *  return 0 on success
*//*-------------------------------------------------------------------*/
int gifti_str2attr_gifti(gifti_image * gim, const char *attr, const char *val)
{
    if( !gim || !attr || !val ) {
        fprintf(stderr,"** GS2AG: bad params (%p,%p,%p)\n",
                (void *)gim, (void *)attr, (void *)val);
        return 1;
    }

    if( G.verb > 2 )
        fprintf(stderr,"++ setting GIFTI attr '%s' from '%s'\n", attr, val);

    if( !strcmp(attr, "Version") ) {
        if( gim->version ) free( gim->version );  /* lose any old copy */
        gim->version = gifti_strdup(val);
    } else if( !strcmp(attr, "NumberOfDataArrays") ) {
        gim->numDA = atol(val);
        if( gim->numDA < 0 ) {
            fprintf(stderr,"** invalid NumberOfDataArrays attribute: %s\n",val);
            gim->numDA = 0;
            return 1;
        }
    } else {
        if( G.verb > 0 )
            fprintf(stderr,"** unknown GIFTI attrib, '%s'='%s'\n",attr,val);
        return 1;
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! This is the main dataset reading routine.  Read a GIFTI dataset
 *  and return the corresponding gifti_image structure.
 *
 *  Reading data is optional, via the read_data flag.
 *  User variables should already be set (via accessor functions).
 *
 *  return an allocated gifti_image struct on success,
 *         NULL on error
*//*-------------------------------------------------------------------*/
gifti_image * gifti_read_image( const char * fname, int read_data )
{
    if( !fname ) {
        fprintf(stderr,"** gifti_read_image: missing filename\n");
        return NULL;
    }

    gxml_set_verb(G.verb);

    return gxml_read_image(fname, read_data, NULL, 0);
}

/*----------------------------------------------------------------------
 *! Similar to gifti_read_data, this function also takes an integer list of
 *  DataArray indices to populate the gifti_image structure with.
 *
 *  The indices are be zero-based, can have repeats and can be in any order.
 *  A simple example to read 3 DA elements (with 2 repeats) might be:
 *
 *    gifti_image * gim;
 *    int           ilist[5] = { 3, 0, 7, 7, 3 };
 *    gim = gifti_read_da_list("my_data.gii", 1, ilist, 5);
 *
 *  return an allocated gifti_image struct on success,
 *         NULL on error
*//*-------------------------------------------------------------------*/
gifti_image * gifti_read_da_list( const char * fname, int read_data,
                                  const int * dalist, int len )
{
    if( !fname ) {
        fprintf(stderr,"** gifti_read_da_list: missing filename\n");
        return NULL;
    }

    gxml_set_verb(G.verb);

    return gxml_read_image(fname, read_data, dalist, len);
}

/*----------------------------------------------------------------------
 *! This is the main dataset writing routine.
 *
 *  User variables should be set before this point.
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_write_image(gifti_image *gim, const char *fname, int write_data)
{
    int errs = 0;

    if( !gim ) {
        fprintf(stderr,"** gifti_write_image, missing gifti_image\n");
        errs++;
    } else if( !fname ) {
        fprintf(stderr,"** gifti_write_image: missing filename\n");
        errs++;
    }

    if( errs ) return 1;

    gxml_set_verb(G.verb);

    return gxml_write_image(gim, fname, write_data);
}


/*----------------------------------------------------------------------
 *! free the gifti_image struct and all its contents
 *
 *  passing NULL (to this and any child function) should be okay
 *
 *  the pointer is garbage after this call
*//*-------------------------------------------------------------------*/
int gifti_free_image( gifti_image * gim )
{
    if( !gim ) {
        if(G.verb > 2) fprintf(stderr,"** free gifti_image w/NULL pointer\n");
        return 1;
    }
    
    if( G.verb > 2 ) fprintf(stderr,"-- freeing gifti_image\n");

    if( gim->version ) { free(gim->version);  gim->version = NULL; }

    (void)gifti_free_nvpairs(&gim->meta);
    (void)gifti_free_LabelTable(&gim->labeltable);
    (void)gifti_free_DataArray_list(gim->darray, gim->numDA);
    (void)gifti_free_nvpairs(&gim->ex_atrs);
    free(gim);

    return 0;
}

/*----------------------------------------------------------------------
 *! free the contents of the gifti_image struct (but not the pointer)
 *
 *  the pointer is garbage after this call
*//*-------------------------------------------------------------------*/
int gifti_free_image_contents( gifti_image * gim )
{
    if( !gim ) {
        if(G.verb > 2) fprintf(stderr,"** GFIC: free w/NULL gifti_image ptr\n");
        return 1;
    }
    
    if( G.verb > 2 ) fprintf(stderr,"-- freeing gifti_image contents\n");

    if( gim->version ) { free(gim->version);  gim->version = NULL; }

    (void)gifti_free_nvpairs(&gim->meta);
    (void)gifti_free_LabelTable(&gim->labeltable);
    (void)gifti_free_DataArray_list(gim->darray, gim->numDA);
    (void)gifti_free_nvpairs(&gim->ex_atrs);

    return 0;
}

/*----------------------------------------------------------------------
 *! free the contents of the nvpairs struct (but not the pointer)
 *
 *  passing NULL is okay
*//*-------------------------------------------------------------------*/
int gifti_free_nvpairs( nvpairs * p )
{
    int c;

    if( !p ) {
        if( G.verb > 3 ) fprintf(stderr,"** free w/NULL nvpairs ptr\n");
        return 1;
    }

    if( G.verb > 3 ) fprintf(stderr,"-- freeing %d nvpairs\n", p->length);

    if( p->name && p->value ) {
        for( c = 0; c < p->length; c++ ) {
            if( p->name[c] ) free(p->name[c]);
            if( p->value[c] ) free(p->value[c]);
        }
        free(p->name);
        free(p->value);
        p->name = NULL;
        p->value = NULL;
    }
    p->length = 0;

    return 0;
}

/*----------------------------------------------------------------------
 *! free the contents of the LabelTable struct (but not the pointer)
 *
 *  passing NULL is okay
*//*-------------------------------------------------------------------*/
int gifti_free_LabelTable( giiLabelTable * T )
{
    int c;

    if( !T ) {
        if(G.verb > 3) fprintf(stderr,"** free w/NULL giiLabelTable ptr\n");
        return 1;
    }

    if(G.verb > 3)
        fprintf(stderr,"-- freeing %d giiLabelTables\n", T->length);

    if( T->index && T->label ) {
        for( c = 0; c < T->length; c++ )
            if( T->label[c] ) free(T->label[c]);
        free(T->index);
        free(T->label);
        T->index = NULL;
        T->label = NULL;
    }
    T->length = 0;

    return 0;
}

/*----------------------------------------------------------------------
 *! free the DataArray list (the array and all its contents)
 *
 *  the darray list pointer is garbage after this call
 *
 *  passing NULL is okay
*//*-------------------------------------------------------------------*/
int gifti_free_DataArray_list(giiDataArray ** darray, int numDA)
{
    int c;

    if( !darray ) {
        if( G.verb > 3 ) fprintf(stderr,"** GFDA: free NULL darray list\n");
        return 1;
    }

    if( G.verb > 3 ) fprintf(stderr,"-- freeing %d giiDataArrays\n", numDA);

    if( numDA < 0 ) return 1;

    for( c = 0; c < numDA; c++ )
        if( gifti_free_DataArray(darray[c]) ) return 1;

    free(darray);

    return 0;
}

/*----------------------------------------------------------------------
 *! free the DataArray struct and all its contents
 *
 *  the DataArray pointer is garbage after this call
 *
 *  passing NULL is okay
*//*-------------------------------------------------------------------*/
int gifti_free_DataArray( giiDataArray * darray )
{
    if( !darray ) {
        if( G.verb > 3 ) fprintf(stderr,"** tried to free NULL darray ptr\n");
        return 1;
    }

    if( G.verb > 3 ) fprintf(stderr,"-- freeing giiDataArray\n");

    if(darray->ext_fname) { free(darray->ext_fname); darray->ext_fname = NULL; }

    (void)gifti_free_nvpairs(&darray->meta);
    if( darray->coordsys) {
        gifti_free_CoordSystem(darray->coordsys);
        darray->coordsys = NULL;
    }
    if( darray->data ) { free(darray->data); darray->data = NULL; }
    (void)gifti_free_nvpairs(&darray->ex_atrs);
    free(darray);

    return 0;
}

/*----------------------------------------------------------------------
 *! free the CoordSystem struct and all its contents
 *
 *  the CoordSystem pointer is garbage after this call
 *
 *  passing NULL is okay
*//*-------------------------------------------------------------------*/
int gifti_free_CoordSystem( giiCoordSystem * cs )
{
    if( !cs ) return 0;

    if( G.verb > 3 ) fprintf(stderr,"-- freeing giiCoordSystem\n");

    if( cs->dataspace ) { free(cs->dataspace); cs->dataspace = NULL; }
    if( cs->xformspace ) { free(cs->xformspace); cs->xformspace = NULL; }

    free(cs);

    return 0;
}

/*----------------------------------------------------------------------
 *! initialize an existing DataArray structure given a list of name=value
 *  attribute pairs
 *
 *  if alen > 0, consider it the length of the attr list
 *  else         process until attr[i] == NULL
 *
 *  if add_to_extras, add any bad attribute pairs to ex_atrs
 *  else              whine about any bad ones and return
*//*-------------------------------------------------------------------*/
int gifti_set_DA_atrs(giiDataArray * da, const char ** attr, int alen,
                      int add_to_extras )
{
    int c, length = alen;

    if( !da || !attr ) {
        if(G.verb>1) fprintf(stderr,"** G_IDFA: bad params (%p,%p)\n",
                             (void *)da,(void *)attr);
        return 1;
    }

    /* if length was not passed, compute it */
    if( length <= 0 ) for( length = 0; attr[length]; length++ ) /* nada */ ;

    if( G.verb > 5 )
        fprintf(stderr,"++ init darray attrs, len %d, ex_atrs = %d\n",
                length, add_to_extras);

    /* insert attributes - if unknown, store with extras */
    for(c = 0; c < length; c += 2 )
        if( gifti_str2attr_darray(da, attr[c],attr[c+1]) ) {
            /* a bad name=value pair, maybe add to ex_atrs */
            if( add_to_extras ) {
                if( gifti_add_to_nvpairs(&da->ex_atrs,attr[c],attr[c+1]) )
                    return 1;
            }
            else {
                if( G.verb > 0 )
                    fprintf(stderr,"** set_darray_atrs, bad pair '%s'='%s'\n",
                            attr[c],attr[c+1]);
                return 1;
            }
        }

    /* and set computed values */

    da->nvals = gifti_darray_nvals(da);
    gifti_datatype_sizes(da->datatype, &da->nbyper, NULL); /* set nbyper */

    return 0;
}

/*----------------------------------------------------------------------
 *! determine whether the given DataArray struct seems valid
 *
 *  if whine is set, print error messages for any failures
 *
 *  return 1, if valid
 *         0, if not
*//*-------------------------------------------------------------------*/
int gifti_valid_DataArray(giiDataArray * da, int whine)
{
    int errs = 0, nbyper;

    if( !da ) {
        if( whine || G.verb > 1 ) fprintf(stderr,"** invalid darray pointer\n");
        return 0;
    }

    if( ! gifti_intent_is_valid(da->intent) ) {
        if( whine || G.verb > 1 )
            fprintf(stderr,"** invalid darray intent code = %d\n", da->intent);
        errs++;
    }

    if( ! gifti_valid_datatype(da->datatype, whine) ) /* message printed */
        errs++;

    /* no checks for ext_fname and ext_offset (until reading) */

    if( da->ind_ord<=GIFTI_IND_ORD_UNDEF || da->ind_ord>GIFTI_IND_ORD_MAX ) {
        if( whine || G.verb > 1 )
            fprintf(stderr,"** invalid darray ind_ord = %d\n", da->ind_ord);
        errs++;
    }

    if( ! gifti_valid_num_dim(da->num_dim, whine) ) /* message printed */
        errs++;

    if( ! gifti_valid_dims(da, whine) ) /* message printed */
        errs++;

    if( da->encoding<=GIFTI_ENCODING_UNDEF || da->encoding>GIFTI_ENCODING_MAX ){
        if( whine || G.verb > 1 )
            fprintf(stderr,"** invalid darray encoding = %d\n", da->encoding);
        errs++;
    }

    if( da->endian<=GIFTI_ENDIAN_UNDEF || da->endian>GIFTI_ENDIAN_MAX ) {
        if( whine || G.verb > 1 )
            fprintf(stderr,"** invalid darray endian = %d\n", da->endian);
        errs++;
    }

    /* of sub-element, only verify giiMetaData */
    if( ! gifti_valid_nvpairs(&da->meta, whine) ) /* message printed */
        errs++;

    if( da->nvals <= 0 ) {
        if( whine || G.verb > 1 )
            fprintf(stderr,"** invalid darray nvals = %u\n",
                    (unsigned)da->nvals );
        errs++;
    }

    if( ! gifti_valid_nbyper(da->nbyper, whine) ) errs++;
    if( ! gifti_valid_nvpairs(&da->ex_atrs, whine) ) errs++;

    /* compare nbyper to what is expected for type */
    errs += gifti_datatype_sizes(da->datatype, &nbyper, NULL);
    if( gifti_valid_nbyper(nbyper, 0) && nbyper != da->nbyper ) {
        if( whine || G.verb > 1 )
            fprintf(stderr,"** nbyper %d, does not match type %s\n",
                    nbyper, gifti_datatype2str(da->datatype));
        errs++;
    }

    if( errs ) return 0;

    return 1;
}

/*----------------------------------------------------------------------
 *! check whether pointers are valid and consistent with length
*//*-------------------------------------------------------------------*/
int gifti_valid_nvpairs(nvpairs * nvp, int whine)
{
    int c;

    if( !nvp ) {
        if( G.verb>0 || whine ) fprintf(stderr,"** invalid nvpairs pointer\n");
        return 0;
    }

    if( nvp->length < 0 ) {
        if( G.verb > 1 || whine )
            fprintf(stderr,"** invalid nvpair length = %d\n", nvp->length);
        return 0;
    }

    if( nvp->length == 0 ) return 1;    /* quick case: valid */

    if( !nvp->name || !nvp->value ){
        if( G.verb > 1 || whine )
            fprintf(stderr,"** invalid nvpair name, value lists = %p, %p\n",
                    (void *)nvp->name, (void *)nvp->value);
        return 0;
    }

    /* quit on first error */
    for( c = 0; c < nvp->length; c++ ) {
        if( ! nvp->name[c] ) {
            if( G.verb > 2 || whine )
                fprintf(stderr,"** invalid nvpair, missing name @ %d\n", c);
            return 0;
        }

        /* value string is not required   25 Feb 2008 */
        if( ! nvp->value[c] && G.verb > 2 )
            fprintf(stderr,"-- missing nvpair value[%d], name %s (is OK)\n",
                    c, nvp->name[c]);
    }

    return 1;
}

/*----------------------------------------------------------------------
 *! check whether pointers are valid and consistent with length
 *
 *  no check is done on the actual indices or labels
*//*-------------------------------------------------------------------*/
int gifti_valid_LabelTable(giiLabelTable * T, int whine)
{
    int c;

    if( !T ) {
        if(G.verb>0||whine) fprintf(stderr,"** invalid LabelTable pointer\n");
        return 0;
    }

    if( T->length < 0 ) {
        if( G.verb > 1 || whine )
            fprintf(stderr,"** invalid LabelTable length = %d\n", T->length);
        return 0;
    }

    if( T->length == 0 ) return 1;    /* quick case: valid */

    if( !T->index || !T->label ){
        if( G.verb > 1 || whine )
            fprintf(stderr,"** invalid nvpair index, label = %p, %p\n",
                    (void *)T->index, (void *)T->label);
        return 0;
    }

    /* quit on first error */
    for( c = 0; c < T->length; c++ ) {
        if( ! T->label[c] ) {
            if( G.verb > 1 || whine )
                fprintf(stderr,"** invalid nvpair label[%d]\n", c);
            return 0;
        }
    }

    return 1;
}

/*----------------------------------------------------------------------
 *! check the bounds on num_dim
*//*-------------------------------------------------------------------*/
int gifti_valid_num_dim(int num_dim, int whine)
{
    if( num_dim <= 0 || num_dim > GIFTI_DARRAY_DIM_LEN ) {
        if( G.verb > 1 || whine )
            fprintf(stderr,"** invalid num_dim = %d\n", num_dim);
        return 0;
    }
    return 1;
}

/*----------------------------------------------------------------------
 *! check that the datatype is in the list
*//*-------------------------------------------------------------------*/
int gifti_valid_datatype(int dtype, int whine)
{
    int c;

    /* check for valid */
    for( c = sizeof(gifti_type_list) / sizeof(gifti_type_ele) - 1; c > 0; c-- )
        if( dtype == gifti_type_list[c].type ) return 1;

    if( whine || G.verb > 1 )
        fprintf(stderr,"** invalid datatype value %d\n", dtype);

    return 0;
}

/*----------------------------------------------------------------------
 *! check that nbyper is one of the values in gifti_type_list
*//*-------------------------------------------------------------------*/
int gifti_valid_nbyper(int nbyper, int whine)
{
    int c;

    /* check for valid */
    for( c = sizeof(gifti_type_list) / sizeof(gifti_type_ele) - 1; c > 0; c-- )
        if( nbyper == gifti_type_list[c].nbyper ) return 1;

    if( whine || G.verb > 1 )
        fprintf(stderr,"** invalid nbyper value %d\n", nbyper);

    return 0;
}

/*----------------------------------------------------------------------
 *! check that dimension values are consistent (and with datatype)
 *  
 *      - num_dim is in range
 *      - each dims[c] is postive (c < num_dim)
 *      - nvals is product of dims
 *      - datatype is valie (required to check nbyper)
 *      - nbyper is correct
*//*-------------------------------------------------------------------*/
int gifti_valid_dims(giiDataArray * da, int whine)
{
    long long vals = 1;
    int       c, nbyper;

    if( !da ) {
        if( G.verb > 0 ) fprintf(stderr,"** GVD: no giiDataArray\n");
        return 0;
    }

    if( ! gifti_valid_num_dim( da->num_dim, whine ) )
        return 0;

    for( c = 0; c < da->num_dim; c++ ) {
        if( da->dims[c] <= 0 ) {
            if( G.verb > 1 || whine )
                fprintf(stderr,"** invalid dims[%d] = %d\n", c, da->dims[c]);
            return 0;
        }

        vals *= da->dims[c];
    }

    /* each DA must have the same length (hmmm, we should check all dims...) */
    if( vals != da->nvals ) {
        if( G.verb > 0 ) {
            fprintf(stderr,"** nvals = %lld does not match %lld for dims[%d]: ",
                    da->nvals, vals, da->num_dim);
            gifti_disp_raw_data(da->dims, DT_INT32, da->num_dim, 1, stderr);
        }
        return 0;
    }

    gifti_datatype_sizes(da->datatype, &nbyper, NULL);
    if( nbyper != da->nbyper ) {
        fprintf(stderr,"** nbyper %d not correct for type %s\n",
                da->nbyper, gifti_datatype2str(da->datatype));
    }

    return 1;
}

/*----------------------------------------------------------------------
 *! set the DataArray attribute, based on the name=value string pair
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_str2attr_darray(giiDataArray * DA, const char *attr,
                                             const char *value)
{
    if( !DA || !attr || !value ) {
        if( G.verb > 0 )
            fprintf(stderr,"** G_S2A_D: bad params (%p,%p,%p)\n",
                    (void *)DA, (void *)attr, (void *)value);
        return 1;
    }

    if(G.verb > 3) fprintf(stderr,"++ setting DA attr '%s'='%s'\n",attr,value);

    if( !strcmp(attr, "Intent") )
        DA->intent = gifti_intent_from_string(value);
    else if( !strcmp(attr, "DataType") )
        DA->datatype = gifti_str2datatype(value);
    else if( !strcmp(attr, "ArrayIndexingOrder") )
        DA->ind_ord = gifti_str2ind_ord(value);
    else if( !strcmp(attr, "Dimensionality") ) DA->num_dim = atoi(value);
    else if( !strcmp(attr, "Dim0") )           DA->dims[0] = atoi(value);
    else if( !strcmp(attr, "Dim1") )           DA->dims[1] = atoi(value);
    else if( !strcmp(attr, "Dim2") )           DA->dims[2] = atoi(value);
    else if( !strcmp(attr, "Dim3") )           DA->dims[3] = atoi(value);
    else if( !strcmp(attr, "Dim4") )           DA->dims[4] = atoi(value);
    else if( !strcmp(attr, "Dim5") )           DA->dims[5] = atoi(value);
    else if( !strcmp(attr, "Encoding") ) 
        DA->encoding = gifti_str2encoding(value);
    else if( !strcmp(attr, "Endian") )
        DA->endian = gifti_str2endian(value);
    else if( !strcmp(attr, "ExternalFileName") )
        DA->ext_fname = gifti_strdup(value);
    else if( !strcmp(attr, "ExternalFileOffset") )
        DA->ext_offset = atoll(value);  /* assumes C99 */
    else {
        if( G.verb > 1 )        /* might go into ex_atrs */
            fprintf(stderr,"** unknown giiDataArray attr, '%s'='%s'\n",
                    G_CHECK_NULL_STR(attr),G_CHECK_NULL_STR(value));
        return 1;
    }

    return 0;
}

/* return 0 (UNDEFINED) on failure */
static int str2list_index( char * list[], int max, const char *str )
{
    int index;
    if( !list || !str ) {
        if( G.verb > 0 ) fprintf(stderr,"** str2list: bad params (%p,%p)\n",
                                 (void *)list, (void *)str);
        return 0;  /* should be *_UNDEFINED */
    }

    for( index = max; index > 0; index-- )
        if( !strcmp(str, list[index]) ) return index;

    return 0;    /* failure */
}

/*----------------------------------------------------------------------
 *! return the index for a GIFTI_IND_ORD_* string
*//*-------------------------------------------------------------------*/
int gifti_str2ind_ord( const char * str )
{
    int rv = str2list_index(gifti_index_order_list,GIFTI_IND_ORD_MAX,str);
    if( rv <= GIFTI_IND_ORD_UNDEF && G.verb > 1 )
        fprintf(stderr,"** bad index order, '%s'\n", str);
    return rv;
}


/*----------------------------------------------------------------------
 *! return the index for a GIFTI_ENDODING_* string
*//*-------------------------------------------------------------------*/
int gifti_str2encoding( const char * str )
{
    int rv = str2list_index(gifti_encoding_list, GIFTI_ENCODING_MAX, str);
    if( rv <= GIFTI_ENCODING_UNDEF && G.verb > 1 )
        fprintf(stderr,"** bad data encoding, '%s'\n", str);
    return rv;
}

/*----------------------------------------------------------------------
 *! return the string at the given index of the given list
 *
 *  This function is meant to index into one of the gifti_*_list arrays,
 *  while being certain that the index is not out of range.
*//*-------------------------------------------------------------------*/
char * gifti_list_index2string(char * list[], int index)
{
    int lsize;    /* list size cannot be computed from the passed pointer */

    if     ( list == gifti_index_order_list )
        lsize = sizeof(gifti_index_order_list)/sizeof(char *);

    else if( list == gifti_encoding_list )
        lsize = sizeof(gifti_encoding_list)/sizeof(char *);

    else if( list == gifti_endian_list )
        lsize = sizeof(gifti_endian_list)/sizeof(char *);

    else {
        fprintf(stderr,"** GLI2S: invalid list\n");
        return "UNKNOWN LIST";
    }

    if( index < 0 || index >= lsize ) {
        if( G.verb > 0) 
            fprintf(stderr,"** GLI2S: index %d out of range {0..%d}\n",
                index,lsize-1);
        return "INDEX OUT OF RANGE";
    }

    /* all that work, just for the expected indexing...  */

    return list[index];
}

/*----------------------------------------------------------------------
 *! return the NIFTI_TYPE_ value corresponding to the given string
*//*-------------------------------------------------------------------*/
int gifti_str2datatype(const char * str)
{
    int len = sizeof(gifti_type_list)/sizeof(gifti_type_ele);
    int c;

    for( c = len - 1; c > 0; c-- )
        if( !strcmp(str, gifti_type_list[c].name) )
            break;

    return gifti_type_list[c].type;
}


/*----------------------------------------------------------------------
 *! return the GIFTI_ENDIAN_ value corresponding to the given string
*//*-------------------------------------------------------------------*/
int gifti_str2endian( const char * str )
{
    int rv = str2list_index(gifti_endian_list, GIFTI_ENDIAN_MAX, str);
    if( rv <= GIFTI_ENCODING_UNDEF && G.verb > 1 )
        fprintf(stderr,"** bad endian, '%s'\n", G_CHECK_NULL_STR(str));
    return rv;
}


/*----------------------------------------------------------------------
 *! return the NIFTI_TYPE_ value string corresponding to the given type
*//*-------------------------------------------------------------------*/
char * gifti_datatype2str(int type)
{
    int len = sizeof(gifti_type_list)/sizeof(gifti_type_ele);
    int c;

    for( c = len - 1; c > 0; c-- )
        if( type == gifti_type_list[c].type)
            break;

    return gifti_type_list[c].name;
}


/*----------------------------------------------------------------------
 *! simply set the struct contents to empty
*//*-------------------------------------------------------------------*/
int gifti_clear_nvpairs(nvpairs * p)
{
    if( !p ) return 1;

    p->length = 0;
    p->name = NULL;
    p->value = NULL;

    return 0;
}

/*----------------------------------------------------------------------
 *! simply set the struct contents to empty
*//*-------------------------------------------------------------------*/
int gifti_clear_LabelTable(giiLabelTable * p)
{
    if( !p ) return 1;

    p->length = 0;
    p->index = NULL;
    p->label = NULL;

    return 0;
}

/*----------------------------------------------------------------------
 *! simply set the struct contents to empty
*//*-------------------------------------------------------------------*/
int gifti_clear_CoordSystem(giiCoordSystem * p)
{
    if( !p ) return 1;

    p->dataspace = NULL;
    p->xformspace = NULL;
    memset(p->xform,0,sizeof(p->xform));

    return 0;
}

/*----------------------------------------------------------------------
 *! add an empty DataArray struct to the gim->darray list
 *
 *  this both reallocates gim->darray and allocates gim->darray[new]
 *
 *  if num_to_add > 0: add that many elements
 *  if defaults      : init each element with default values
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_add_empty_darray(gifti_image * gim, int num_to_add)
{
    giiDataArray * dptr;
    int            c, ntot, nnew = num_to_add > 0 ? num_to_add : 1;

    if( !gim ) return 1;

    if(G.verb > 3)fprintf(stderr,"++ alloc darray[%d] (+%d)\n",gim->numDA,nnew);

    /* allocate for additional pointers */
    ntot = gim->numDA + nnew;
    gim->darray = (giiDataArray **)realloc(gim->darray,
                                           ntot*sizeof(giiDataArray *));

    if( !gim->darray ) {
        fprintf(stderr,"** failed realloc darray, len %d\n", ntot);
        gim->numDA = 0;
        return 1;
    }

    /* allocate the actual giiDataArray structs, inserted at the end */
    for( c = 0; c < nnew; c++ ) {
        dptr = (giiDataArray *)calloc(1, sizeof(giiDataArray));
        if( !dptr ) {
            fprintf(stderr,"** failed to alloc DA element #%d\n",gim->numDA);
            return 1;   /* leave allocated list as is */
        }
        gim->darray[gim->numDA] = dptr;
        gim->numDA++;

        /* clear any non-value attributes/elements */
        gifti_clear_DataArray(dptr);
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! add the given name=value pair to the nvpairs struct
 *
 *  this allocates memory for the p->name and p->value arrays, along
 *  with duplicating the passed name/value strings
*//*-------------------------------------------------------------------*/
int gifti_add_to_nvpairs(nvpairs * p, const char * name, const char * value)
{
    int index;

    if( !p || !name || !value ) {
        if(G.verb > 1) fprintf(stderr,"** GATN: bad params(%p,%p,%p)\n",
                               (void *)p, (void *)name, (void *)value);
        return 1;
    }

    p->length++;
    p->name = (char **)realloc(p->name, p->length * sizeof(char *));
    p->value = (char **)realloc(p->value, p->length * sizeof(char *));

    if( !p->name || !p->value ) {
        fprintf(stderr,"** GATN: failed to realloc %d pointers\n",p->length);
        return 1;
    } else if ( G.verb > 3 )
        fprintf(stderr,"++ add_nvp [%d]: '%s', '%s'\n", p->length,
                name ? name : "NULL", value ? value : "NULL");

    index = p->length - 1;
    p->name[index] = gifti_strdup(name);
    p->value[index] = gifti_strdup(value);

    if( !p->name[index] || !p->value[index] ) {
        fprintf(stderr,"** GATN: failed to copy pair '%s'='%s'\n",name,value);
        return 1;
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! display the contents of the nvpairs struct
*//*-------------------------------------------------------------------*/
int gifti_disp_nvpairs(const char * mesg, const nvpairs * p)
{
    int c;

    if( mesg ) { fputs(mesg, stderr); fputc(' ', stderr); }

    if( !p ){ fputs("disp: nvpairs = NULL\n", stderr); return 1; }

    fprintf(stderr,"nvpairs struct, len = %d :\n", p->length);

    for(c = 0; c < p->length; c++ )
        fprintf(stderr,"    nvpair: '%s' = '%s'\n",
                G_CHECK_NULL_STR(p->name[c]), G_CHECK_NULL_STR(p->value[c]));
    if( p->length > 0 ) fputc('\n', stderr);

    return 0;
}

/*----------------------------------------------------------------------
 *! display the contents of the LabelTable struct
*//*-------------------------------------------------------------------*/
int gifti_disp_LabelTable(const char * mesg, const giiLabelTable * p)
{
    int c;

    if( mesg ) { fputs(mesg, stderr); fputc(' ', stderr); }

    if( !p ){ fputs("disp: giiLabelTable = NULL\n", stderr); return 1; }

    fprintf(stderr,"giiLabelTable struct, len = %d :\n", p->length);

    for(c = 0; c < p->length; c++ )
        fprintf(stderr,"    index %d, label '%s'\n",
                p->index[c], G_CHECK_NULL_STR(p->label[c]));
    if( p->length > 0 ) fputc('\n', stderr);

    return 0;
}

/*----------------------------------------------------------------------
 *! display the contents of the CoordSystem struct
*//*-------------------------------------------------------------------*/
int gifti_disp_CoordSystem(const char * mesg, const giiCoordSystem * p)
{
    int c1, c2;

    if( mesg ) { fputs(mesg, stderr); fputc(' ', stderr); }

    if( !p ){ fputs("disp: giiCoordSystem = NULL\n", stderr); return 1; }

    fprintf(stderr,"giiCoordSystem struct\n"
                   "    dataspace  = %s\n"
                   "    xformspace = %s\n",
                   G_CHECK_NULL_STR(p->dataspace),
                   G_CHECK_NULL_STR(p->xformspace));
    for( c1 = 0; c1 < 4; c1++ )
    {
        fprintf(stderr,"    xform[%d] :", c1);
        for( c2 = 0; c2 < 4; c2++ )
            fprintf(stderr,"  %f", p->xform[c1][c2]);
        fputc('\n', stderr);
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! display the contents of the DataArray struct
 *
 *  if 'subs' is set, display the contents of the sub-structures
*//*-------------------------------------------------------------------*/
int gifti_disp_DataArray(const char * mesg, const giiDataArray * p, int subs)
{
    fprintf(stderr,"--------------------------------------------------\n");

    if( mesg ) { fputs(mesg, stderr); fputc(' ', stderr); }

    if( !p ){ fputs("disp: giiDataArray = NULL\n", stderr); return 1; }

    fprintf(stderr,"giiDataArray struct\n"
                   "    intent   %4d = %s\n"
                   "    datatype   %2d = %s\n"
                   "    ind_ord    %2d = %s\n"
                   "    num_dim       = %d\n"
                   "    dims          = %d, %d, %d, %d, %d, %d\n"
                   "    encoding   %2d = %s\n"
                   "    endian     %2d = %s\n"
                   "    ext_fname     = %s\n"
                   "    ext_offset    = %lld\n"
                , p->intent,
                gifti_intent_to_string(p->intent),
                p->datatype, gifti_datatype2str(p->datatype),
                p->ind_ord,
                gifti_list_index2string(gifti_index_order_list, p->ind_ord),
                p->num_dim,
                p->dims[0], p->dims[1], p->dims[2],
                p->dims[3], p->dims[4], p->dims[5],
                p->encoding,
                gifti_list_index2string(gifti_encoding_list, p->encoding),
                p->endian,
                gifti_list_index2string(gifti_endian_list, p->endian),
                G_CHECK_NULL_STR(p->ext_fname), p->ext_offset
           );

    if( subs ) gifti_disp_nvpairs("darray->meta", &p->meta);
    if( subs ) gifti_disp_CoordSystem("darray->coordsys", p->coordsys);
                
    fprintf(stderr,"    data       = %s\n"
                   "    nvals      = %u\n"
                   "    nbyper     = %d\n",
                p->data ? "<set>" : "NULL", (unsigned)p->nvals, p->nbyper);

    if( subs ) gifti_disp_nvpairs("darray->ex_atrs", &p->ex_atrs);
    fprintf(stderr,"--------------------------------------------------\n");

    return 0;
}

/*----------------------------------------------------------------------
 *! display the contents of the gifti_image struct
 *
 *  if 'subs' is set, display the contents of the sub-structures, such
 *  as all of the DataArray elements
*//*-------------------------------------------------------------------*/
int gifti_disp_gifti_image(const char * mesg, const gifti_image * p, int subs)
{
    fprintf(stderr,"==================================================\n");

    if( mesg ) { fputs(mesg, stderr); fputc(' ', stderr); }

    if( !p ){ fputs("disp: gifti_image = NULL\n", stderr); return 1; }

    fprintf(stderr,"gifti_image struct\n"
                   "    version    = %s\n"
                   "    numDA      = %d\n",
                   G_CHECK_NULL_STR(p->version), p->numDA);

    if( subs ) {
        char buf[32];
        int  c;

        gifti_disp_nvpairs("gim->meta", &p->meta);
        gifti_disp_LabelTable("gim->labeltable", &p->labeltable);
        for( c = 0; c < p->numDA; c++ ) {
            sprintf(buf, "gim->darray[%d]", c);
            gifti_disp_DataArray(buf, p->darray[c], subs);
        }
    }

    fprintf(stderr,"gifti_image struct\n"
                   "    swapped    = %d\n"
                   "    compressed = %d\n", p->swapped, p->compressed);

    if( subs ) gifti_disp_nvpairs("gim->ex_atrs" , &p->ex_atrs);
    else fprintf(stderr," -- darray totals: %lld MB\n",gifti_gim_DA_size(p,1));

    fprintf(stderr,"==================================================\n");

    return 0;
}

/*----------------------------------------------------------------------
 *! compute the number of bytes summed over all DataArray elements
 *
 *  if in_mb is set, return the (rounded) number of megabytes
 *  (i.e. 17 is 17 MB)
*//*-------------------------------------------------------------------*/
long long gifti_gim_DA_size(const gifti_image * p, int in_mb)
{
    long long bytes = 0;
    int       c;

    if( !p ) return -1;
    if( !p->darray || p->numDA <= 0 ) return 0;

    for( c = 0; c < p->numDA; c++ ) {
        if( ! p->darray[c]->data ) continue;  /* skip anything without data */
        if( p->darray[c]->nvals <= 0 || p->darray[c]->nbyper <= 0 ) {
            fprintf(stderr,"** have data[%d], but nvals = %lld, nbyper = %d\n",
                    c, p->darray[c]->nvals, p->darray[c]->nbyper);
            return 0;
        }
        bytes += p->darray[c]->nvals * p->darray[c]->nbyper;
    }

    if( bytes <= 0LL ) return 0;

    if( in_mb ) bytes = (bytes + (1<<19) ) >> 20;  /* round to nearest MB */

    return bytes;
}

/*----------------------------------------------------------------------
 *! given a datatype, return the corresponding bytes per value and swapsize
 *
 *  nbyper and swapsize are filled only if the pointers are set
*//*-------------------------------------------------------------------*/
int gifti_datatype_sizes(int datatype, int *nbyper, int *swapsize)
{
    int c;

    for( c = sizeof(gifti_type_list) / sizeof(gifti_type_ele) - 1; c > 0; c-- )
        if( datatype == gifti_type_list[c].type ) {
            if( nbyper ) *nbyper = gifti_type_list[c].nbyper;
            if( swapsize ) *swapsize = gifti_type_list[c].swapsize;
            return 0;
        }

    if( G.verb > 0 ) fprintf(stderr,"** GDS with bad datatype %d\n", datatype);
    if( nbyper ) *nbyper = 0;
    if( swapsize ) *swapsize = 0;

    return 1;
}

/*----------------------------------------------------------------------
 *! compute the total number of data values in a DataArray element
*//*-------------------------------------------------------------------*/
long long gifti_darray_nvals(giiDataArray * da)
{
    long long ndim = 1;
    int       c;

    if(!da){ fprintf(stderr,"** GDND, no ptr\n"); return 0; }

    if( ! gifti_valid_num_dim(da->num_dim, 0) ) {
        fprintf(stderr,"** giiDataArray has illegal num_dim = %d\n",
                da->num_dim);
        return 0;
    }

    for( c = 0; c < da->num_dim; c++ ) ndim *= da->dims[c];

    if( ndim <= 0 ) {
        gifti_disp_DataArray("** bad Dim list in ", da, 0);
        return 0;
    }

    return ndim;
}

/*----------------------------------------------------------------------
 *! find giiDataArray element #index of the given intent
*//*-------------------------------------------------------------------*/
giiDataArray * gifti_find_DA(gifti_image * gim, int intent, int index)
{
    int c, nfound;

    if( !gim || !gifti_intent_is_valid(intent) || index < 0 ) {
        fprintf(stderr,"** find_DA: bad inputs (%p, %d, %d)\n",
                (void*)gim, intent, index);
        return NULL;
    }

    if ( !gim-> darray ) return NULL;

    for( c = 0, nfound = 0; c < gim->numDA; c++ )
        if( gim->darray[c] && gim->darray[c]->intent == intent ) {
            if( nfound == index )
                return gim->darray[c];  /* success */
            nfound++;   /* else, increment counter and keep looking */
        }

    return NULL;
}

/*----------------------------------------------------------------------
 *! return an allocated list of giiDataArray pointers of the given intent
 *
 *  'list' should be freed or taken
*//*-------------------------------------------------------------------*/
int gifti_find_DA_list(gifti_image * gim, int intent,
                       giiDataArray *** list, int * len)
{
    int c, nfound;

    if( !gim || !gifti_intent_is_valid(intent) || !list || !len ) {
        fprintf(stderr,"** find_DA: bad inputs (%p, %d, %p, %p)\n",
                (void *)gim, intent, (void *)list, (void *)len);
        return 1;
    }

    if ( !gim->darray ) return 1;

    /* create one big enough to hold everything */
    *len = gim->numDA;
    *list = (giiDataArray **)calloc(*len, sizeof(giiDataArray *));
    if( !*list ) {
        fprintf(stderr,"** find_DA_list: failed to alloc %d ptrs\n",*len);
        *len = 0;
        return 1;
    }

    for( c = 0, nfound = 0; c < gim->numDA; c++ )
        if( gim->darray[c] && gim->darray[c]->intent == intent )
            (*list)[nfound++] = gim->darray[c];

    /* if we didn't find any, nuke list, but do not return an error */
    if( nfound == 0 ) { free(*list); *list = NULL; *len = 0; return 0; }

    /* otherwise, reallocate a smaller list */
    if( nfound < *len ) {
        *len  = nfound;
        *list = (giiDataArray**)realloc(*list,*len*sizeof(giiDataArray*));
        if( !*list ) {
            fprintf(stderr,"** find_DA_list: failed realloc of %d ptrs\n",*len);
            *len = 0;
            return 1;
        }
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! given metadata name, return the corresponding value (or NULL)
 *
 *  no allocation is done here
*//*-------------------------------------------------------------------*/
char * gifti_get_meta_value(nvpairs * nvp, char * name)
{
    int c;

    if( !nvp || !name ) {
        if( G.verb > 3 )
          fprintf(stderr,"** get_meta_value: NULL input (%p, %p)\n",
                (void*)nvp, name);
        return NULL;
    }

    if( G.verb > 5 )
        fprintf(stderr,"-- G_get_meta_value: looking for name = '%s'\n", name);

    if ( !nvp->name || !nvp->value || nvp->length <= 0 ) {
        if( G.verb > 3 )
            fprintf(stderr,"-- G_get_meta_value: no name/value array\n");
        return NULL;
    }

    for( c = 0; c < nvp->length; c++ )
        if( !strcmp(nvp->name[c], name) )
            break;  /* found */

    if( c >= nvp->length ) return NULL;

    if( G.verb > 3 )
        fprintf(stderr,"++ found meta '%s'='%s'\n",nvp->name[c],nvp->value[c]);

    return nvp->value[c];
}

/*----------------------------------------------------------------------
 *! return the number of 'rows' and 'columns' of a DataArray element
 *
 *  define rows to be the number of nodes, which should be the slowest
 *  changing element, depending on the index order (kuru kuru pa)
*//*-------------------------------------------------------------------*/
int gifti_DA_rows_cols(giiDataArray * da, long long * rows, long long * cols)
{
    *rows = da->dims[0];  /* init */
    *cols = 1;
                                                                                
    if( da->num_dim == 1 ) return 0;  /* use default */
                                                                                
    if( da->ind_ord == GIFTI_IND_ORD_ROW_MAJOR ) {
        /* treat Dim[0] as nodes (they change most slowly) */
        *rows = da->dims[0];
        *cols = (*rows) ? da->nvals / *rows : 1;    /* be safe */
    } else {
        if( ! gifti_valid_num_dim(da->num_dim, 1) ){
            fprintf(stderr,"** cannot assign DA_rows_cols");
            return 1;
        }

        *rows = da->dims[da->num_dim-1];  /* take highest index */
        *cols = (*rows > 0) ? da->nvals / *rows : 1;
    }
                                                                                
    return 0;
}

/*----------------------------------------------------------------------
 *! print the gifti library history string
*//*-------------------------------------------------------------------*/
void gifti_disp_lib_hist(void)
{
   int c, len = sizeof(gifti_history)/sizeof(char *);
   for( c = 0; c < len; c++ )
       fputs(gifti_history[c], stdout);
}

/*----------------------------------------------------------------------
 *! print the gifti library version string
*//*-------------------------------------------------------------------*/
void gifti_disp_lib_version(void)
{
    printf("%s, compiled %s\n", gifti_version, __DATE__);
}

/*----------------------------------------------------------------------
 *! return the gifti library version string
*//*-------------------------------------------------------------------*/
char * gifticlib_version(void)
{
    return gifti_version;
}

/*----------------------------------------------------------------------
 *! print the gifti DTD URL
*//*-------------------------------------------------------------------*/
void gifti_disp_dtd_url(void)
{
    printf("The GIFTI Document Type Definition (DTD) is at:\n    %s\n",
           GIFTI_XML_DTD_SOURCE);
}

/*----------------------------------------------------------------------
 *! display data in hexidecimal, on one line
 *
 *  if mesg is set, print the message first
 *  if fp is not set, print to stdout
*//*-------------------------------------------------------------------*/
int gifti_disp_hex_data(const char *mesg, const void *data, int len, FILE *fp)
{
    const char * dp = (const char *)data;
    FILE       * stream;
    int          c;
    
    stream = fp ? fp : stdout;

    if( !data || len < 1 ) return -1;

    if( mesg ) fputs(mesg, stream);

    for( c = 0; c < len; c++ )
        fprintf(stream, " %02x", dp[c]);

    return 0;
}

/*----------------------------------------------------------------------
 *! swap sets of 2-byte values
*//*-------------------------------------------------------------------*/
int gifti_swap_2bytes(void * data, long long nsets)
{
    char    * cp1 = (char *)data, * cp2;
    char      tval;
    long long c;
    
    for( c = 0; c < nsets; c++ ) {
        cp2 = cp1 + 1;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp1 += 2;
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! swap sets of 4-byte values
*//*-------------------------------------------------------------------*/
int gifti_swap_4bytes(void * data, long long nsets)
{
    char    * cp0 = (char *)data, * cp1, * cp2;
    char      tval;
    long long c;
    
    for( c = 0; c < nsets; c++ ) {
        cp1 = cp0;
        cp2 = cp0 + 3;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp1++;  cp2--;
        tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
        cp0 += 4;
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! swap sets of N-byte values
 *  
 *  if N < 2         : just return
 *  if N = 2 or N = 4: call explicit function for that size (speed)
 *  else             : swap in a loop
*//*-------------------------------------------------------------------*/
int gifti_swap_Nbytes(void * data, long long nsets, int swapsize)
{
    char    * cp0, * cp1, * cp2;
    char      tval;
    long long c;
    int       offset;      /* swapsize - 1, for speed */

    if( ! data || nsets < 0 || swapsize < 0 ) {
        fprintf(stderr,"** swap_Nbytes: bad params (%p,%lld,%d)\n",
                (void *)data, nsets, swapsize);
        return 1;
    }

    if     ( swapsize  < 2 ) return 0;  /* nothing to do */
    else if( swapsize == 2 ) return gifti_swap_2bytes(data, nsets);
    else if( swapsize == 4 ) return gifti_swap_4bytes(data, nsets);
    
    /* peform a swap */
    cp0 = (char *)data;
    offset = swapsize-1;  /* just for speed */

    for( c = 0; c < nsets; c++ ) {
        cp1 = cp0;
        cp2 = cp0 + offset;
        while( cp2 > cp1 ) {
            tval = *cp1;  *cp1 = *cp2;  *cp2 = tval;
            cp1++;  cp2--;
        }
        cp0 += swapsize;
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! return the current CPU endian: GIFTI_ENDIAN_BIG or _LITTLE
*//*-------------------------------------------------------------------*/
int gifti_get_this_endian(void)
{
   int    one = 1;
   char * cp = (char *)&one;

   if( *cp ) return GIFTI_ENDIAN_LITTLE;

   return GIFTI_ENDIAN_BIG;
}

/*----------------------------------------------------------------------
 *! if endian does not match that of this CPU, swap the data bytes
 *
 *  return whether bytes were swapped
*//*-------------------------------------------------------------------*/
int gifti_check_swap(void * data, int endian, long long nsets, int swapsize)
{
    if( !data || nsets < 0 || swapsize < 0 ) {
        fprintf(stderr,"** check_swap: bad params (%p,%lld, %d)\n",
                (void *)data, nsets, swapsize);
        return 0;
    } else if ( endian <= GIFTI_ENDIAN_UNDEF || endian > GIFTI_ENDIAN_MAX ) {
        fprintf(stderr, "** check_swap: invalid endian %d\n", endian);
        return 0;
    }

    /* if endian is the same as this one, just return */
    if( endian == gifti_get_this_endian() ) {
        if( G.verb > 2 )
            fprintf(stderr,"-- darray no swap needed : %lld sets of %d bytes\n",
                    nsets, swapsize);
        return 0;
    }

    if( G.verb > 2 )
        fprintf(stderr,"++ darray swap: %lld sets of %d bytes\n",
                nsets, swapsize);

    /* do the swap */
    (void)gifti_swap_Nbytes(data, nsets, swapsize);

    return 1;
}

/*---------------------------------------------------------------------*/
/*! Given a NIFTI_INTENT string, such as "NIFTI_INTENT_NODE_INDEX",
 *  return the corresponding integral intent code.  The intent code is
 *  the macro value defined in nifti1.h.
 *
 *  return 0 on failure (NIFTI_INTENT_NONE)
*//*-------------------------------------------------------------------*/
int gifti_intent_from_string( const char * name )
{
    int tablen = sizeof(gifti_intent_list)/sizeof(gifti_intent_ele);
    int c;

    if( !name ) return 0;

    for( c = tablen-1; c > 0; c-- )
        if( !strcmp(name, gifti_intent_list[c].name) )
            break;

    return gifti_intent_list[c].code;
}


/*---------------------------------------------------------------------*/
/*! Given a NIFTI_TYPE value, such as NIFTI_TYPE_INT16, return the
 *  corresponding macro label as a string.  The dtype code is the
 *  macro value defined in nifti1.h.
*//*-------------------------------------------------------------------*/
char * gifti_intent_to_string( int code )
{
    int tablen = sizeof(gifti_intent_list)/sizeof(gifti_intent_ele);
    int c;

    for( c = tablen-1; c > 0; c-- )
        if( gifti_intent_list[c].code == code )
            break;

    return gifti_intent_list[c].name;
}


/*---------------------------------------------------------------------*/
/*! Return whether the given code is a valid NIFTI_INTENT code.
*//*-------------------------------------------------------------------*/
int gifti_intent_is_valid( int code )
{
    int tablen = sizeof(gifti_intent_list)/sizeof(gifti_intent_ele);
    int c;

    for( c = tablen-1; c > 0; c-- )
        if( gifti_intent_list[c].code == code )
            break;

    return( c > 0 );
}


/*---------------------------------------------------------------------*/
/*! duplicate the given string
*//*-------------------------------------------------------------------*/
char * gifti_strdup( const char * src )
{
    char * newstr;
    int    len;

    if( !src ) return NULL;  /* easy case */

    len = strlen(src) + 1;

    newstr = (char *)malloc(len * sizeof(char));
    if( !newstr ) {
        fprintf(stderr,"** failed gifti_strdup, len = %d\n", len);
        return NULL;
    }

    strcpy(newstr, src);

    return newstr;
}

/*---------------------------------------------------------------------*/
/*! duplicate the given giiDataArray struct, optionally including data
 *
 *  Allocate for a new struct, and fill it so that the contents are
 *  identical (sans pointers) to orig.  Sub-structure arrays will need
 *  to be allocated, also.
 *
 *  If get_data is not set, gnew->data will be left as NULL.
 *  
 *  return the address of the newly allocated structure
*//*-------------------------------------------------------------------*/
giiDataArray * gifti_copy_DataArray(const giiDataArray * orig, int get_data)
{
    giiDataArray * gnew;

    if( ! orig ){ fprintf(stderr,"** copy_DA: input is NULL\n"); return NULL; }

    if(G.verb > 5) fprintf(stderr,"++ copying giiDataArray...\n");

    gnew = (giiDataArray *)calloc(1, sizeof(giiDataArray));
    if(!gnew){fprintf(stderr,"** copy_DA, failed to alloc DA\n"); return NULL;}

    /* cheat - start by copying the entire struct contents */
    *gnew = *orig;

    /* copy any pointer data or structures */
    gnew->ext_fname = gifti_strdup(orig->ext_fname);
    gifti_copy_nvpairs(&gnew->meta, &orig->meta);
    gnew->coordsys = gifti_copy_CoordSystem(orig->coordsys);

    /* maybe the needy user wants data, too */
    if(orig->data && get_data) {
        if(G.verb > 5) fprintf(stderr,"++ copy_DA, adding data\n");
        gnew->data = malloc(gnew->nvals * gnew->nbyper);
        if(!gnew->data)       /* continue? */
            fprintf(stderr,"** copy DA, failed to alloc %lld bytes for data\n",
                    gnew->nvals * gnew->nbyper);
        memcpy(gnew->data, orig->data, gnew->nvals * gnew->nbyper);
    } else
        gnew->data = NULL;

    /* last and certainly least, ex_atrs */
    gifti_copy_nvpairs(&gnew->ex_atrs, &orig->ex_atrs);

    return gnew;
}

/*---------------------------------------------------------------------*/
/*! dupliate the giiCoordSystem struct (passing NULL is okay)
*//*-------------------------------------------------------------------*/
giiCoordSystem * gifti_copy_CoordSystem(const giiCoordSystem * src)
{
    giiCoordSystem * csnew;
    int              r, c;

    if( !src ) return NULL;     /* this may be okay */

    if( G.verb > 6 ) fprintf(stderr,"++ copy_CS\n");

    csnew = (giiCoordSystem *)malloc(sizeof(giiCoordSystem));
    if( !csnew ){ fprintf(stderr,"** copy_CS: failed alloc\n"); return NULL; }

    csnew->dataspace  = gifti_strdup(src->dataspace);
    csnew->xformspace = gifti_strdup(src->xformspace);

    for( r = 0; r < 4; r++ )
        for( c = 0; c < 4; c++ )
            csnew->xform[r][c] = src->xform[r][c];

    return csnew;
}

/*---------------------------------------------------------------------*/
/*! dupliate the contents of the giiLabelTable struct
*//*-------------------------------------------------------------------*/
int gifti_copy_LabelTable(giiLabelTable * dest, const giiLabelTable * src)
{
    int c;

    if( !src || !dest ) {
        fprintf(stderr,"** copy_LabelTable: bad params (%p,%p)\n",
                (void *)src, (void *)dest);
        return 1;
    }

    if( G.verb > 6 ) fprintf(stderr,"++ copy_LT\n");

    /* quick case: empty table */
    if( src->length <= 0 ) {
        dest->length = 0;
        dest->index = NULL;
        dest->label = NULL;
        return 0;
    }

    /* otherwise, allocate and copy lists */
    dest->length = src->length;

    dest->index = (int *)malloc(dest->length * sizeof(int));
    dest->label = (char **)malloc(dest->length * sizeof(char *));

    /* check for failure */
    if( !dest->index || !dest->label ) {
        fprintf(stderr,"** failed to dup label arrays of length %d\n",
                dest->length);
        gifti_free_LabelTable(dest);
        return 1;
    }

    /* copy indices */
    for( c = 0; c < dest->length; c++ )
        dest->index[c] = src->index[c];

    /* copy labels */
    for( c = 0; c < dest->length; c++ )
        dest->label[c] = gifti_strdup(src->label[c]);

    return 0;
}

/*---------------------------------------------------------------------*/
/*! dupliate the contents of one nvpairs structure into an empty one
 *
 *  return 0 on success
*//*-------------------------------------------------------------------*/
int gifti_copy_nvpairs(nvpairs * dest, const nvpairs * src)
{
    if( !dest || !src ){
        fprintf(stderr,"** copy_NVP, bad params (%p,%p)\n",
                (void*)dest,(void*)src);
        return 1;
    }

    if( G.verb > 6 ) fprintf(stderr,"++ copy_nvp, length %d\n", src->length);

    /* check for a simple case */
    if( src->length <= 0 || !src->name || !src->value ) {
        dest->length = 0;
        dest->name = dest->value = NULL;
        return 0;
    }

    /* else, copy the lists */
    dest->length = src->length;
    dest->name   = gifti_copy_char_list(src->name, src->length);
    dest->value  = gifti_copy_char_list(src->value, src->length);

    return 0;
}

/*---------------------------------------------------------------------*/
/*! dupliate the list of strings
*//*-------------------------------------------------------------------*/
char ** gifti_copy_char_list(char ** list, int len)
{
    char ** newlist = NULL;
    int     c;

    if( !list || len <= 0 ) return NULL;

    newlist = (char **)malloc(len * sizeof(char *));
    if( !newlist ) {
        fprintf(stderr,"** copy_char_list fails for %d pointers\n", len);
        return NULL;
    }

    for( c = 0; c < len; c++)                   /* big finish */
        newlist[c] = gifti_strdup(list[c]);

    return newlist;
}


/*---------------------------------------------------------------------*/
/*! print raw data (nvals of type 'type') to the given file stream
 *
 *  possibly write a trailing newline
*//*-------------------------------------------------------------------*/
int gifti_disp_raw_data(const void * data, int type, int nvals, int newline,
                        FILE * stream)
{
    FILE * fp = stream ? stream : stdout;
    char * dp, fbuf[64];
    int    c, size;

    gifti_datatype_sizes(type, &size, NULL);   /* get nbyper */
    if( size == 0 ) {
        fprintf(stderr,"** GDRD: cannot print with size 0, type %d\n", type);
        return 1;
    }

    for( c = 0, dp = (char *)data; c < nvals; c++, dp += size ) {
        switch( type ) {
            case NIFTI_TYPE_INT8:
                fprintf(fp, "%d", *(char *)dp);
                break;
            case NIFTI_TYPE_INT16:
                fprintf(fp, "%d", *(short *)dp);
                break;
            case NIFTI_TYPE_INT32:
                fprintf(fp, "%d", *(int *)dp);
                break;
            case NIFTI_TYPE_INT64:
                fprintf(fp, "%lld", *(long long *)dp);
                break;
            case NIFTI_TYPE_UINT8:
                fprintf(fp, "%u", *(unsigned char *)dp);
                break;
            case NIFTI_TYPE_UINT16:
                fprintf(fp, "%u", *(unsigned short *)dp);
                break;
            case NIFTI_TYPE_UINT32:
                fprintf(fp, "%u", *(unsigned int *)dp);
                break;
            case NIFTI_TYPE_UINT64:
                fprintf(fp, "%llu", *(unsigned long long *)dp);
                break;
            case NIFTI_TYPE_FLOAT32:
                sprintf(fbuf,"%f", *(float *)dp);
                gifti_clear_float_zeros(fbuf);
                fprintf(fp, "%s", fbuf);
                break;
            case NIFTI_TYPE_FLOAT64:
                sprintf(fbuf,"%f", *(double *)dp);
                gifti_clear_float_zeros(fbuf);
                fprintf(fp, "%s", fbuf);
                break;
            default:
                fprintf(stderr,"** Gdisp_raw_data: unknown type %d\n", type);
                return 1;
        }
        if( c < nvals - 1 ) fputc(' ', fp);
    }

    if ( newline ) fputc('\n', fp);

    return 0;
}

/*----------------------------------------------------------------------
 *! remove trailing zeros from string of printed float
 *  return  1 if something was cleared
 *          0 if not
*//*-------------------------------------------------------------------*/
int gifti_clear_float_zeros( char * str )
{
   char * dp, * valp;
   int    len;

   if( !str || !*str ) return 0;  /* that was easy    */

   dp  = strchr(str, '.');        /* find '.'         */
   if( !dp ) return 0;            /* no '.', easy too */

   len = strlen(dp);

   /* never clear what is just to the right of '.' */
   for( valp = dp+len-1; (valp > dp+1) && (*valp==' ' || *valp=='0'); valp-- )
       *valp = '\0';     /* clear, so we don't worry about break conditions */

   if( valp < dp + len - 1 ) return 1;

   return 0;
}

/*----------------------------------------------------------------------
 *! set all DataArray attributes of the given name to the given value
*//*-------------------------------------------------------------------*/
int gifti_set_atr_in_DAs(gifti_image *gim, const char *name, const char *value,
                         const int * dalist, int len)
{
    int c, ind;

    if( !gim || !name || !value ) {
        fprintf(stderr,"** set_DA_atrs: bad params (%p,%p,%p)\n",
                (void *)gim, (void *)name, (void *)value);
        return 1;
    }

    if( !gim->darray ) return 0;    /* just leave */

    if( dalist && len > 0 ) {       /* set atrs for those DA in dalist */
        if( ! gifti_valid_int_list(dalist, len, 0, gim->numDA-1, 1) )
            return 1;

        /* good to go */
        for( c = 0; c < len; c++ ) {
            ind = dalist[c];    /* for clarity */
            if( ! gim->darray[ind] ) continue;  /* trioanoid */
            if( gifti_str2attr_darray(gim->darray[ind], name, value) ) {
                if( G.verb > 1 )
                    fprintf(stderr,"** bad DA attr '%s'='%s'\n",name,value);
                return 1;
            }
        }

        if( G.verb > 2 )
            fprintf(stderr,"++ set atrs in %d DAs, '%s'='%s'\n",len,name,value);

        return 0;
    }

    /* else continue, do all of them */

    for( c = 0; c < gim->numDA; c++ ) {
        if( ! gim->darray[c] ) continue;  /* trioanoid */
        if( gifti_str2attr_darray(gim->darray[c], name, value) ) {
            if(G.verb>1)fprintf(stderr,"** bad DA attr '%s'='%s'\n",name,value);
            return 1;
        }
    }

    if( G.verb > 4 )
        fprintf(stderr,"++ set attr in all DAs, '%s'='%s'\n", name, value);

    return 0;
}

/*----------------------------------------------------------------------
 *! set MetaData name/value pairs in all DAs in list (or all in gim)
*//*-------------------------------------------------------------------*/
int gifti_set_DA_meta( gifti_image *gim, const char *name, const char *value,
                       const int * dalist, int len, int replace )
{
    int c, ind;

    if( !gim || !name || !value ) {
        fprintf(stderr,"** set_DA_meta: bad params (%p,%p,%p)\n",
                (void *)gim, (void *)name, (void *)value);
        return 1;
    }

    if( !gim->darray ) return 0;    /* just leave */

    if( dalist && len > 0 ) {       /* set meta for those DA in dalist */
        if( ! gifti_valid_int_list(dalist, len, 0, gim->numDA-1, 1) )
            return 1;

        /* good to go */
        for( c = 0; c < len; c++ ) {
            ind = dalist[c];    /* for clarity */
            if( ! gim->darray[ind] ) continue;  /* trioanoid */
            if( gifti_add_to_meta(&gim->darray[ind]->meta,name,value,replace) )
                return 1;
        }

        if( G.verb > 2 )
            fprintf(stderr,"++ set meta in %d DAs, '%s'='%s'\n",len,name,value);

        return 0;
    }

    /* else continue, do all of them */

    for( c = 0; c < gim->numDA; c++ ) {
        if( ! gim->darray[c] ) continue;  /* trioanoid */
        if( gifti_add_to_meta(&gim->darray[c]->meta,name,value,replace) )
            return 1;
    }

    if( G.verb > 4 )
        fprintf(stderr,"++ set MetaData in all DAs, '%s'='%s'\n", name, value);

    return 0;
}

/*----------------------------------------------------------------------
 *! allocate and return a newly created gifti_image_struct
 *
 *  if numDA > 0, allocate and initialize gim->darray
 *  if intent is a valid NIFTI_INTENT code, set it in all DataArray elements
 *  if dtype is a valid NIFTI_TYPE, set it in all DA elements
 *  if dims is set (MUST be of length 6), set the dims in all DA elements
 *  if alloc_data, allocate zero-filled data in each DA element
 *
 *  note that if numDA <= 0, the function returns an empty gifti_image
*//*-------------------------------------------------------------------*/
gifti_image * gifti_create_image( int numDA, int intent, int dtype, int ndim,
                                  const int * dims, int alloc_data )
{
    gifti_image * gim;
    int           c, errs = 0;

    if( G.verb > 1 ) {  /* maybe start with some chit-chat */
        fprintf(stderr,"++ creating gifti_image with %d DA elements\n", numDA);
        if( G.verb > 2 ) {
            fprintf(stderr,"     intent[%d] = %s, dtype[%d] = %s,\n"
                           "     alloc_data = %d, ndim = %d, dims: ",
                    intent, gifti_intent_to_string(intent),
                    dtype,  gifti_datatype2str(dtype), alloc_data, ndim);
            if( !dims ) fprintf(stderr,"<empty>\n");
            else gifti_disp_raw_data(dims, DT_INT32, GIFTI_DARRAY_DIM_LEN,
                                     1, stderr);
        }
    }

    /* basic step - create empty image (with a version string) */
    gim = (gifti_image *)calloc(1, sizeof(gifti_image));
    if(!gim){ fprintf(stderr,"** failed to alloc gifti_image\n"); return NULL; }

    gifti_clear_gifti_image(gim);
    gim->version = gifti_strdup(GIFTI_XML_VERSION);

    if( numDA <= 0 ) return gim;        /* done */

    /* apply numDA, which is incremented in add_empty_darray() */
    gim->numDA = 0;
    if( gifti_add_empty_darray(gim, numDA) ) {  /* then cannot continue */
        gifti_free_image(gim);
        return NULL;
    }

    /* init all to defaults */
    for( c = 0; c < gim->numDA; c++ )
        errs += gifti_set_DA_defaults(gim->darray[c]);

    /* and fill in any other pieces */

    if( gifti_intent_is_valid(intent) )
        errs += gifti_set_atr_in_DAs(gim,"Intent",
                                     gifti_intent_to_string(intent), NULL, 0);

    if( gifti_valid_datatype(dtype, 1) ) {
        errs += gifti_set_atr_in_DAs(gim,"DataType", gifti_datatype2str(dtype),
                                     NULL, 0);
        errs += gifti_update_nbyper(gim);
    }

    if( dims && ndim >= 0 ) errs += gifti_set_dims_all_DA(gim, ndim, dims);

    /* don't try this if there are errors */
    if( !errs && alloc_data ) errs += gifti_alloc_DA_data(gim, NULL, 0);

    if( errs ) {  /* then fail out */
        gifti_free_image(gim);      
        return NULL;
    }

    return gim;
}

/*----------------------------------------------------------------------
 *! allocate nvals*nbyper bytes of (zero-filled) data in each DataArray
 *  (in dalist or simply in gim)
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_alloc_DA_data(gifti_image * gim, const int * dalist, int len)
{
    giiDataArray * da;
    long long      nbytes, ntot = 0;
    int            c, index, nset = 0, use_list, numDA;

    if( !gim || !gim->darray || gim->numDA <= 0 ) return 0;

    /* decide whether to use dalist or allocate all data */
    use_list = gifti_valid_int_list(dalist, len, 0, gim->numDA-1, 0);

    if( use_list && G.verb > 2 )
        fprintf(stderr,"++ allocating data for %s\n",
                use_list ? "DA in list" : "all DAs");

    if( use_list && DA_data_exists(gim, dalist, len) ) {
        fprintf(stderr,"** data already exists for some DAs in list\n");
        return 1;
    }

    numDA = use_list ? len : gim->numDA;
    for( c = 0; c < numDA; c++ ) {
        index = use_list ? dalist[c] : c;  /* choose appropriate DA index */
        da = gim->darray[index];           /* set convenient pointer */

        if( ! da ) continue;

        /* if dimensions do not make sense, fail out */
        if( ! gifti_valid_dims(da, G.verb > 0) ) return 1;

        if( da->nvals < 0 || da->nbyper < 0 ) {
            fprintf(stderr,"** bad nvals, nbyper in DA[%d]\n",index);
            return 1;
        }
        nbytes = da->nvals * da->nbyper;
        if( nbytes <= 0 ) continue;  /* skip empty data */

        ntot += nbytes;          /* compute total bytes */
        nset++;

        da->data = calloc(nbytes, sizeof(char));
        if( !da->data ) {
            fprintf(stderr,"** gifti_alloc_DA_data: failed on DA %d of %d\n"
                           "     %lld bytes (%lld total)\n",
                           index, numDA, nbytes, ntot);
            return 1;
        }
    }

    if( G.verb > 3)
        fprintf(stderr,"++ alloc'd %lld bytes in %d DA elements\n", ntot, nset);

    return 0;
}

/*----------------------------------------------------------------------
 *! set num_dims, dims and nvals in every DataArray element
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_set_dims_all_DA(gifti_image * gim, int ndim, const int * dims)
{
    long long nvals;
    int       c, d, nset = 0;

    if(!gim || ndim < 0 || ndim > GIFTI_DARRAY_DIM_LEN || !dims) {
        fprintf(stderr,"** SDA_DA: bad params (%p, %d, %p)\n",
                (void *)gim, ndim, (void *)dims);
        return 1;
    }

    if( !gim->darray || gim->numDA == 0 ) return 0;

    /* first compute nvals */
    for( d = 0, nvals = 1; d < ndim; d++) nvals *= dims[d];
    if( nvals <= 0 ) {
        fprintf(stderr,"** GSDA_DA: malformed dims[%d]: ", ndim);
        gifti_disp_raw_data(dims, DT_INT32, GIFTI_DARRAY_DIM_LEN, 1, stderr);
        return 1;
    }
    if( ndim == 0 ) nvals = 0;

    /* insert num_dim and fill dims (pad with 0s) */
    for( c = 0; c < gim->numDA; c++ ) {
        if( !gim->darray[c] ) continue;  /* paranoid */
        gim->darray[c]->num_dim = ndim;
        for( d = 0; d < ndim; d++ )
            gim->darray[c]->dims[d] = dims[d];
        for( /* continue */ ; d < GIFTI_DARRAY_DIM_LEN; d++ )
            gim->darray[c]->dims[d] = 0;
        gim->darray[c]->nvals = nvals;
        nset++;
    }

    if(G.verb > 3) {
        fprintf(stderr,"++ set dims in %d of %d DA elements to: ",
                nset, gim->numDA);
        gifti_disp_raw_data(dims, DT_INT32, GIFTI_DARRAY_DIM_LEN, 1, stderr);
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! update nbyper/swapsize for all DataArray elements
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_update_nbyper(gifti_image * gim)
{
    giiDataArray * da;
    int            c, errs = 0;

    if( !gim ) return 1;

    if( !gim->darray || gim->numDA == 0 ) return 0;

    for( c = 0; c < gim->numDA; c++ ){
        da = gim->darray[c];    /* just to look cleaner */
        if( !da ) continue;
        errs += gifti_datatype_sizes(da->datatype, &da->nbyper, NULL);
    }

    return errs;
}

/*----------------------------------------------------------------------
 *! fill DataArray element with default values
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_set_DA_defaults(giiDataArray * da)
{
    int c;

    if(!da) { fprintf(stderr,"** NULL in set_DA_defaults\n"); return 1; }

    if( G.verb > 6 ) fprintf(stderr,"-- setting DA defaults\n");

    gifti_clear_DataArray(da);                  /* start with empty struct */

    /* and fill with any non-NULL, non-zero values */

    da->intent = NIFTI_INTENT_NONE;
    da->datatype = NIFTI_TYPE_FLOAT32;
    da->ind_ord = GIFTI_IND_ORD_ROW_MAJOR;
    da->num_dim = 1;                            /* one value per node */

    for( c = 0; c < GIFTI_DARRAY_DIM_LEN; c++ ) da->dims[c] = 0;

    da->encoding = GIFTI_ENCODING_B64BIN;       /* zlib may not be available */
    da->endian = gifti_get_this_endian();
    da->ext_offset = 0;

    da->nvals = 0;
    da->nbyper = 0;
    gifti_datatype_sizes(da->datatype, &da->nbyper, NULL);

    return 0;
}

/*----------------------------------------------------------------------
 *! clear the DataArray element
*//*-------------------------------------------------------------------*/
int gifti_clear_DataArray(giiDataArray * da)
{
    if(!da) { fprintf(stderr,"** NULL in clear_DataArray\n"); return 1; }

    if( G.verb > 5 ) fprintf(stderr,"-- clearing DataArray\n");

    memset(da, 0, sizeof(giiDataArray));

    da->ext_fname = NULL;
    gifti_clear_nvpairs(&da->meta);
    da->coordsys = NULL;
    da->data = NULL;
    gifti_clear_nvpairs(&da->ex_atrs);

    return 0;
}

/*----------------------------------------------------------------------
 *! simply clear all contents of the passed gifti_image
 *  (being explicit with pointers)
 *
 *  return 0 on success
 *         1 on error
*//*-------------------------------------------------------------------*/
int gifti_clear_gifti_image(gifti_image * gim)
{
    if(!gim) { fprintf(stderr,"** NULL in clear_gifti_image\n"); return 1; }

    if( G.verb > 5 ) fprintf(stderr,"-- clearing gifti_image\n");

    /* set the version and clear all pointers */
    memset(gim, 0, sizeof(gim));

    gim->version = NULL;
    gifti_clear_nvpairs(&gim->meta);
    gifti_clear_LabelTable(&gim->labeltable);
    gim->darray = NULL;
    gifti_clear_nvpairs(&gim->ex_atrs);

    return 0;
}

/*----------------------------------------------------------------------
 *! read a dataset, just for numDA
 *
 *  may write faster gxml function for this, if it seems important
*//*-------------------------------------------------------------------*/
int gifti_read_dset_numDA(const char * fname)
{
    gifti_image * gim;
    int           numDA;

    if( !fname ) {
        fprintf(stderr,"** NULL to gifti_read_dset_numDA\n");
        return -1;
    }

    if( G.verb > 2 ) fprintf(stderr,"++ read dset numDA, file '%s'\n",fname);

    gim = gifti_read_da_list(fname, 0, NULL, 0);

    if( !gim ) return -1;       /* errors already printed */

    numDA = gim->numDA;

    if(G.verb > 1)
        fprintf(stderr,"++ read dset numDA, file '%s', numDA = %d\n",
                fname, numDA);

    gifti_free_image(gim);      /* lose dataset and return */

    return numDA;
}

/*----------------------------------------------------------------------
 *! return whether the list values are from min to max
*//*-------------------------------------------------------------------*/
int gifti_valid_int_list(const int *list, int len, int min, int max, int whine)
{
    int c;

    if( !list || len <= 0 ) return 0;

    for( c = 0; c < len; c++ )
        if( list[c] < min || list[c] > max ) {
            if( whine )
                fprintf(stderr,"** bad list index [%d] = %d, not in [%d,%d]\n",
                    c, list[c], min, max);
            return 0;
        }

    return 1;
}

/* return whether any DAs in list have data */
static int DA_data_exists(gifti_image * gim, const int * dalist, int len)
{
    int length, uselist = 0;
    int c, ind;

    if( !dalist || len <= 0 ) { /* then scan all DA elements */
        length = gim->numDA;
        if( length <= 0 ) return 0;
    } else if( !gifti_valid_int_list(dalist, len, 0, gim->numDA-1, 1) ) {
        return 0;
    } else {
        uselist = 1;
        length = len;
    }

    for( c = 0; c < length; c++ ) {
        ind = uselist ? dalist[c] : c;
        if( gim->darray[ind] && gim->darray[ind]->data )
            return 1;
    }

    return 0;
}

/*----------------------------------------------------------------------
 *! add the name=value pair to the MetaData lists
 *
 *  if the name already exists, fail, unless 'replace' is set
 *
 *  return 0 on success, 1 on error
*//*-------------------------------------------------------------------*/
int gifti_add_to_meta( giiMetaData * md, const char * name, const char * value,
                       int replace )
{
    int c;

    if( !md || !name || !value ) return 1;

    if( G.verb > 5 )
        fprintf(stderr,"++ GA2M: name '%s', value '%s', replace = %d\n",
                       name, value, replace);

    /* see if 'name' is already here */
    for( c = 0; c < md->length; c++ )
    {
        if( !md->name[c] && G.verb > 2 ) {
            fprintf(stderr,"** G MD[%d]: no name to check for replacement\n",c);
            continue;
        }

        if( !strcmp(md->name[c], name) ) {      /* a match, apply and return */
            if( !md->value[c] && G.verb > 2 ) {
                fprintf(stderr,"** G MD[%d]: no value to replace\n",c);
                md->value[c] = gifti_strdup(value);
                return 0;
            }

            if( replace ) {
                if( G.verb > 5 ) fprintf(stderr,"   (add via REPLACE)\n");
                if( md->value[c] ) free(md->value[c]);
                md->value[c] = gifti_strdup(value);
                return 0;
            } else {
                fprintf(stderr,"** G_add_to_meta: name '%s', already exists\n",
                               name);
                return 1;
            }
        }
    }

    /* name is new, so just add it */

    if( G.verb > 5 ) fprintf(stderr,"   (adding new entry)\n");

    md->length++;
    md->name = (char **)realloc(md->name, md->length * sizeof(char *));
    md->value = (char **)realloc(md->value, md->length * sizeof(char *));

    if( !md->name || !md->value ) {
        fprintf(stderr,"** GA2M:failed to realloc %d MD pointers\n",md->length);
        md->length = 0;
        return 1;
    }

    md->name[md->length-1] = gifti_strdup(name);
    md->value[md->length-1] = gifti_strdup(value);

    if( ! md->name[md->length-1] || ! md->value[md->length-1] )
        return 1;

    return 0;
}

/*----------------------------------------------------------------------
 *! check for validity of the gifti_image (including sub-structures)
 *
 *  if whine is set, complain about any errors
 *
 *  return 1 if valid, 0 otherwise
*//*-------------------------------------------------------------------*/
int gifti_valid_gifti_image( gifti_image * gim, int whine )
{
    int c, errs = 0;

    if( !gim ) {
        if(whine) fprintf(stderr,"** invalid: gifti_image ptr is NULL\n");
        return 0;
    }

    if( G.verb > 3 ) fprintf(stderr,"-- checking for valid gifti_image...\n");

    if( gim->numDA < 0 ) {
        if(whine) fprintf(stderr,"** invalid: numDA = %d\n", gim->numDA);
        errs ++;
    }

    if( !gim->version || strcmp(gim->version, GIFTI_XML_VERSION) ) {
        if(whine) fprintf(stderr, "** invalid: gim version = %s\n",
                          G_CHECK_NULL_STR(gim->version));
        errs++;
    }

    if( ! gifti_valid_nvpairs(&gim->meta, whine) ) errs ++;

    if( ! gifti_valid_LabelTable(&gim->labeltable, whine) ) errs ++;

    for( c = 0; c < gim->numDA; c++ ) {
        if( G.verb > 5 ) fprintf(stderr,"-- checking DA[%d]\n", c);
        if( ! gifti_valid_DataArray(gim->darray[c], whine) ) {
            if( G.verb > 3 ) fprintf(stderr,"-- DA[%d] has errors\n",c);
            errs++;
        } else if( G.verb > 4 )
            fprintf(stderr,"-- DA[%d] is VALID\n",c);
    }

    /* no check on swapped or compressed */

    if( ! gifti_valid_nvpairs(&gim->ex_atrs, whine) ) errs ++;

    if( G.verb > 2 ) {
        fprintf(stderr,"-- gifti_image: errors = %d", errs);
        if( errs ) fprintf(stderr," (INVALID)\n");
        else       fprintf(stderr," (VALID)\n");
    }

    if( errs ) return 0;
    else       return 1;
}

/*---------------------------------------------------------------------*/
/*! return whether data exists
 *
 *  - darray, each darray[i] and darray[i]->data must be set
 *  
 *  return 1 if true, 0 otherwise
*//*-------------------------------------------------------------------*/
int gifti_image_has_data(const gifti_image * gim)
{
    int c;

    if( !gim || !gim->darray || gim->numDA <= 0 ) return 0;

    for( c = 0; c < gim->numDA; c++ )
        if( !gim->darray[c] ) {
            if(G.verb > 3) fprintf(stderr,"** gim missing data at ind %d\n",c);
            return 0;
        }

    return 1;
}

/*---------------------------------------------------------------------*/
/*! duplicate the given gifti_image struct, optionally including data
 *
 *  Allocate and copy all contents of the gifti_image structure and
 *  sub-structures.  If copy_data is not set, all data pointers within
 *  DataArray elements will be left as NULL.
 *  
 *  return a pointer to the newly allocated structure
*//*-------------------------------------------------------------------*/
gifti_image * gifti_copy_gifti_image(const gifti_image * gold, int copy_data)
{
    gifti_image * gnew;
    int           c, errs = 0;   /* check for errors at each step */

    if( !gold ){ fprintf(stderr,"** copy_gim: input is NULL\n"); return NULL; }

    if(G.verb > 3) fprintf(stderr,"++ copying gifti_image (%s data)...\n",
                           copy_data ? "with" : "without" );

    gnew = (gifti_image *)calloc(1, sizeof(gifti_image));
    if(!gnew){fprintf(stderr,"** copy_gim, failed alloc\n"); return NULL;}

    /* copy one piece at a time */
    gnew->numDA = gold->numDA;
    gnew->version = gifti_strdup(gold->version);

    errs = gifti_copy_nvpairs(&gnew->meta, &gold->meta);
    if(!errs) errs=gifti_copy_LabelTable(&gnew->labeltable, &gold->labeltable);

    /* if needed, create and copy the DataArray list */
    if( !errs && gold->darray && gold->numDA > 0 ) {
        gnew->darray=(giiDataArray**)malloc(gold->numDA*sizeof(giiDataArray*));
        if( !gnew->darray ) {
            fprintf(stderr,"** copy_gim: failed to alloc darray of len %d\n",
                    gold->numDA);
            errs = 1;
        }

        for( c = 0; !errs && c < gold->numDA; c++ ) {
            gnew->darray[c] = gifti_copy_DataArray(gold->darray[c], copy_data);
            if( !gnew->darray[c]) errs++;
        }
    }

    /* and copy the extras */
    if( !errs ) {
        gnew->swapped = gold->swapped;
        gnew->compressed = gold->compressed;
        errs = gifti_copy_nvpairs(&gnew->ex_atrs, &gold->ex_atrs);
    }

    /* on failure, blow everything away */
    if( errs ) { gifti_free_image(gnew); return NULL; }

    return gnew;
}

/*---------------------------------------------------------------------*/
/*! convert all data to NIFTI_TYPE_FLOAT32
 *
 *  for each DataArray
 *      if data exists, convert it (free old, allocate new)
 *      else, leave as NULL
*//*-------------------------------------------------------------------*/
int gifti_convert_to_float(gifti_image * gim)
{
    giiDataArray * da;
    void         * olddata;
    int            oldtype, newtype = NIFTI_TYPE_FLOAT32; /* for future? */
    int            c, nbyper, oldnbyper;

    if( !gim ) return 1;

    if( !gim->darray || gim->numDA <= 0 ) {
        if( G.verb > 1 ) fprintf(stderr,"-- no darray to convert to float\n");
        return 0;
    }

    if( G.verb > 1 ) fprintf(stderr,"++ converting gifti_image to float\n");

    /* first, check for problems in any DA */
    for( c = 0; c < gim->numDA; c++ ) {
        da = gim->darray[c];
        if( !da ) continue;

        /* if the data has an unknown type, panic into error */
        if( da->data && !gifti_valid_datatype(da->datatype, 0) ) {
            fprintf(stderr,"** unknown dtype %d, cannot convert to floats\n",
                    da->datatype);
            return 1;
        }

        /* dimension information must be consistent */
        if( !gifti_valid_dims(da, 1) ) return 1;
    }

    /* will update nbyper in each DA, below */
    gifti_datatype_sizes(newtype, &nbyper, NULL);

    /* all is well, update all DA elements */
    for( c = 0; c < gim->numDA; c++ ) {
        da = gim->darray[c];
        if( !da ) continue;

        oldtype = da->datatype;
        oldnbyper = da->nbyper;

        if( oldtype == newtype ) {
            if(G.verb > 3) fprintf(stderr,"++ convert DA[%d] already type %s\n",
                                   c, gifti_datatype2str(newtype));
            continue; /* no change needed */
        }

        if(G.verb > 3)
            fprintf(stderr,"++ convert DA[%d] from %s to %s\n", c,
                    gifti_datatype2str(oldtype), gifti_datatype2str(newtype));

        if( oldtype == newtype ) continue; /* no change needed */

        /* start trashing DataArray: adjust datatype and nbyper */
        da->datatype = newtype;
        da->nbyper = nbyper;

        /* if there is no data, we are done with this DA */
        if( !da->data ) {
            if( G.verb > 4 ) fprintf(stderr,"-- no data to convert\n");
            continue;
        }

        /* store old data pointer, and allocate new data space */
        olddata = da->data;
        da->data = NULL;    /* so alloc_DA_data doesn't whine */
        if( gifti_alloc_DA_data(gim, &c, 1) ) return 1;

        /* copy the data, and nuke the old stuff */
        if( copy_data_as_float(da->data,newtype,olddata,oldtype,da->nvals) ) {
            /* undo this copy and return */
            free(da->data);
            da->data = olddata; da->datatype = oldtype; da->nbyper = oldnbyper;
            return 1;
        }

        free(olddata);  /* done with this, data has been copied */
    }

    return 0;
}

/* copy old data to float array */
static int copy_data_as_float(void * dest, int dtype, void * src, int stype,
                              long long nvals)
{
    float     * dptr = (float *)dest;
    long long   c;

    if( !dest || !src ) {
        fprintf(stderr,"** copy_data_as_float: missing pointers\n");
        return 1;
    }

    if( dtype != NIFTI_TYPE_FLOAT32 ) {
        fprintf(stderr,"** can't copy to float with dest type %d\n", dtype);
        return 1;
    }

    switch( stype ) {
        default: {
            fprintf(stderr,"** copy2float: can't handle src type %d\n",stype);
            return 1;
        }

        case NIFTI_TYPE_INT8: {
            char * sptr = (char *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_INT16: {
            short * sptr = (short *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_INT32: {
            int * sptr = (int *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_INT64: {
            long long * sptr = (long long *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_UINT8: {
            unsigned char * sptr = (unsigned char *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_UINT16: {
            unsigned short * sptr = (unsigned short *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_UINT32: {
            unsigned int * sptr = (unsigned int *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_UINT64: {
            unsigned long long * sptr = (unsigned long long *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }

        case NIFTI_TYPE_FLOAT32: {
            float * sptr = (float *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = sptr[c];
            break;
        }

        case NIFTI_TYPE_FLOAT64: {
            double * sptr = (double *)src;
            for( c = 0; c < nvals; c++ ) dptr[c] = (float)sptr[c];
            break;
        }
    }

    return 0;
}
