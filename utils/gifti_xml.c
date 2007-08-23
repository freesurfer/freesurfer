/* Author: Rick Reynolds reynoldr@mail.nih.gov */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "expat.h"
#include "gifti.h"
#include "gifti_xml.h"

#ifndef XML_STATUS_ERROR
#define XML_STATUS_ERROR 0
#endif

#define GXML_BSIZE 2048

/* local prototypes */
static int  append_to_cdata     (gxml_data *, const char *, int);
static int  append_to_data      (gxml_data *, const char *, int);
static int  append_to_xform     (gxml_data *, const char *, int);
static int  decode_ascii        (gxml_data*, char*, int, int, void*, int*,int*);
static int  ename2type          (const char *);
static int  epush               (gxml_data *, int, const char *, const char **);
static int  epop                (gxml_data *, int, const char *);
static void init_gxml_data      (gxml_data *, int);
static int  partial_buf_size    (int);

static int  push_gifti          (gxml_data *, const char **);
static int  push_meta           (gxml_data *);
static int  push_md             (gxml_data *);
static int  push_name           (gxml_data *, const char **);
static int  push_value          (gxml_data *, const char **);
static int  push_LT             (gxml_data *, const char **);
static int  push_label          (gxml_data *, const char **);
static int  push_darray         (gxml_data *, const char **);
static int  push_cstm           (gxml_data *);
static int  push_data           (gxml_data *);
static int  push_dspace         (gxml_data *);
static int  push_xspace         (gxml_data *);
static int  push_xform          (gxml_data *, const char **);
static int  push_cdata          (gxml_data *, const char **);
static int  reset_xml_buf       (gxml_data *, char **, int *);

static void show_attrs          (gxml_data *,int,const char **);
static void show_depth          (int, int, FILE *);
static void show_enames         (FILE *);
static int  show_stack          (char *, gxml_data *);

static int  stack_is_valid      (gxml_data *);
static int  whitespace_len      (const char *, int);
static int  update_partial_buffer(char **, int *, int);


static MetaData * find_current_MetaData(gxml_data *, int);

static void XMLCALL cb_start_ele    (void *, const char *, const char **);
static void XMLCALL cb_end_ele      (void *, const char *);
static void XMLCALL cb_char         (void *, const char *, int);
static void XMLCALL cb_instr        (void *, const char *, const char *);
static void XMLCALL cb_comment      (void *, const char *);
static void XMLCALL cb_cdata_start  (void *);
static void XMLCALL cb_cdata_end    (void *);
static void XMLCALL cb_default      (void *, const char *, int);
static void XMLCALL cb_xml_dec      (void *, const char *, const char * , int);
static void XMLCALL cb_start_doctype(void *, const char *, const char *,
                                     const char *, int);
static void XMLCALL cb_end_doctype  (void *);
static void XMLCALL cb_elem_dec     (void *, const char *, XML_Content *);
static XML_Parser init_xml_parser   (void *);

/* writing functions */
static int  gxml_write_gifti(gxml_data *, FILE *);
static int  gxml_write_preamble(gxml_data *, FILE *);

static int  ewrite_cdata_ele        (int, const char *, int, FILE *);
static int  ewrite_coordsys         (gxml_data *, CoordSystem *, FILE *);
static int  ewrite_data             (gxml_data *, DataArray *, FILE *);
static int  ewrite_data_line        (void *, int, int, int, int, FILE *);
static int  ewrite_double_line      (double *, int, int, FILE *);
static int  ewrite_int_attr         (const char *, int, int, int, FILE *);
static int  ewrite_str_attr         (const char*, const char*, int, int, FILE*);
static int  ewrite_darray           (gxml_data *, DataArray *, FILE *);
static int  ewrite_ex_atrs          (gxml_data *, nvpairs *, int, int, FILE *);
static int  ewrite_LT               (gxml_data *, LabelTable *, FILE *);
static int  ewrite_meta             (gxml_data *, MetaData *, FILE *);

/* these should match GXML_ETYPE_* defines */
static char * enames[GXML_MAX_ELEN] = {
    "Invalid", "GIFTI", "MetaData", "MD", "Name", "Value", "LabelTable",
    "Label", "DataArray", "CoordinateSystemTransformMatrix", "Data",
    "DataSpace", "TransformedSpace", "MatrixData", "CDATA"
};

/* ---------------------------------------------------------------------- */
/* GIFTI XML global struct and access functions */
static gxml_data GXD = {
    1,          /* verb, default to 1 (0 means quiet)         */
    1,          /* flag: whether to store data                */
    GXML_BSIZE, /* buf_size, allocated for XML parsing        */
    0,          /* depth, stack depth, shaken, not stirred    */
    0,          /* errors, number of encountered errors       */
    0,          /* skip depth (at positive depth, 0 is clear) */
    {0},        /* stack, ints, max depth GXML_MAX_DEPTH      */

    0,          /* dind, index into decoded data array        */
    0,          /* doff, offset into current data buffer      */
    0,          /* clen, length of current CDATA string       */
    0,          /* xlen, length of current xform buffer       */
    0,          /* dlen, length of current Data buffer        */
    NULL,       /* cdata, CDATA char pointer                  */
    NULL,       /* xdata, xform buffer pointer                */
    NULL,       /* ddata, Data buffer pointer                 */
    NULL        /* gifti_image *, for results                 */
};


/* note: the buffer needs to be large enough to contain any contiguous
         piece of (CDATA?) text, o.w. it will require parsing in pieces */
gifti_image * gxml_read_image( const char * fname, int read_data )
{
    gxml_data  * xd = &GXD;     /* point to global struct */
    XML_Parser   parser;
    FILE       * fp;
    char       * buf = NULL;
    int          bsize;    /* be sure it doesn't change at some point */
    int          done = 0, blen;
    int          pcount = 1;
 
    init_gxml_data(xd, 0); /* reset non-user variables */
    xd->dstore = read_data;  /* store for global access */

    if( !fname ) {
        fprintf(stderr,"** gxml_read_image: missing filename\n");
        return NULL;
    }

    fp = fopen(fname, "r");
    if( !fp ) {
        fprintf(stderr,"** failed to open GIFTI xml file '%s'\n", fname);
        return NULL;
    }

    /* create a new buffer */
    bsize = 0;
    if( reset_xml_buf(xd, &buf, &bsize) ) { fclose(fp); return NULL; }

    if(xd->verb > 2) fprintf(stderr,"-d reading gifti image '%s'\n", fname);
    if(xd->verb > 2) fprintf(stderr,"-d using %d byte XML buffer\n",bsize);
    if(xd->verb > 3) show_enames(stderr);

    /* allocate return structure */
    xd->gim = (gifti_image *)calloc(1,sizeof(gifti_image));
    if( !xd->gim ) {
        fprintf(stderr,"** failed to alloc initial gifti_image\n");
        free(buf);
        return NULL;
    }

    /* create parser, init handlers */
    parser = init_xml_parser((void *)xd);

    while( !done )
    {
        if( reset_xml_buf(xd, &buf, &bsize) )
        { gifti_free_image(xd->gim); xd->gim = NULL; break; }  /* fail out */

        blen = (int)fread(buf, 1, bsize, fp);
        done = blen < sizeof(buf);

        if(xd->verb > 4) fprintf(stderr,"-d XML_Parse # %d\n", pcount);
        pcount++;
        if( XML_Parse(parser, buf, blen, done) == XML_STATUS_ERROR) {
            fprintf(stderr,"** %s at line %u\n",
                    XML_ErrorString(XML_GetErrorCode(parser)),
                    (unsigned int)XML_GetCurrentLineNumber(parser));
            gifti_free_image(xd->gim);
            xd->gim = NULL;
            break;
        }
    }

    if(xd->verb > 1) {
        if(xd->gim)
            fprintf(stderr,"-d gifti image '%s', success "
                           "(%d DA elements = %d MB)\n",
                    fname, xd->gim->numDA, gifti_gim_DA_size(xd->gim,1));
        else fprintf(stderr,"** gifti image '%s', failure\n", fname);
    }

    fclose(fp);
    if( buf       ) free(buf);          /* parser buffer */
    if( xd->xdata ) free(xd->xdata);    /* xform matrix buffer */
    if( xd->ddata ) free(xd->ddata);    /* Data buffer */
    XML_ParserFree(parser);

    return xd->gim;
}


/* return 0 on success */
int gxml_write_image(gifti_image * gim, const char * fname, int write_data)
{
    gxml_data * xd = &GXD;     /* point to global struct */
    FILE      * fp;

    if( !gim ) {
        fprintf(stderr,"** GXML write: no gifti_image\n");
        return 1;
    } else if ( !fname ) {
        fprintf(stderr,"** GXML write: no filename\n");
        return 1;
    }

    if(GXD.verb > 1) {
        fprintf(stderr,"+d writing gifti image (%s data) to '%s'",
                write_data?"with":"no", fname);
        if( write_data )
            fprintf(stderr," (%d DA elements = %d MB)",
                    gim->numDA, gifti_gim_DA_size(xd->gim,1));
        fputc('\n', stderr);
    }

    init_gxml_data(xd, 0);    /* reset non-user variables */
    xd->dstore = write_data;  /* store for global access */
    xd->skip = 3;  /* rcr - how to decide this 'spaces per indent' */
    xd->ddata = strdup(fname);
    xd->dlen = strlen(fname);
    xd->gim = gim;

    fp = fopen(fname, "w");
    if( !fp ) {
        fprintf(stderr,"** failed to open '%s' for gifti write\n", fname);
        return 1;
    }

    (void)gxml_write_gifti(xd, fp);

    if( xd->ddata ) free(xd->ddata);
    if( xd->xdata ) free(xd->xdata);

    fclose(fp);

    return 0;
}


/* maybe these can be enhanced tomorrow, tomorrow, tomorrow... */
int gxml_set_verb( int val ){ GXD.verb = val; return 0; }
int gxml_get_verb( void    ){ return GXD.verb; }

int gxml_set_dstore( int val ){ GXD.dstore = val; return 0; }
int gxml_get_dstore( void    ){ return GXD.dstore; }

/* buf_size is applied only at main reading time, for now */
int gxml_set_buf_size( int val ){ GXD.buf_size = val; return 0; }
int gxml_get_buf_size( void    ){ return GXD.buf_size; }


static void init_gxml_data( gxml_data * dp, int doall )
{
    if( doall ) {       /* user modifiable */
        dp->verb = 1;
        dp->dstore = 1;
        dp->buf_size = GXML_BSIZE;
    }

    dp->errors = 0;
    dp->skip = 0;
    dp->depth = 0;
    memset(dp->stack, 0, sizeof(dp->stack));

    dp->dind = 0;
    dp->doff = 0;
    dp->clen = 0;
    dp->xlen = 0;
    dp->dlen = 0;
    dp->cdata = NULL;
    dp->xdata = NULL;
    dp->ddata = NULL;
    dp->gim = NULL;
}

/* ---------------------------------------------------------------------- */

static void show_depth( int depth, int show, FILE * fp )
{
    if( show ) fprintf(fp, "%*s %02d ", 3*depth, " ", depth);
    else       fprintf(fp, "%*s    ", 3*depth, " ");
}

static void show_enames( FILE * fp )
{
    int c;
    fprintf(fp, "-------------------------------\n"
                "+d ename list :\n");
    for( c = 0; c <= GXML_ETYPE_LAST; c++ )
        fprintf(fp,"    %02d : %s\n", c, enames[c]);
    fprintf(fp, "-------------------------------\n");
}

static int ename2type( const char * name )
{
    int etype;
    for( etype = GXML_ETYPE_LAST; etype > GXML_ETYPE_INVALID; etype-- )
        if( !strcmp(name, enames[etype]) )
            break;
    return etype;
}

/* name should be null terminated */
static int epush( gxml_data * xd, int etype, const char * ename,
                                                  const char ** attr )
{
    if( xd->depth < 0 || xd->depth > GXML_MAX_DEPTH ) {
        fprintf(stderr,"** push: stack depth %d out of [0,%d] range\n",
                xd->depth, GXML_MAX_DEPTH);
        xd->errors++;
        return 1;
    }

    if( xd->verb > 2 ) {       /* maybe we want to print something */
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr,"+d push %02d: '%s'\n", etype, enames[etype]);
    }

    xd->stack[xd->depth] = etype;
    xd->depth++;

    /* if we are in a skip block, do nothing but monitor stack */
    if( xd->skip ) {
        if( xd->verb > 1 )
            fprintf(stderr,"-d skip=%d, depth=%d, skipping element '%s'\n",
                    xd->skip, xd->depth, ename);
        return 0;
    }

    /* determine whether we should enter a skip block */
    if( etype == GXML_ETYPE_INVALID ) {
        if(xd->verb > 0)
            fprintf(stderr,"** pushed invalid element, '%s', skip depth %d\n",
                    ename, xd->depth);
        xd->skip = xd->depth;
        return 1;
    }

    if ( xd->verb > 4 ) show_stack("+d ", xd);
    if( !stack_is_valid(xd) ) return 1;

    /* call appropriate XML processing function */
    switch( etype ) {
        case GXML_ETYPE_GIFTI      : return push_gifti (xd, attr);
        case GXML_ETYPE_META       : return push_meta  (xd);
        case GXML_ETYPE_MD         : return push_md    (xd);
        case GXML_ETYPE_NAME       : return push_name  (xd, attr);
        case GXML_ETYPE_VALUE      : return push_value (xd, attr);
        case GXML_ETYPE_LABELTABLE : return push_LT    (xd, attr);
        case GXML_ETYPE_LABEL      : return push_label (xd, attr);
        case GXML_ETYPE_DATAARRAY  : return push_darray(xd, attr);
        case GXML_ETYPE_CSTM       : return push_cstm  (xd);
        case GXML_ETYPE_DATA       : return push_data  (xd);
        case GXML_ETYPE_DATASPACE  : return push_dspace(xd);
        case GXML_ETYPE_XFORMSPACE : return push_xspace(xd);
        case GXML_ETYPE_MATRIXDATA : return push_xform (xd, attr);
        case GXML_ETYPE_CDATA      : return push_cdata (xd, attr);
        default: /* drop through */
            break;
    }

    fprintf(stderr,"** epush, unknow type '%s'\n",enames[etype]);
    return 1;
}

/* initialize the gifti_element and set attributes */
static int push_gifti(gxml_data * xd, const char ** attr )
{
    gifti_image *     gim;
    int               c;
    if( !xd ) return 1;
    if( !attr ) return 0;

    /* be explicit with pointers (struct should be clear) */
    gim = xd->gim;
    gim->version = NULL;
    clear_nvpairs(&gim->meta);
    clear_LabelTable(&gim->labeltable);
    gim->darray = NULL;
    clear_nvpairs(&gim->ex_atrs);

    for(c = 0; attr[c]; c+= 2 )
        if( gifti_str2attr_gifti(gim, attr[c], attr[c+1]) )
            if( gifti_add_to_nvpairs(&gim->ex_atrs,attr[c],attr[c+1]) )
                return 1;

    if( xd->verb > 2 ) fprintf(stderr,"+d set %d GIFTI attr(s)\n",c/2);
    if( xd->verb > 3 ) gifti_disp_gifti_image("push:", gim, 0);

    return 0;
}

/* simply verify that we have not been here before */
static int push_meta(gxml_data * xd)
{
    MetaData * md = find_current_MetaData(xd, 0);  /* MD is 1 below MetaData */

    if( md->length != 0 || md->name || md->value ) {
        fprintf(stderr,"** push meta: already initialized??\n");
        return 1;
    }

    return 0;
}

/* find the parent struct, and return its meta field */
static MetaData * find_current_MetaData(gxml_data * xd, int cdepth)
{
    DataArray * da;
    MetaData  * md;
    int         da_ind, parent;
    if( !xd || cdepth < 0 || xd->depth < (2+cdepth) ) {
        fprintf(stderr,"FMeta: bad params (%p,%d)\n",xd,cdepth);
        return NULL;
    }

    /* find the appropriate parent struct */
    parent = xd->stack[xd->depth-2-cdepth];
    if( parent == GXML_ETYPE_GIFTI )
        md = &xd->gim->meta;
    else if( parent == GXML_ETYPE_DATAARRAY ) {
        if( !xd->gim->darray ) {
            fprintf(stderr,"** FMeta: gim->darry not initialized\n");
            return NULL;
        }
        da_ind = xd->gim->numDA-1;
        da = xd->gim->darray[da_ind];
        if( !da ) {
            fprintf(stderr,"** FMeta: gim->darry[%d] not initialized\n",da_ind);
            return NULL;
        }
        md = &da->meta;
    } else {
        fprintf(stderr,"** FMeta: child of invalid parent '%s'\n",
                enames[parent]);
        return NULL;
    }

    return md;
}

/* we will add a pair, so update length and allocate pointers */
static int push_md(gxml_data * xd)
{
    MetaData * md = find_current_MetaData(xd, 1);  /* MD is 1 below MetaData */
    if( !md ) return 1;  /* error were printed */

    md->length++;
    md->name = (char **)realloc(md->name, md->length * sizeof(char *));
    md->value = (char **)realloc(md->value, md->length * sizeof(char *));

    if( !md->name || !md->value ) {
        fprintf(stderr,"** failed to realloc %d MD pointers\n",md->length);
        md->length = 0;
        return 1;
    }

    /* and clear the new pointers */
    md->name[md->length-1] = NULL;
    md->value[md->length-1] = NULL;

    return 0;
}

/* set cdata to the current meta->name address, and clear it */
static int push_name(gxml_data * xd, const char ** attr)
{
    MetaData * md = find_current_MetaData(xd, 2);  /* name is 2 below Meta */
    if( !md ) return 1;

    xd->cdata = &md->name[md->length-1];  /* use cdata to fill */
    *xd->cdata = NULL;                    /* init to empty */
    xd->clen = 0;

    return 0;
}

/* set cdata to the current meta->value address, and clear it */
static int push_value(gxml_data * xd, const char ** attr)
{
    MetaData * md = find_current_MetaData(xd, 2);  /* name is 2 below Meta */
    if( !md ) return 1;

    xd->cdata = &md->value[md->length-1];  /* use cdata to fill */
    *xd->cdata = NULL;                     /* init to empty */
    xd->clen = 0;

    return 0;
}

/* initialize the gifti_element and set attributes */
static int push_LT(gxml_data * xd, const char ** attr)
{
    LabelTable * lt = &xd->gim->labeltable;
    if( lt->length || lt->index || lt->label ) {
        fprintf(stderr,"** multiple LabelTables?\n");
    }

    return 0;
}

/* initialize the gifti_element and set attributes */
static int push_label(gxml_data * xd, const char ** attr)
{
    LabelTable * lt = &xd->gim->labeltable;

    lt->length++;
    lt->index = (int *)realloc(lt->index, lt->length * sizeof(int));
    lt->label = (char **)realloc(lt->label, lt->length * sizeof(char *));

    /* set index from the attributes */
    if( !attr[0] || strcmp(attr[0],"Index"))
        lt->index[lt->length-1] = 0;
    else
        lt->index[lt->length-1] = atoi(attr[1]);

    xd->cdata = lt->label + (lt->length-1); /* addr of newest (char *) */
    *xd->cdata = NULL;                      /* init to empty */
    xd->clen = 0;

    return 0;
}

/* initialize the gifti_element and set attributes */
static int push_darray(gxml_data * xd, const char ** attr)
{
    DataArray * da;
    int         buf_size;

    if( gifti_add_empty_darray(xd->gim) ) return 1;

    da = xd->gim->darray[xd->gim->numDA-1];  /* get new pointer */

    /* fill the struct from the attributes */
    if( gifti_init_darray_from_attrs(da, attr) ) return 1;

    /* make a request to potentially update the XML buffer size */
    if( da->nvals>0 && da->nbyper>0 ) {
        buf_size = partial_buf_size(da->nvals*da->nbyper);
        if( buf_size != xd->buf_size ) {
            if( xd->verb > 2 )
                fprintf(stderr,"+d update XML buf size, %d to %d (for %u)\n",
                    xd->buf_size, buf_size, (unsigned)da->nvals*da->nbyper);
            xd->buf_size = buf_size;
        }
    }

    if( xd->verb > 4 ) gifti_disp_DataArray("push:", da, 0);

    return 0;
}

/* verify the elements are clear */
static int push_cstm(gxml_data * xd)
{
    DataArray * da = xd->gim->darray[xd->gim->numDA-1];  /* get new pointer */

    da->coordsys = (CoordSystem *)malloc(sizeof(CoordSystem));
    clear_CoordSystem(da->coordsys);

    return 0;
}

/* verify the processing buffer space, alloc data space */
static int push_data(gxml_data * xd)
{
    DataArray * da = xd->gim->darray[xd->gim->numDA-1];  /* get cur pointer */

    xd->dind = 0;       /* init for filling */
    xd->doff = 0;

    /* partial buffer for processing space */
    if( update_partial_buffer(&xd->ddata, &xd->dlen, da->nbyper*da->nvals) )
        return 1;

    /* allocate space for data */
    if( da->nvals <= 0 || da->nbyper <= 0 ) {
        fprintf(stderr,"** PD: bad vals,bytes = %u, %d\n",
                (unsigned)da->nvals,da->nbyper);
        return 1;
    }

    da->data = calloc(da->nvals, da->nbyper);
    if( ! da->data ) {
        fprintf(stderr,"** PD: failed to alloc %u bytes for darray[%d]\n",
                (unsigned)da->nvals*da->nbyper, xd->gim->numDA-1);
        return 1;
    } else if ( xd->verb > 3 )
        fprintf(stderr,"+d PD: alloc %u bytes for darray[%d]\n",
                (unsigned)da->nvals*da->nbyper, xd->gim->numDA-1);

    return 0;
}

/* point cdata to the correct location and init */
static int push_dspace(gxml_data * xd)
{
    if( !xd->gim->darray[xd->gim->numDA-1]->coordsys ) {
        fprintf(stderr,"** found dataspace without coordsys, skipping...\n");
        xd->skip = xd->depth;
        return 1;
    }

    xd->cdata = &xd->gim->darray[xd->gim->numDA-1]->coordsys->dataspace;
    *xd->cdata = NULL;                      /* init to empty */
    xd->clen = 0;
    return 0;
}

/* point cdata to the correct location and init */
static int push_xspace(gxml_data * xd)
{
    if( !xd->gim->darray[xd->gim->numDA-1]->coordsys ) {
        fprintf(stderr,"** found xformspace without coordsys, skipping...\n");
        xd->skip = xd->depth;
        return 1;
    }

    xd->cdata = &xd->gim->darray[xd->gim->numDA-1]->coordsys->xformspace;
    *xd->cdata = NULL;                      /* init to empty */
    xd->clen = 0;
    return 0;
}

/* verify the processing buffer space */
static int push_xform(gxml_data * xd, const char ** attr)
{
    if( !xd->gim->darray[xd->gim->numDA-1]->coordsys ) {
        fprintf(stderr,"** found xform without coordsys, skipping...\n");
        xd->skip = xd->depth;
        return 1;
    }

    /* just make sure we have a text buffer to work with */
    if( !xd->xdata || xd->xlen <= 0 ) {
        xd->xlen = 2048;
        xd->xdata = (char *)malloc(xd->xlen * sizeof(char));
        if( !xd->xdata ) {
            fprintf(stderr,"** cannot alloc %d bytes for xform\n",xd->xlen);
            return 1;
        }
    }

    xd->dind = 0;       /* init for filling */
    xd->doff = 0;

    return 0;
}

/* if we are in a char append state, append */
static int push_cdata(gxml_data * xd, const char ** attr)
{
    return 0;
}

static int epop( gxml_data * xd, int etype, const char * ename )
{
    if( xd->verb > 3 )
    {
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr,"+d pop %02d : '%s'\n", etype, enames[etype]);
    }

    xd->cdata = NULL;                   /* clear fields for future use */
    xd->clen = 0;

    if( xd->skip == xd->depth ) {       /* just completed skip element */
        if( xd->verb > 1 )
            fprintf(stderr,"-d popping skip element '%s' at depth %d\n",
                    ename, xd->depth);
        xd->skip = 0;  /* clear skip level */
    } else {    /* may peform pop action for this element */
        switch( etype ) {
/*        case GXML_ETYPE_GIFTI      : return push_gifti (xd, attr);
            case GXML_ETYPE_META       : return push_meta  (xd);
            case GXML_ETYPE_MD         : return push_md    (xd);
            case GXML_ETYPE_NAME       : return push_name  (xd, attr);
            case GXML_ETYPE_VALUE      : return push_value (xd, attr);
            case GXML_ETYPE_LABELTABLE : return push_LT    (xd, attr);
            case GXML_ETYPE_LABEL      : return push_label (xd, attr);
            case GXML_ETYPE_DATAARRAY  : return push_darray(xd, attr);
            case GXML_ETYPE_CSTM       : return push_cstm  (xd, attr);
            case GXML_ETYPE_DATA       : return push_data  (xd, attr);
            case GXML_ETYPE_DATASPACE  : return push_dspace(xd, attr);
            case GXML_ETYPE_XFORMSPACE : return push_xspace(xd, attr);
            case GXML_ETYPE_CDATA      : return push_cdata (xd, attr); */
            default: /* drop through */
                break;
            case GXML_ETYPE_MATRIXDATA :
            case GXML_ETYPE_DATA       :        /* clear dind/doff? */
                break;

            case GXML_ETYPE_GIFTI      :
                if(xd->verb > 4) gifti_disp_gifti_image("pop:",xd->gim,1);
                break;
        }
    }

    xd->depth--;

    if( xd->depth < 0 || xd->depth > GXML_MAX_DEPTH ) {
        fprintf(stderr,"** pop: stack depth %d out of [0,%d] range\n",
                xd->depth, GXML_MAX_DEPTH);
        xd->errors++;
        return -1;
    }

    return 0;
}

/* return the number of bytes of leading whitespace, up to a max of len */
static int whitespace_len(const char * str, int len)
{
    int c;
    if( !str || !*str || len < 1 ) return 0;
    for( c = 0; c < len; c++ )
        if( !isspace(str[c]) ) return c;

    return len;
}

static void show_attrs(gxml_data * xd, int etype, const char ** attr)
{
    int count;
    show_depth(xd->depth, 1, stderr);
    fprintf(stderr, ": element %s\n", enames[etype]);
    for( count = 0; attr[count]; count += 2 ){
        show_depth(xd->depth, 0, stderr);
        fprintf(stderr,"      attr: %s='%s'\n", attr[count], attr[count+1]);
    }
}


static void XMLCALL cb_start_ele(void *udata, const char *ename,
                                              const char **attr)
{
    gxml_data * xd = (gxml_data *)udata;
    int         etype;

    etype = ename2type(ename);
    if( xd->verb > 3 ) show_attrs(xd, etype, attr);

    /* process attributes and push() */

    (void)epush(xd, etype, ename, attr);
}

/* if ending Data, clear prev_end_check */
static void XMLCALL cb_end_ele(void *udata, const char * ename)
{
    gxml_data * xd = (gxml_data *)udata;

    epop(xd, ename2type(ename), ename);
}

/* May divide Data, but apparently not attributes, perhaps because
   the Data section is longer than the buffer is wide.

   if in Data:
        if prev_end_check
            if( ! char_is_whitespace(first) && concat_is_number() )
                concatencate as number to adjust previous number
                verify rest is whitespace
                return  (we don't expect to start a new number)
            else
                apply previous number
        if( !char_is_whitespace(last) )
            store trailing non-space in concat_buf
        else
            prev_end_check = 0
        apply number (though it may change later)
*/
static void XMLCALL cb_char(void *udata, const char * cdata, int length)
{
    gxml_data * xd = (gxml_data *)udata;
    const char * str = cdata;
    int          len = length, wlen = 0, parent;

    if( xd->skip > 0 ) {
        if(xd->verb > 2) fprintf(stderr,"-d skipping char [%d]\n",len);
        return;
    }

    /* act based on the parent type */
    parent = xd->stack[xd->depth-1];
    if( parent == GXML_ETYPE_CDATA ) parent = xd->stack[xd->depth-2];

    if( parent != GXML_ETYPE_DATA ) wlen = whitespace_len(str,length);

    switch( parent ) {
        case GXML_ETYPE_DATA       :
            (void)append_to_data(xd, cdata, length);
            break;
        case GXML_ETYPE_MATRIXDATA :
            (void)append_to_xform(xd, cdata, length);
            break;

        case GXML_ETYPE_GIFTI      :
        case GXML_ETYPE_META       :
        case GXML_ETYPE_MD         :
        case GXML_ETYPE_LABELTABLE :
        case GXML_ETYPE_DATAARRAY  :
        case GXML_ETYPE_CSTM       :
            if( wlen != length && xd->verb ) {
                fprintf(stderr,"** invalid chars under %s: '%.*s'\n",
                        enames[parent], length, cdata);
            }
            break;

        case GXML_ETYPE_NAME       :
        case GXML_ETYPE_VALUE      :
        case GXML_ETYPE_LABEL      :
        case GXML_ETYPE_DATASPACE  :
        case GXML_ETYPE_XFORMSPACE :
            if( xd->verb > 4 )
                fprintf(stderr,"+d append cdata, parent %s\n",enames[parent]);
            (void)append_to_cdata(xd, cdata, length);
            break;

        case GXML_ETYPE_CDATA      :
            fprintf(stderr,"** CDATA is the parent of CDATA???\n");
            return;

        default: /* drop through */
            fprintf(stderr,"** unknown parent of char: %d\n", parent);
            return;
    }

    if( wlen == length ) {      /* if it is all whitespace */
        if( xd->verb < 5 ) return;
        str = "whitespace";     /* just note the whitespace */
        len = strlen(str);
    }

    if( xd->verb > 4 ) {
        show_depth(xd->depth, 1, stderr);
        if( parent == GXML_ETYPE_DATA && len > 40 ) len = 40;
        fprintf(stderr, "char[%d]: %.*s\n", length, len, str);
    }
}

static int append_to_cdata(gxml_data * xd, const char * cdata, int len)
{
    int offset;
    if( !xd || !cdata || len <= 0 ) {
        fprintf(stderr,"** A2CD, bad params (%p,%p,%d)\n",xd, cdata, len);
        return 1;
    }
    if( !*xd->cdata ) {
        offset = 0;
        xd->clen = len + 1;  /* first time, alloc for null */
    }
    else {
        offset = xd->clen - 1;
        xd->clen += len;
    }

    if( xd->verb > 4 )
        fprintf(stderr,"+d a2cdata, len %d, clen %d, data '%.*s'\n",
                len, xd->clen, len, cdata);

    *xd->cdata = (char *)realloc(*xd->cdata, xd->clen*sizeof(char));
    if(!*xd->cdata) {
        fprintf(stderr,"** A2CD, failed to realloc %d bytes\n",xd->clen);
        return 1;
    }

    memcpy(*xd->cdata + offset, cdata, len);    /* append the new data */
    (*xd->cdata)[xd->clen-1] = '\0';            /* and null terminate */

    return 0;
}


/* this must go to the data of the latest darray struct */
static int append_to_data(gxml_data * xd, const char * cdata, int len)
{
    static int  mod_prev = 0;
    DataArray * da = xd->gim->darray[xd->gim->numDA-1];
    char      * dptr;
    char      * cptr;
    int         rem_vals, rem_len = len, copy_len, unused;
    int         type = da->datatype;

    if( !da || !xd->dlen || !xd->ddata || xd->dind < 0 ) {
        fprintf(stderr,"** A2D: bad setup (%p,%d,%p,%d)\n",
                da, xd->dlen, xd->ddata, xd->dind);
        return 1;
    } else if( da->encoding != GIFTI_ENCODING_ASCII ) {
        fprintf(stderr,"** not ASCII data encoding (%s)\n",
                gifti_encoding_list[da->encoding]);
        return 1;
    } else if( !da->data ) {
        fprintf(stderr,"** A2D: no data allocated\n");
        return 1;
    } else if( xd->verb > 4 )
        fprintf(stderr,"+d appending %d bytes to data\n",len);

    /* if there is only whitespace, blow outta here */
    if( whitespace_len(cdata, len) == len ) { xd->doff = 0; return 0; }

    /* Copy cdata to local buffer in pieces, for null termination and for
       storage of trailing characters that may need to be processed again
       (after more characters are read by the parser).                    */
    while( rem_len > 0 ) {
        /*--- prepare intermediate buffer ---*/

        /* point to the current location */
        cptr = (char *)cdata + len - rem_len;

        /* if we're looking at whitespace, any unused data is garbage */
        if( isspace(*cptr)) xd->doff = 0;

        /* decide how many bytes to copy (avail space w/max of rem_len) */
        copy_len = xd->dlen - xd->doff - 1;
        if( copy_len > rem_len ) {
            unused = copy_len - rem_len;  /* unused at end of buffer */
            copy_len = rem_len;
        } else unused = 0;

        /* copy it to our buffer and null terminate */
        memcpy(xd->ddata+xd->doff, cptr, copy_len);
        xd->ddata[xd->doff+copy_len] = '\0';

        /*--- process the ascii data ---*/

        /* note how many values remain to be computed */
        rem_vals = da->nvals - xd->dind;
        if(xd->verb > 5)
            fprintf(stderr,"-d %d vals left at offset %d, nbyper %d\n",
                    rem_vals, xd->dind, da->nbyper);

        if( xd->dind == 0 ) mod_prev = 0;       /* nothing to modify at first */
        dptr = (char *)da->data + (xd->dind)*da->nbyper;
        xd->doff = decode_ascii(xd,
                        xd->ddata,              /* data source */
                        xd->doff+copy_len,      /* data length */
                        type,                   /* data type */
                        dptr,                   /* starting destination */
                        &rem_vals,              /* nvals to read */
                        &mod_prev               /* can we mod previous val */
                        );

        /*--- check results --- */
        if( xd->doff < 0 ) { xd->doff = 0; return 1; } /* error */
        if( xd->doff >= xd->dlen - 1 ) {
            if(xd->verb) fprintf(stderr,"** A2D: failed to process buffer\n");
            fprintf(stderr,"** rem = %d\n", xd->doff);
            xd->doff = 0;        /* blow away the buffer and continue */
        }

        /*--- adjust intermediate buffer ---*/

        /* move any unused bytes to the beginning (last doff, before unused) */
        if( xd->doff > 0 ) {
            if( xd->verb > 5 )
                fprintf(stderr,"+d A2D: move %d bytes from %d (blen %d)\n",
                    xd->doff, xd->dlen - unused - xd->doff, xd->dlen);
            /* (subtract unused+1, since 1 bytes is saved for null */
            memmove(xd->ddata, xd->ddata+xd->dlen -(unused+1) - xd->doff,
                    xd->doff);
            if( xd->verb > 6 )
                fprintf(stderr,"   bytes are '%.*s'\n",xd->doff,
                        (char *)xd->ddata);
        }

        /* adjust rem_len for next time */
        rem_len -= copy_len;
        xd->dind = da->nvals - rem_vals;  /* note remaining values */
    }

    return 0;
}

/* this must go to the xform of the latest darray struct */
/* (process as 1-D array) */
static int append_to_xform(gxml_data * xd, const char * cdata, int len)
{
    static int  mod_prev = 0;
    DataArray * da = xd->gim->darray[xd->gim->numDA-1];
    double    * dptr;
    char      * cptr;
    int         rem_vals, rem_len = len, copy_len, unused;
    int         type = gifti_str2datatype("NIFTI_TYPE_FLOAT64"); /* double */

    if( !da || !xd->xlen || !xd->xdata || xd->dind < 0 ) {
        fprintf(stderr,"** A2X: bad setup (%p,%d,%p,%d)\n",
                da, xd->xlen, xd->xdata, xd->dind);
        return 1;
    } else if( xd->verb > 4 )
        fprintf(stderr,"+d appending %d bytes to xform\n",len);

    /* if there is only whitespace, blow outta here */
    if( whitespace_len(cdata, len) == len ) { xd->doff = 0; return 0; }

    /* Copy cdata to local buffer in pieces, for null termination and for
       storage of trailing characters that may need to be processed again
       (after more characters are read by the parser).                    */
    while( rem_len > 0 ) {
        /*--- prepare intermediate buffer ---*/

        /* point to the current location */
        cptr = (char *)cdata + len - rem_len;

        /* if we're looking at whitespace, any unused data is garbage */
        if( isspace(*cptr)) xd->doff = 0;

        /* decide how many bytes to copy (avail space w/max of rem_len) */
        copy_len = xd->xlen - xd->doff - 1;
        if( copy_len > rem_len ) {
            unused = copy_len - rem_len;  /* unused at end of buffer */
            copy_len = rem_len;
        } else unused = 0;

        /* copy it to our buffer and null terminate */
        memcpy(xd->xdata+xd->doff, cptr, copy_len);
        xd->xdata[xd->doff+copy_len] = '\0';

        /* note how many values remain to be computed */
        rem_vals = 16 - xd->dind;

        /*--- process the ascii data ---*/
        if( xd->dind == 0 ) mod_prev = 0;       /* nothing to modify at first */
        dptr = (double *)da->coordsys->xform + (xd->dind);  /* as array */
        xd->doff = decode_ascii(xd,
                        xd->xdata,              /* data source */
                        xd->doff+copy_len,      /* data length */
                        type,                   /* data type */
                        dptr,                   /* starting destination */
                        &rem_vals,              /* nvals to read */
                        &mod_prev               /* can we mod previous val */
                        );

        /*--- check results --- */
        if( xd->doff < 0 ) { xd->doff = 0; return 1; } /* error */
        if( xd->doff >= xd->xlen - 1 ) {
            if(xd->verb) fprintf(stderr,"** A2X: failed to process buffer\n");
            fprintf(stderr,"** rem = %d\n", xd->doff);
            xd->doff = 0;        /* blow away the buffer and continue */
        }

        /*--- adjust intermediate buffer ---*/

        /* move any unused bytes to the beginning (last doff, before unused) */
        if( xd->doff > 0 ) {
            if( xd->verb > 5 )
                fprintf(stderr,"+d A2X: move %d bytes from %d (blen %d)\n",
                        xd->doff, xd->dlen - unused - xd->doff, xd->dlen);
                /* (subtract unused+1, since 1 bytes is saved for null */
                memmove(xd->xdata, xd->xdata+xd->xlen -(unused+1) -xd->doff,
                        xd->doff);
            if( xd->verb > 6 )
                fprintf(stderr,"   bytes are '%.*s'\n",xd->doff,
                        (char *)xd->ddata);
        }

        /* adjust rem_len for next time */
        rem_len -= copy_len;
        xd->dind = 16 - rem_vals;  /* note remaining values */
    }

    return 0;
}

/* given: source pointer, remaining length, nvals desired, dest loc and type
          (cdata is null-terminated)
   output: nvals processed, new dest location 
   return: nbytes that may still need to processed (< 0 on error)

   read failure happens only when no characters are processed
*/
static int decode_ascii(gxml_data * xd, char * cdata, int cdlen, int type,
                        void * dptr, int * nvals, int * mod_prev)
{
    char * p1, *p2;     /* for strtoX */
    char * prev;        /* for remain */
    double dval;        /* for strtod */
    long   lval;        /* for strtol */
    int    remain = 0;  /* use bytes remaining */
    int    vals = 0;

    /* if reprocessing, maybe let the user know */
    if( xd->verb > 4)
        fprintf(stderr,"-d DA: type %s, len %d, nvals %d\n",
                gifti_datatype2str(type),cdlen,*nvals);

    if( xd->doff > 0 && *mod_prev ) {
        if( xd->verb > 4)
            fprintf(stderr,"+d DA: re-proc '%.*s' from '%.*s'...\n",
                    xd->doff, cdata, xd->doff+15, cdata);
        vals--;  /* back up */
    }

    switch( type ) {
        default : 
            fprintf(stderr,"** decode_ascii cannot decode type %d\n",type);
            return -1;
        case 2: {       /* NIFTI_TYPE_UINT8 */
            unsigned char * ptr = (unsigned char *)dptr;
            p1 = cdata;
            prev = p1;
            while( vals < *nvals && p1 ) {
                lval = strtol(p1, &p2, 10);   /* try to read next value */
                if( p1 == p2 ) break;   /* nothing read, terminate loop */
                prev = p1;              /* store old success ptr */
                p1 = p2;                /* move to next posn */
                ptr[vals] = lval;       /* assign new value  */
                if(xd->verb>6)fprintf(stderr,"  v %d (%ld)",ptr[vals],lval);
                vals++;                 /* count new value   */
            }
            if(xd->verb > 6) fputc('\n', stderr);
            break;
        }
        case 4: {       /* NIFTI_TYPE_INT16 */
            short * ptr = (short *)dptr;
            p1 = cdata;
            prev = p1;
            while( vals < *nvals && p1 ) {
                lval = strtol(p1, &p2, 10);   /* try to read next value */
                if( p1 == p2 ) break;   /* nothing read, terminate loop */
                prev = p1;              /* store old success ptr */
                p1 = p2;                /* move to next posn */
                ptr[vals] = lval;       /* assign new value  */
                if(xd->verb>6)fprintf(stderr,"  v %d (%ld)",ptr[vals],lval);
                vals++;                 /* count new value   */
            }
            if(xd->verb > 6) fputc('\n', stderr);
            break;
        }
        case 8: {       /* NIFTI_TYPE_INT32 */
            int * ptr = (int *)dptr;
            p1 = cdata;
            prev = p1;
            while( vals < *nvals && p1 ) {
                lval = strtol(p1, &p2, 10);   /* try to read next value */
                if( p1 == p2 ) break;   /* nothing read, terminate loop */
                prev = p1;              /* store old success ptr */
                p1 = p2;                /* move to next posn */
                ptr[vals] = lval;       /* assign new value  */
                if(xd->verb>6)fprintf(stderr,"  v %d (%ld)",ptr[vals],lval);
                vals++;                 /* count new value   */
            }
            if(xd->verb > 6) fputc('\n', stderr);
            break;
        }
        case 16: {      /* NIFTI_TYPE_FLOAT32 */
            float * ptr = (float *)dptr;
            p1 = cdata;
            prev = p1;
            while( vals < *nvals && p1 ) {
                dval = strtod(p1, &p2); /* try to read next value */
                if( p1 == p2 ) break;   /* nothing read, terminate loop */
                prev = p1;              /* store old success ptr */
                p1 = p2;                /* move to next posn */
                ptr[vals] = dval;       /* assign new value  */
                if(xd->verb>6)fprintf(stderr,"  v %f (%f)",ptr[vals],dval);
                vals++;                 /* count new value   */
            }
            if(xd->verb > 6) fputc('\n', stderr);
            break;
        }
        case 64: {      /* NIFTI_TYPE_FLOAT64 */
            double * ptr = (double *)dptr;
            p1 = cdata;
            prev = p1;
            while( vals < *nvals && p1 ) {
                dval = strtod(p1, &p2); /* try to read next value */
                if( p1 == p2 ) break;   /* nothing read, terminate loop */
                prev = p1;              /* store old success ptr */
                p1 = p2;                /* move to next posn */
                ptr[vals] = dval;       /* assign new value  */
                if(xd->verb>6)fprintf(stderr,"  v %f (%f)",ptr[vals],dval);
                vals++;                 /* count new value   */
            }
            if(xd->verb > 6) fputc('\n', stderr);
            break;
        }
        case 512: {     /* NIFTI_TYPE_UINT16 */
            unsigned short * ptr = (unsigned short *)dptr;
            p1 = cdata;
            prev = p1;
            while( vals < *nvals && p1 ) {
                lval = strtol(p1, &p2, 10);   /* try to read next value */
                if( p1 == p2 ) break;   /* nothing read, terminate loop */
                prev = p1;              /* store old success ptr */
                p1 = p2;                /* move to next posn */
                ptr[vals] = lval;       /* assign new value  */
                if(xd->verb>6)fprintf(stderr,"  v %d (%ld)",ptr[vals],lval);
                vals++;                 /* count new value   */
            }
            if(xd->verb > 6) fputc('\n', stderr);
            break;
        }
    }

    /* update the number of values processed */
    if( vals > 0 ) (*nvals) -= vals;

    /* ponder remaining: if *p1 is space, look from there, else from prev */
    if( p1 ){
        if( isspace(*p1) ) {
            remain = cdlen - (p1 - cdata);
            *mod_prev = 0;
        }
        else if( prev ) {
            remain = cdlen - (prev - cdata);
            *mod_prev = 1;  /* still looking at previous val */
        }
    }

    /* if only whitespace left, ignore */
    if( whitespace_len(cdata + (cdlen-remain), remain) == remain )
        remain = 0;

    if(xd->verb > 6) fprintf(stderr,"-d DA: remain = %d\n", remain);

    return remain;
}

static void XMLCALL cb_instr(void *udata, const char *target, const char *data)
{
    gxml_data * xd = (gxml_data *)udata;
    if( xd->verb > 3 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "instr: %s='%s'\n",target,data);
    }
}

static void XMLCALL cb_comment(void *udata, const char * str)
{
    gxml_data * xd = (gxml_data *)udata;
    if( xd->verb > 3 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "comment: '%s'\n",str);
    }
}

static void XMLCALL cb_cdata_start(void *udata)
{
    gxml_data * xd = (gxml_data *)udata;

    if( xd->verb > 3 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "cdata_start\n");
    }
    (void)epush(xd, GXML_ETYPE_CDATA, enames[GXML_ETYPE_CDATA], NULL);
}

static void XMLCALL cb_cdata_end(void *udata)
{
    gxml_data * xd = (gxml_data *)udata;
    epop(xd, GXML_ETYPE_CDATA, enames[GXML_ETYPE_CDATA]);
}

static void XMLCALL cb_default(void *udata, const char * str, int length)
{
    gxml_data * xd = (gxml_data *)udata;
    int wlen = whitespace_len(str,length);
    int len = length;

    if( len == wlen )
    {
        if( xd->verb < 4 ) return;

        str = "whitespace";     /* just note the whitespace */
        len = strlen(str);
    }

    if( xd->verb > 2 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "default [%d]: '%.*s'\n",length,len,str);
    }
}

static void XMLCALL cb_xml_dec(void *udata, const char * ver,
                               const char * enc, int standalone)
{
    gxml_data * xd = (gxml_data *)udata;
    if( xd->verb > 2 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "xmldec ver = %s, enc = %s, standalone = %d\n",
                ver,enc,standalone);
    }
}

static void XMLCALL cb_start_doctype(void *udata, const char * doctype,
                const char * sysid, const char * pubid, int has_subset )
{
    gxml_data * xd = (gxml_data *)udata;
    if( xd->verb > 2 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "start_doctype, dt='%s', sid='%s',pid='%s', sub=%d\n",
               doctype, sysid, pubid, has_subset);
    }
}

static void XMLCALL cb_end_doctype(void *udata)
{
    gxml_data * xd = (gxml_data *)udata;
    if( xd->verb > 2 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr, "end_doctype\n");
    }
}

static void XMLCALL cb_elem_dec(void *udata, const char * ename,
                                             XML_Content * content)
{
    gxml_data * xd = (gxml_data *)udata;
    if( xd->verb > 2 ){
        show_depth(xd->depth, 1, stderr);
        fprintf(stderr,"%s: type=%d, quant=%d, name=%s, numc=%d, cp=%p\n",
                ename, content->type, content->quant, content->name,
                content->numchildren, content->children);
    }
}

static XML_Parser init_xml_parser( void * user_data )
{
    XML_Parser parser;

    parser = XML_ParserCreate(NULL);
    XML_SetUserData(parser, user_data);
    XML_SetStartElementHandler(parser, cb_start_ele);
    XML_SetEndElementHandler(parser, cb_end_ele);
    XML_SetCharacterDataHandler(parser, cb_char);
    XML_SetProcessingInstructionHandler(parser, cb_instr);
    XML_SetCommentHandler(parser, cb_comment);
    XML_SetStartCdataSectionHandler(parser, cb_cdata_start);
    XML_SetEndCdataSectionHandler(parser, cb_cdata_end);
    XML_SetDefaultHandler(parser, cb_default);
    XML_SetXmlDeclHandler(parser, cb_xml_dec);
    XML_SetStartDoctypeDeclHandler(parser, cb_start_doctype);
    XML_SetEndDoctypeDeclHandler(parser, cb_end_doctype);
    XML_SetElementDeclHandler(parser, cb_elem_dec);

    if( GXD.verb > 3 ) fprintf(stderr,"-d parser initialized\n");

    return parser;
}


static int show_stack(char * mesg, gxml_data * xd)
{
    int c;
    if( !xd ) return 1;
    if( mesg ) fputs(mesg, stderr);
    fprintf(stderr,"stack[%d]", xd->depth);
    for( c = 0; c < xd->depth; c++ )
        fprintf(stderr," : %s", enames[xd->stack[c]]);
    fputc('\n', stderr);
    return 0;
}

static int stack_is_valid(gxml_data * xd)
{
    int valid, etype, parent, bad_parent;

    if( xd->depth  < 0 ) return 0;
    if( xd->depth == 0 ) return 1;

    etype = xd->stack[xd->depth-1];         /* depth is at least 1 */

    /* process depth 1 separately, so we can assume a parent later */
    if( xd->depth == 1 ) {
        if( etype != GXML_ETYPE_GIFTI ) {
            show_stack("** invalid element on ", xd);
            return 0;
        }
        return 1;
    }

    /* verify proper parent (or invalid type) */
    valid = 1;          
    bad_parent = 0;
    parent = xd->stack[xd->depth-2];    /* depth is at least 2 */
    switch( etype ) {
        default:
        case GXML_ETYPE_INVALID:
        case GXML_ETYPE_GIFTI:   /* should only be at depth 1 */
            valid = 0;
            break;
        case GXML_ETYPE_META:
            if( parent != GXML_ETYPE_GIFTI &&
                parent != GXML_ETYPE_DATAARRAY )       bad_parent = 1;
            break;
        case GXML_ETYPE_MD:
            if( parent != GXML_ETYPE_META )            bad_parent = 1;
            break;
        case GXML_ETYPE_NAME:
            if( parent != GXML_ETYPE_MD )              bad_parent = 1;
            break;
        case GXML_ETYPE_VALUE:
            if( parent != GXML_ETYPE_MD )              bad_parent = 1;
            break;
        case GXML_ETYPE_LABELTABLE:
            if( parent != GXML_ETYPE_GIFTI )           bad_parent = 1;
            break;
        case GXML_ETYPE_LABEL:
            if( parent != GXML_ETYPE_LABELTABLE )      bad_parent = 1;
            break;
        case GXML_ETYPE_DATAARRAY:
            if( parent != GXML_ETYPE_GIFTI )           bad_parent = 1;
            break;
        case GXML_ETYPE_CSTM:
            if( parent != GXML_ETYPE_DATAARRAY )       bad_parent = 1;
            break;
        case GXML_ETYPE_DATA:
            if( parent != GXML_ETYPE_DATAARRAY )       bad_parent = 1;
            break;
        case GXML_ETYPE_DATASPACE:
            if( parent != GXML_ETYPE_CSTM )            bad_parent = 1;
            break;
        case GXML_ETYPE_XFORMSPACE:
            if( parent != GXML_ETYPE_CSTM )            bad_parent = 1;
            break;
        case GXML_ETYPE_MATRIXDATA:
            if( parent != GXML_ETYPE_CSTM )            bad_parent = 1;
            break;
        case GXML_ETYPE_CDATA:
            if( parent != GXML_ETYPE_NAME      &&
                parent != GXML_ETYPE_VALUE     &&
                parent != GXML_ETYPE_DATASPACE &&
                parent != GXML_ETYPE_XFORMSPACE )      bad_parent = 1;
            break;
    }

    /* possibly print a message if the stack looks bad */
    if( bad_parent && GXD.verb )
        fprintf(stderr,"** %s: bad parent '%s'\n",enames[etype],enames[parent]);
    if( (!valid || bad_parent) && GXD.verb > 1 ) show_stack("** invalid ", xd);

    return valid;
}

/* if bsize is no longer correct, update it and realloc the buffer */
static int reset_xml_buf(gxml_data * xd, char ** buf, int * bsize)
{
    if( *bsize == xd->buf_size ) {
        if( xd->verb > 3 )
            fprintf(stderr,"-d buffer kept at %d bytes\n", *bsize);
        return 0;
    }

    if( xd->verb > 2 )
        fprintf(stderr,"+d update buf, %d to %d bytes\n",*bsize,xd->buf_size);

    *bsize = xd->buf_size;
    *buf = (char *)realloc(*buf, *bsize * sizeof(char));

    if( ! *buf ) {
        fprintf(stderr,"** failed to alloc %d bytes of xml buf!\n", *bsize);
        *bsize = 0;
        return 1;
    }

    return 0;
}

/* decide how big a processing buffer should be
   (either for a small xform matrix or a Data element)
*/
static int partial_buf_size(int nbytes)
{
    if( nbytes <= 64*1024 )      return nbytes;
    if( nbytes <= 10*1024*1024 ) return nbytes / 10;

    return 1024*1024;
}


static int update_partial_buffer(char ** buf, int * blen, int bytes)
{
    int bsize = partial_buf_size(bytes);

    if( !buf || !blen || bytes <= 0 ) {
        fprintf(stderr,"** UPB: bad params (%p,%p,%d)\n", buf, blen, bytes);
        return 1;
    }

    /* just make sure we have a text buffer to work with */
    if( *buf || *blen != bsize ) {
        if( GXD.verb > 2 )
            fprintf(stderr,"+d UPB, alloc %d bytes (from %d) for buffer\n",
                    bsize, bytes);
        *buf = (char *)realloc(*buf, bsize * sizeof(char));
        if( !*buf ) {
            fprintf(stderr,"** UPB: cannot alloc %d bytes for buffer\n",bsize);
            return 1;
        }
        *blen = bsize;
    }

    return 0;
}

static int gxml_write_gifti(gxml_data * xd, FILE * fp)
{
    gifti_image * gim = xd->gim;

    int c, offset;
    int first = 1;  /* first attr to print? */

    if( !gim || !fp ) return 1;

    if( xd->verb > 2 )
        fprintf(stderr,"+d gifti image, numDA = %d, size = %d MB\n",
                gim->numDA, gifti_gim_DA_size(gim,1));

    gxml_write_preamble(xd, fp);
    fprintf(fp,"<%s",enames[GXML_ETYPE_GIFTI]);
    if(gim->version){ fprintf(fp," Version=\"%s\"", gim->version); first = 0; }

    /* add any extra attributes */
    offset = strlen(enames[GXML_ETYPE_GIFTI]) + 2;
    ewrite_ex_atrs(xd, &gim->ex_atrs, offset, first, fp);
    fputs(">\n",fp);

    xd->depth++;
    ewrite_meta(xd, &gim->meta, fp);
    ewrite_LT(xd, &gim->labeltable, fp);

    /* write the DataArray */
    if(!gim->darray) {
        if( xd->verb > 0 )
            fprintf(stderr,"** gifti_image '%s', missing darray\n",xd->ddata);
    } else {
        for( c = 0; c < gim->numDA; c++ )
            ewrite_darray(xd, gim->darray[c], fp);
    }
    
    xd->depth--;
    fprintf(fp,"</%s>\n",enames[GXML_ETYPE_GIFTI]);

    return 0;
}

static int ewrite_darray(gxml_data * xd, DataArray * da, FILE * fp)
{
    int  spaces = xd->skip * xd->depth;
    int  offset, c;
    char dimstr[5] = "Dim0";

    if( xd->verb > 3 ) fprintf(stderr,"+d write DataArray\n");

    if( !da ) return 0;

    offset = strlen(enames[GXML_ETYPE_DATAARRAY]) + 2 + spaces;
    fprintf(fp, "%*s<DataArray", spaces, " ");

    /* print attributes */
    ewrite_str_attr("Category", gifti_category_list[da->category],offset,1,fp);
    ewrite_str_attr("DataType", gifti_datatype2str(da->datatype),offset,0,fp);
    ewrite_str_attr("DataLocation",gifti_dataloc_list[da->location],
                    offset,0,fp);
    ewrite_str_attr("ArrayIndexingOrder",gifti_index_order_list[da->ind_ord],
                    offset,0,fp);
    ewrite_int_attr("Dimensionality", da->num_dim, offset, 0, fp);
    for( c = 0; c < da->num_dim; c++ ) {
        ewrite_int_attr(dimstr, da->dims[c], offset, 0, fp);
        dimstr[3]++;  /* too devious??  iterate '0', '1', ... */
    }
    ewrite_str_attr("Encoding", gifti_encoding_list[da->encoding],offset,0,fp);
    ewrite_str_attr("Endian", gifti_endian_list[da->endian],offset,0,fp);
    fprintf(fp, ">\n");

    /* write sub-elements */
    xd->depth++;
    ewrite_meta(xd, &da->meta, fp);
    ewrite_coordsys(xd, da->coordsys, fp);
    ewrite_data(xd, da, fp);
    xd->depth--;

    fprintf(fp, "%*s</DataArray>\n", spaces, " ");

    return 0;
}


/* this depends on ind_ord, how to write out lines */
static int ewrite_data(gxml_data * xd, DataArray * da, FILE * fp)
{
    int c, spaces = xd->skip * xd->depth;
    int rows, cols;

    if( !da ) return 0;         /* okay, may not exist */

    if( xd->verb > 3 ) fprintf(stderr,"+d write Data\n");
    fprintf(fp, "%*s<%s>\n", spaces, " ", enames[GXML_ETYPE_DATA]);

    if( xd->dstore ) {
        gifti_DA_rows_cols(da, &rows, &cols);  /* product will be nvals */
        for(c = 0; c < rows; c++ )
            ewrite_data_line(da->data, da->datatype,c,cols,spaces+xd->skip,fp);
    }

    fprintf(fp, "%*s</%s>\n", spaces, " ", enames[GXML_ETYPE_DATA]);
    return 0;
}


static int ewrite_coordsys(gxml_data * xd, CoordSystem * cs, FILE * fp)
{
    int c, spaces = xd->skip * xd->depth;

    if( !cs ) return 0;         /* okay, may not exist */

    if( xd->verb > 3 ) fprintf(stderr,"+d write CoordSystem\n");

    fprintf(fp, "%*s<%s>\n", spaces, " ", enames[GXML_ETYPE_CSTM]);
    spaces += xd->skip;

    ewrite_cdata_ele(GXML_ETYPE_DATASPACE, cs->dataspace, spaces, fp);
    ewrite_cdata_ele(GXML_ETYPE_XFORMSPACE, cs->xformspace, spaces,fp);

    fprintf(fp, "%*s<MatrixData>\n", spaces, " ");
    for(c = 0; c < 4; c++ )
        ewrite_double_line(cs->xform[c], 4, spaces+xd->skip, fp);
    fprintf(fp, "%*s</MatrixData>\n", spaces, " ");

    spaces -= xd->skip;
    fprintf(fp, "%*s</%s>\n", spaces, " ", enames[GXML_ETYPE_CSTM]);

    return 0;
}


/* rcr - review format strings */
static int ewrite_data_line(void * data, int type, int row, int cols,
                            int space, FILE * fp)
{
    int c;
    if( !data || row < 0 || cols <= 0 || !fp ) return 1;

    fprintf(fp, "%*s", space, " ");
    switch( type ) {
        default : 
            fprintf(stderr,"** write_data_line, unknown type %d\n",type);
            return -1;
        case 2: {       /* NIFTI_TYPE_UINT8 */
            unsigned char * ptr = (unsigned char *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%u ", ptr[c]);
            break;
        }
        case 4: {       /* NIFTI_TYPE_INT16 */
            short * ptr = (short *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%d ", ptr[c]);
            break;
        }
        case 8: {       /* NIFTI_TYPE_INT32 */
            int * ptr = (int *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%d ", ptr[c]);
            break;
        }
        case 16: {      /* NIFTI_TYPE_FLOAT32 */
            float * ptr = (float *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%f ", ptr[c]);
            break;
        }
        case 32: {      /* NIFTI_TYPE_COMPLEX64 */
            float * ptr = (float *)data + row * cols;
            for(c = 0; c < 2*cols; c+=2)fprintf(fp, "%f %f   ",ptr[c],ptr[c+1]);
            break;
        }
        case 64: {      /* NIFTI_TYPE_FLOAT64 */
            double * ptr = (double *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%f ", ptr[c]);
            break;
        }
        case 128: {     /* NIFTI_TYPE_RGB24 */
            unsigned char * ptr = (unsigned char *)data + row * cols;
            for( c = 0; c < 3*cols; c+=3 )
                fprintf(fp, "%u %u %u   ", ptr[c], ptr[c+1], ptr[c+2]);
            break;
        }
        case 256: {     /* NIFTI_TYPE_INT8 */
            char * ptr = (char *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%d ", ptr[c]);
            break;
        }
        case 512: {     /* NIFTI_TYPE_UINT16 */
            unsigned short * ptr = (unsigned short *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%u ", ptr[c]);
            break;
        }
        case 768: {     /* NIFTI_TYPE_UINT32 */
            unsigned int * ptr = (unsigned int *)data + row * cols;
            for( c = 0; c < cols; c++ ) fprintf(fp, "%u ", ptr[c]);
            break;
        }
        case 1024: {    /* NIFTI_TYPE_INT64 */
            /* rcr - do we need to check #defines? */
            break;
        }
        case 1280: {    /* NIFTI_TYPE_UINT64 */
            /* rcr - do we need to check #defines? */
            break;
        }
        case 1536: {    /* NIFTI_TYPE_FLOAT128 */
            /* rcr - do we need to check #defines? */
            break;
        }
        case 1792: {    /* NIFTI_TYPE_COMPLEX128 */
            double * ptr = (double *)data + row * cols;
            for(c = 0; c < 2*cols; c+=2)fprintf(fp, "%f %f   ",ptr[c],ptr[c+1]);
            break;
        }
        case 2048: {    /* NIFTI_TYPE_COMPLEX256 */
            /* rcr - do we need to check #defines? */
            break;
        }
    }

    fputc('\n', fp);

    return 0;
}


static int ewrite_double_line(double * data, int nvals, int space, FILE * fp)
{
    int c;
    if( !data || nvals <= 0 || !fp ) return 1;

    fprintf(fp, "%*s", space, " ");
    for( c = 0; c < nvals; c++ )        /* duplicate trailing space for diff */
        fprintf(fp, "%f ", data[c]);
    fputc('\n', fp);

    return 0;
}


static int ewrite_cdata_ele(int ele, const char * cdata, int spaces, FILE * fp)
{
    int index = ele;

    if(ele < 0 || ele > GXML_MAX_ELEN) index = 0;  /* be safe */

    fprintf(fp, "%*s<%s><![CDATA[%s]]></%s>\n",
            spaces, " ", enames[index], cdata, enames[index]);

    return 0;
}

static int ewrite_LT(gxml_data * xd, LabelTable * lt, FILE * fp)
{
    int c, spaces = xd->skip * xd->depth;

    if( xd->verb > 3 ) fprintf(stderr,"+d write LabelTable\n");

    if( !lt || lt->length == 0 || !lt->index || !lt->label ) {
        fprintf(fp, "%*s<LabelTable/>\n", spaces, " ");
        return 0;
    }

    fprintf(fp, "%*s<LabelTable>\n", spaces, " ");
    for( c = 0; c < lt->length; c++ ) {
        if( !lt->label[c] ) {
            if(xd->verb > 1) fprintf(stderr,"** label[%d] unset\n", c);
            continue;
        }

        fprintf(fp,"%*s<Label", spaces+xd->skip, " ");
        if( lt->index[c] ) fprintf(fp, " Index='%d'", lt->index[c]);
        fprintf(fp,">%s</Label>\n", lt->label[c]);
    }
    fprintf(fp, "%*s</LabelTable>\n", spaces, " ");

    return 0;
}


static int ewrite_meta(gxml_data * xd, MetaData * md, FILE * fp)
{
    int c, spaces = xd->skip * xd->depth;

    if( xd->verb > 3 ) fprintf(stderr,"+d write MetaData\n");

    if( !md || md->length == 0 || !md->name || !md->value ) {
        fprintf(fp, "%*s<MetaData/>\n", spaces, " ");
        return 0;
    }

    fprintf(fp, "%*s<MetaData>\n", spaces, " ");
    for( c = 0; c < md->length; c++ ) {
        if( !md->name[c] || !md->value[c] ) {
            if(xd->verb > 1) fprintf(stderr,"** MD[%d] unset\n", c);
            continue;
        }

        fprintf(fp,"%*s<MD>\n", spaces+xd->skip, " ");

        fprintf(fp,"%*s<Name><![CDATA[%s]]></Name>\n",
                spaces+2*xd->skip, " ", md->name[c]);
        fprintf(fp,"%*s<Value><![CDATA[%s]]></Value>\n",
                spaces+2*xd->skip, " ", md->value[c]);

        fprintf(fp,"%*s</MD>\n", spaces+xd->skip, " ");
    }
    fprintf(fp, "%*s</MetaData>\n", spaces, " ");

    return 0;
}


/* print a list of attributes, indented to the same level */
static int ewrite_ex_atrs(gxml_data * xd, nvpairs * nvp, int offset,
                          int first, FILE * fp)
{
    int c, spaces = xd->skip * xd->depth + offset;

    if(xd->verb > 2) fprintf(stderr,"+d write %d ex_atr's\n", nvp->length);

    for( c = 0; c < nvp->length; c++ ) {
        ewrite_str_attr(nvp->name[c], nvp->value[c], spaces, first, fp);
        if( first ) first = 0;
    }

    return 0;
}


static int ewrite_int_attr(const char *name, int value, int spaces,
                           int first, FILE * fp)
{
    fprintf(fp, "%s%*s%s=\"%d\"",
            (first) ? "" : "\n",        /* maybe a newline   */
            (first) ?  1 : spaces, " ", /* 1 or many spaces  */
            name, value);
    return 0;
}


static int ewrite_str_attr(const char * name, const char * value, int spaces,
                           int first, FILE * fp)
{
    fprintf(fp, "%s%*s%s=\"%s\"",
            (first) ? "" : "\n",        /* maybe a newline   */
            (first) ?  1 : spaces, " ", /* 1 or many spaces  */
            name, value);
    return 0;
}


static int gxml_write_preamble(gxml_data * xd, FILE * fp)
{
    char version[]  = "1.0";     /* rcr - move to header */
    char encoding[] = "UTF-8";
    char dtd[]      = "http://brainmap.wustl.edu/gifti/gifti.dtd";

    fprintf(fp, "<?xml version=\"%s\" encoding=\"%s\"?>\n", version, encoding);
    fprintf(fp, "<!DOCTYPE GIFTI SYSTEM \"%s\">\n", dtd);

    return 0;
}
