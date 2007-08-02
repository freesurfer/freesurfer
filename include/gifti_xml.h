/* Author: Rick Reynolds reynoldr@mail.nih.gov */

#ifndef _GIFTI_XML_H_
#define _GIFTI_XML_H_

#define GXML_MAX_DEPTH 10    /* maximum stack depth */
#define GXML_MAX_ELEN  128   /* maximum element length */

/* ----------------------------------------------------------------------
   element      depths  parent(s)       children
   -------      ------  --------------  -----------------------
   GIFTI        0                       MetaData, LabelTable, DataArray
   MetaData     1       GIFTI           MD
                2       DataArray
   MD           2,+1    MetaData        Name, Value
   Name         3,+1    MD              CDATA/char
   Value        3,+1    MD              CDATA/char
   LabelTable   1       GIFTI           Label
   Label        2       LabelTable      CDATA/char

   DataArray    1       GIFTI           MetaData, CSTM, Data
   CSTM         2       DataArray       DataSpace, TransformedSpace, MatrixData
   Data         2       DataArray
   DataSpace    3       CSTM            CDATA/char
   TransformedSpace  3  CSTM            CDATA/char
   MatrixData   3       CSTM            char


   CDATA        4,+1    Name, Value     char
                4       DataSpace       char
                4       TransformedSpace char
   char         any     any             whitespace
                5       CDATA

   -- other objects to handle --
   XML declaration:     version, encoding, standalone
   DOCTYPE:             type=GIFTI, sid=.../gifti.dtd, pid, sub
   default:
   ----------------------------------------------------------------------
*/

/* this list must match enames, and is ordered via the above comment */
#define GXML_ETYPE_INVALID      0
#define GXML_ETYPE_GIFTI        1      /* GIFTI element            */
#define GXML_ETYPE_META         2      /* MetaData element         */
#define GXML_ETYPE_MD           3      /* MD element               */
#define GXML_ETYPE_NAME         4      /* Name element             */
#define GXML_ETYPE_VALUE        5      /* Value element            */
#define GXML_ETYPE_LABELTABLE   6      /* LabelTable element       */
#define GXML_ETYPE_LABEL        7      /* Label element            */
#define GXML_ETYPE_DATAARRAY    8      /* DataArray element        */
#define GXML_ETYPE_CSTM         9      /* CSTM element             */
#define GXML_ETYPE_DATA        10      /* Data element             */
#define GXML_ETYPE_DATASPACE   11      /* DataSpace element        */
#define GXML_ETYPE_XFORMSPACE  12      /* TransformedSpace element */
#define GXML_ETYPE_MATRIXDATA  13      /* MatrixData element       */
#define GXML_ETYPE_CDATA       14      /* CDATA element            */
#define GXML_ETYPE_LAST        14      /* should match last entry  */


typedef struct {
    int            verb;            /* verbose level                */
    int            dstore;          /* flag: store data             */
    int            buf_size;        /* for XML buffer               */
    int            errors;          /* number of errors encountered */
    int            skip;            /* stack depth to skip          */
    int            depth;           /* current stack depth          */
    int            stack[GXML_MAX_DEPTH+1]; /* stack of etypes      */

    int            dind;            /* index into decode array      */
    int            doff;            /* offset into data buffer      */
    int            clen;            /* length of current CDATA      */
    int            xlen;            /* length of xform buffer       */
    int            dlen;            /* length of Data buffer        */
    char        ** cdata;           /* pointer to current CDATA     */
    char         * xdata;           /* xform buffer (free at end)   */
    char         * ddata;           /* Data buffer (free at end)    */
    gifti_image  * gim;             /* pointer to returning image   */
} gxml_data;

/* protos */

/* main interface */
gifti_image * gxml_read_image (const char * fname, int read_data);
int           gxml_write_image(gifti_image * gim, const char * fname,
                               int write_data);

int   gxml_set_verb        ( int val );
int   gxml_get_verb        ( void    );
int   gxml_set_dstore      ( int val );
int   gxml_get_dstore      ( void    );
int   gxml_set_buf_size    ( int val );
int   gxml_get_buf_size    ( void    );


#endif /* _GIFTI_XML_H_ */
