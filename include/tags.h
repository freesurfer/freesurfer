#ifndef TAGS_H
#define TAGS_H

#define TAG_COLORTABLE         1
#define TAG_USEREALRAS         2
#define TAG_CMDLINE            3

#define TAG_GCAMORPH_GEOM      10
#define TAG_GCAMORPH_TYPE      11
#define TAG_GCAMORPH_LABELS    12

#define TAG_SURF_GEOM          20

#define TAG_MGH_XFORM          30

int TAGreadStart(FILE *fp, long long *plen) ;
int TAGwriteStart(FILE *fp, int tag, long long *phere, long long len) ;
int TAGwriteEnd(FILE *fp, long long there) ;
int TAGskip(FILE *fp, int tag, long long len) ;
int TAGmakeCommandLineString(int argc, char **argv, char *cmd_line) ;
int TAGwriteCommandLine(FILE *fp, char *cmd_line) ;
int TAGwrite(FILE *fp, int tag, void *buf, long long len) ;
 
#endif
