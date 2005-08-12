#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fio.h"
#include "error.h"
#include "tags.h"

int
TAGskip(FILE *fp, int tag, long long len)
{
#if 1
	unsigned char *buf ;
	int ret ;

	buf = (unsigned char *)calloc(len, sizeof(unsigned char)) ;
	ret = fread(buf, sizeof(unsigned char), len, fp) ;
	free(buf) ;
	return(ret) ;
#else
	return(fseek(fp, len, SEEK_CUR)) ;  // doesn't work for gzipped files
#endif
}
int
TAGreadStart(FILE *fp, long long *plen)
{
	int  tag ;
	
	tag = freadInt(fp) ;
	if (feof(fp))
		return(0) ;
	if (tag == TAG_OLD_MGH_XFORM)
	{
		*plen = (long long)freadInt(fp) ;  // sorry - backwards compatibility with Tosa's stuff
		*plen = *plen -1 ; // doesn't include null
	}
	else
		*plen = freadLong(fp) ;
	
	return(tag) ;
}

int
TAGwriteStart(FILE *fp, int tag, long long *phere, long long len)
{
	long here ;

	fwriteInt(tag, fp) ;

	here = ftell(fp) ;
	*phere = (long long)here ;
	fwriteLong(len, fp) ;
	
	return(NO_ERROR) ;
}

int
TAGwrite(FILE *fp, int tag, void *buf, long long len)
{
	long long here ;

	TAGwriteStart(fp, TAG_CMDLINE, &here, len) ;
	fwrite(buf, sizeof(char), len, fp) ;
	TAGwriteEnd(fp, here) ;
	return(NO_ERROR) ;
}

int
TAGwriteEnd(FILE *fp, long long there)
{
	long long here ;

	here = ftell(fp) ;
#if 0
	fseek(fp, there, SEEK_SET) ;
	fwriteLong((long long)(here-(there+sizeof(long long))), fp) ;
	fseek(fp, here, SEEK_SET) ;
#endif

	return(NO_ERROR) ;
}

int
TAGmakeCommandLineString(int argc, char **argv, char *cmd_line)
{
	int i ;

	cmd_line[0] = 0 ;
	for (i = 0 ; i < argc ; i++)
	{
		strcat(cmd_line, argv[i]) ;
		strcat(cmd_line, " ") ;
	}
	return(NO_ERROR) ;
}

int
TAGwriteCommandLine(FILE *fp, char *cmd_line)
{
	long long here ;

	TAGwriteStart(fp, TAG_CMDLINE, &here, strlen(cmd_line)+1) ;
	fprintf(fp, "%s", cmd_line) ;
	TAGwriteEnd(fp, here) ;
	return(NO_ERROR) ;
}

