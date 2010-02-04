#!/bin/bash
echo "const char *$1Str = \"\\" > $1.cl.h
cat $1.cl |  sed 's/\( *\)$//' | sed 's/$/ \\n\\/' >> $1.cl.h
echo "\";" >> $1.cl.h
