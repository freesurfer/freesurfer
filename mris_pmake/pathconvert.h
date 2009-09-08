/***************************************************************************
 *   Copyright (C) 2006 by Rudolph Pienaar / Christian Haselgrove          *
 *   {ch|rudolph}@nmr.mgh.harvard.edu                                      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
 NAME

 pathconvert.h

 DESCRIPTION

 'pathconvert.h' declares two functions, 'abs2rel' and 'rel2abs'
 that convert directory path specifications in *nix.

 These functions were originally defined by Shigio Yamaguchi are
 copyright 1999 Tama Communications Corporation and are released
 under the LGPL.

*/

#ifndef __PATHCONVERT_H__
#define __PATHCONVERT_H__

#ifdef __cplusplus
extern  "C" {
#endif


#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>


  char    *abs2rel __P((const char *, const char *, char *, size_t));
  char    *rel2abs __P((const char *, const char *, char *, size_t));

#ifdef __cplusplus
}
#endif


#endif // __PATHCONVERT_H__
