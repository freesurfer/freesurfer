/***************************************************************************
 *   Copyright (C) 2004 by Rudolph Pienaar / Christian Haselgrove          *
 *    Center for Morphometric Analysis       *
 *    Massachusetts General Hospital        *
 *    Building 149, 13th St.         *
 *    Charlestown, MA 02129         *
 *    {ch|rudolph}@nmr.mgh.harvard.edu      *
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

using namespace std;

#ifndef __SURFACE_H__
#define __SURFACE_H__

#ifdef __cplusplus
extern  "C" {
#endif

#include "mri.h"
#include "mrisurf.h"

  void  mark_geodesic(  MRIS   *surf,
                        int   vno_i,
                        int   vno_f,
                        int   mark);


#endif //__SURFACE_H__

#ifdef __cplusplus
}
#endif
