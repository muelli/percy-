//  Percy++
//  Copyright 2007 Ian Goldberg <iang@cs.uwaterloo.ca>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of version 2 of the GNU General Public License as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  There is a copy of the GNU General Public License in the COPYING file
//  packaged with this plugin; if you cannot find it, write to the Free
//  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
//  02111-1307  USA

#ifndef __RR_ROOTS_H__
#define __RR_ROOTS_H__

#include "ZZ_pXY.h"

NTL_CLIENT

// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  The global ZZ_pContext should already be set to
// p1 * p2, where p2 is prime, and p1 is either prime or 1.
// This routine may also return some spurious values.
vec_ZZ_pX findroots(const ZZ_pXY &P, int degreebound, ZZ p1, ZZ p2);

#endif
