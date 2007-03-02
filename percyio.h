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

#ifndef __PERCYIO_H__
#define __PERCYIO_H__

#include <ZZ.h>

NTL_CLIENT

// Read and write little-endian 4-byte unsigned ints, 2-byte unsigned shorts,
// and ZZs.
void percy_write_le_int(ostream &os, unsigned int val);
void percy_write_le_short(ostream &os, unsigned short val);
void percy_write_ZZ(ostream &os, ZZ val);
void percy_read_le_int(istream &is, unsigned int &val);
void percy_read_le_short(istream &is, unsigned short &val);
void percy_read_ZZ(istream &is, ZZ &val);

#endif
