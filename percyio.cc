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

#include "percyio.h"

// Read and write little-endian 4-byte unsigned ints, 2-byte unsigned shorts,
// and ZZs.
void percy_write_le_int(ostream &os, unsigned int val)
{
    unsigned char buf[4];

    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;
    buf[2] = (val >> 16) & 0xff;
    buf[3] = (val >> 24) & 0xff;

    os.write((char *)buf, 4);
}

void percy_write_le_short(ostream &os, unsigned short val)
{
    unsigned char buf[2];

    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;

    os.write((char *)buf, 2);
}

void percy_write_ZZ(ostream &os, ZZ val)
{ 
    // Output the length
    long len = NumBytes(val);
    if (len > 65535) {
	std::cerr << "ZZ is too long (" << len << " bytes)\n";
	len = 65535;
    }
    unsigned char *buf = new unsigned char[len];
    BytesFromZZ(buf, val, len);
    percy_write_le_short(os, len);
    os.write((char *)buf, len);
    delete[] buf;
}

void percy_read_le_int(istream &is, unsigned int &val)
{
    unsigned char buf[4];

    is.read((char *)buf, 4);

    val = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
}

void percy_read_le_short(istream &is, unsigned short &val)
{
    unsigned char buf[2];

    is.read((char *)buf, 2);

    val = buf[0] | (buf[1] << 8);
}

void percy_read_ZZ(istream &is, ZZ &val)
{
    // Input the length
    unsigned short len;
    percy_read_le_short(is, len);
    unsigned char *buf = new unsigned char[len];
    is.read((char *)buf, len);
    ZZFromBytes(val, buf, len);
    delete[] buf;
}
