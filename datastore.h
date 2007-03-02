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

#ifndef __DATASTORE_H__
#define __DATASTORE_H__
#include <ZZ.h>
#include "percyparams.h"

NTL_CLIENT

class DataStore {
public:
    virtual ~DataStore() = 0;

    // Retrieve the cth word of the jth block of the database as a ZZ
    // Both c and j are 0-based.
    virtual ZZ get_word(unsigned int c, unsigned int j) = 0;
};

class FileDataStore : public DataStore {
protected:
    int dbfd;
    unsigned char *database;
    size_t totbytes;
    off_t offset;
    unsigned int bytes_per_word, words_per_block, num_blocks;

public:
    // Create a DataStore backed by a file.
    FileDataStore(const char *filename, const PercyParams &params,
	    bool tau_independent = false);
    ~FileDataStore();
 
    // Retrieve the cth word of the jth block of the database as a ZZ
    // Both c and j are 0-based.
    ZZ get_word(unsigned int c, unsigned int j);
};

#endif
