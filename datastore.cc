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

#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include "datastore.h"
#include "percyio.h"

// You have to define even pure virtual destructors in C++.  You just do.
DataStore::~DataStore() {}

FileDataStore::FileDataStore(const char *filename, const PercyParams &params,
	bool tau_independent)
{
    offset = 0;

    if (tau_independent) {
	// Process the PIRD header
	ifstream dbhead(filename);

	// We'll turn these into exceptions later
	char headbuf[6];
	dbhead.read(headbuf, 6);
	if (memcmp(headbuf, "PIRD\x01\x00", 6)) {
	    std::cerr << "Split database not in PIRD format.\n";
	    exit(1);
	}
	ZZ dbmodulus;
	percy_read_ZZ(dbhead, dbmodulus);

	if (!params.modulus_match(dbmodulus)) {
	    std::cerr << "Incorrect modulus for split database.\n";
	    exit(1);
	}

	offset = dbhead.tellg();
	dbhead.close();
    }

    unsigned int bytes_per_word = tau_independent ? params.modulus_bytes()
	: params.bytes_per_word();
    unsigned int words_per_block = params.words_per_block();
    unsigned int num_blocks = params.num_blocks();

    // Open the file so we can mmap it
    dbfd = open(filename, O_RDONLY);

    if (dbfd < 0) {
	perror("open");
	exit(1);
    }
    totbytes = bytes_per_word * words_per_block * num_blocks;
    struct stat st;
    fstat(dbfd, &st);
    if (st.st_size < (off_t)totbytes+offset) {
        fprintf(stderr, "Database too small!\n");
        exit(1);
    }
    database = (unsigned char *)mmap(NULL, totbytes+offset, PROT_READ,
	    MAP_SHARED, dbfd, 0);
    if (database == MAP_FAILED) {
	perror("mmap");
	exit(1);
    }
    
    this->bytes_per_word = bytes_per_word;
    this->words_per_block = words_per_block;
    this->num_blocks = num_blocks;
}

FileDataStore::~FileDataStore()
{
    munmap(database, totbytes);
    close(dbfd);
}

// Retrieve the cth word of the jth block of the database as a ZZ.
// Both c and j are 0-based.
ZZ FileDataStore::get_word(unsigned int c, unsigned int j)
{
    ZZ ret;
    ZZFromBytes(ret, database + offset +
	    (j*words_per_block+c)*(bytes_per_word), bytes_per_word);
    return ret;
}
