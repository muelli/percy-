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

#include <vec_ZZ_p.h>
#include "percyserver.h"

NTL_CLIENT

// Handle a request.
void PercyServer::handle_request(PercyParams &params, istream &is,
	ostream &os)
{
    // Set the appropriate modulus
    unsigned int modulus_bytes;
    if (params.hybrid()) {
	params.mod_modulussq();
	modulus_bytes = params.modulussq_bytes();
    } else {
	params.mod_modulus();
	modulus_bytes = params.modulus_bytes();
    }

    // Read some values from the params
    bool hybrid_protection = params.hybrid();
    unsigned int words_per_block = params.words_per_block();
    unsigned int num_blocks = params.num_blocks();

    // Read the input vector, which is a sequence of num_blocks entries,
    // each of length modulus_bytes.
    vec_ZZ_p inputvector;

    inputvector.SetLength(num_blocks);

    for(unsigned int i = 0; i < num_blocks; ++i) {
	unsigned char inputdata[modulus_bytes];
	memset(inputdata, 0, modulus_bytes);
	ZZ inputz;
	is.read((char *)inputdata, modulus_bytes);
	ZZFromBytes(inputz, inputdata, modulus_bytes);
	inputvector[i] = to_ZZ_p(inputz);
    }

    // Compute the output vector and send it back to the client

    for (unsigned int c = 0; c < words_per_block; ++c) {
	ZZ_p reswrd = compute_one(hybrid_protection, num_blocks,
		inputvector, c);
	unsigned char bytes[modulus_bytes];
	BytesFromZZ(bytes, rep(reswrd), modulus_bytes);
	os.write((char *)bytes, modulus_bytes);
    }
    os.flush();
}

ZZ_p PercyServer::compute_one(bool hybrid_protection, unsigned int num_blocks,
	const vec_ZZ_p &inputvector, unsigned int c)
{
    ZZ_p value;

    value = hybrid_protection ? 1 : 0;

    for (unsigned int j = 0; j < num_blocks; ++j) {
	// The cth word of the jth block
	ZZ wrd = datastore.get_word(c, j);

	if (hybrid_protection) {
	    value = value * power(inputvector[j], wrd);
	} else {
	    value = value + inputvector[j] * to_ZZ_p(wrd);
	}
    }

    if (byzantine) {
	// Produce a *consistent* incorrect value for maximal client
	// confusion.
	value += 1;
    }
    return value;
}
