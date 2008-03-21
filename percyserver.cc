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

#include <sys/time.h>
#include <time.h>
#include <vec_ZZ_p.h>
#include "percyserver.h"

NTL_CLIENT

// Handle a request.
void PercyServer::handle_request_GF28(PercyParams &params, istream &is,
	ostream &os)
{
    // Read some values from the params
    unsigned int words_per_block = params.words_per_block();
    unsigned int num_blocks = params.num_blocks();

    // Read the number of queries
    unsigned char nq[2];
    is.read((char *)nq, 2);
    unsigned short num_queries = (nq[0] << 8) | nq[1];

    // For each query, read the input vector, which is a sequence of
    // num_blocks entries, each of length 1 byte.
    GF28_Element inputvector[num_queries*num_blocks];
    GF28_Element outputvector[words_per_block*num_queries];
    memset(outputvector, '\0', words_per_block*num_queries);

    is.read((char *)inputvector, num_queries*num_blocks);

    const FileDataStore* ds = static_cast<FileDataStore*>(&datastore);

    const GF28_Element *data = (const GF28_Element*)(ds->get_data());

    // Compute the output vector and send it back to the client

    if (num_queries > 1) {
#if 0
	for (unsigned int j = 0; j < num_blocks; ++j) {
	    GF28_Element *outp = outputvector;
	    const GF28_Element *inv = inputvector + j;
	    for (unsigned int c = 0; c < words_per_block; ++c) {
		GF28_Element d = data[j * words_per_block + c];
		const GF28_Element *multrow = GF28_mult_table[d];
		const GF28_Element *invq = inv;
		for(int q = 0; q < num_queries; ++q) {
		    outp[q] ^= multrow[*invq];
		    invq += num_blocks;
		}
		outp += num_queries;
	    }
	}
#else
	for (unsigned int j = 0; j < num_blocks; ++j) {
	    GF28_Element *outp = outputvector;
	    const GF28_Element *inv = inputvector + j;
	    const GF28_Element *datap = data+(j * words_per_block);
	    for (unsigned int c = 0; c < words_per_block; ++c) {
		const GF28_Element *multrow = GF28_mult_table[*(datap++)];
		const GF28_Element *invq = inv;
		for(int q = 0; q < num_queries; ++q) {
		    outp[q] ^= multrow[*invq];
		    invq += num_blocks;
		}
		outp += num_queries;
	    }
	}
#endif
    } else {
	//struct timeval ts, te;
	//gettimeofday(&ts, NULL);
	const GF28_Element *block = data;
	for (unsigned int j = 0; j < num_blocks; ++j) {
	    const GF28_Element *multrow = GF28_mult_table[inputvector[j]];
	    const GF28_Element *blockc = block;
	    GF28_Element *oc = outputvector;
	    GF28_Element *oc_end = oc + (words_per_block & ~7);
	    while (oc < oc_end) {
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
		*(oc++) ^= multrow[*(blockc++)];
	    }
	    for (unsigned int c = 0; c < (words_per_block & 7); ++c) {
		*(oc++) ^= multrow[*(blockc++)];
	    }
	    block += words_per_block;
	}
	//gettimeofday(&te, NULL);
	//int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
	//fprintf(stderr, "%d.%3d msec computation\n", td/1000, td%1000);
    }

    if (byzantine) {
	for (unsigned int i=0; i<words_per_block*num_queries; ++i) {
	    outputvector[i]++;
	}
    }

    os.write((char *)outputvector, words_per_block*num_queries);
    os.flush();
}

// Handle a request.
void PercyServer::handle_request(PercyParams &params, istream &is,
	ostream &os)
{
    if (params.is_gf28()) {
	handle_request_GF28(params, is, os);
	return;
    }

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

    // Read the number of queries
    unsigned char nq[2];
    is.read((char *)nq, 2);
    unsigned short num_queries = (nq[0] << 8) | nq[1];

    // For each query, read the input vector, which is a sequence of
    // num_blocks entries, each of length modulus_bytes.
    vec_ZZ_p inputvector[num_queries];

    for(int q = 0; q < num_queries; ++q) {
	inputvector[q].SetLength(num_blocks);

	for(unsigned int i = 0; i < num_blocks; ++i) {
	    unsigned char inputdata[modulus_bytes];
	    memset(inputdata, 0, modulus_bytes);
	    ZZ inputz;
	    is.read((char *)inputdata, modulus_bytes);
	    ZZFromBytes(inputz, inputdata, modulus_bytes);
	    inputvector[q][i] = to_ZZ_p(inputz);
	}
    }

    // Compute the output vector and send it back to the client

    for (unsigned int c = 0; c < words_per_block; ++c) {
	ZZ_p reswrds[num_queries];
	compute_one(reswrds, hybrid_protection, num_blocks,
		num_queries, inputvector, c);
	unsigned char bytes[modulus_bytes];
	for(int q = 0; q < num_queries; ++q) {
	    BytesFromZZ(bytes, rep(reswrds[q]), modulus_bytes);
	    os.write((char *)bytes, modulus_bytes);
	}
    }
    os.flush();
}

void PercyServer::compute_one(ZZ_p *value, bool hybrid_protection,
	unsigned int num_blocks, unsigned short num_queries,
	const vec_ZZ_p *inputvector, unsigned int c)
{
    for(int q = 0; q < num_queries; ++q) {
	value[q] = hybrid_protection ? 1 : 0;
    }

    for (unsigned int j = 0; j < num_blocks; ++j) {
	// The cth word of the jth block
	ZZ wrd = datastore.get_word(c, j);

	if (hybrid_protection) {
	    for(int q = 0; q < num_queries; ++q) {
		value[q] = value[q] * power(inputvector[q][j], wrd);
	    }
	} else {
	    for(int q = 0; q < num_queries; ++q) {
		value[q] = value[q] + inputvector[q][j] * to_ZZ_p(wrd);
	    }
	}
    }

    if (byzantine) {
	// Produce a *consistent* incorrect value for maximal client
	// confusion.
	for(int q = 0; q < num_queries; ++q) {
	    value[q] += 1;
	}
    }
}
