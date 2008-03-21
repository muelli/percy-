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

#include <iostream>
#include <fstream>
#include <vec_vec_ZZ_p.h>
#include <ZZ_pX.h>
#include "recover.h"
#include "percyclient.h"

// Generate t-private (t+1)-of-l shares of a given secret value.
static void genshares(unsigned short t, unsigned short l,
    const vec_ZZ_p &indices, vec_ZZ_p &values, unsigned int secret)
{
    // Pick a random polynomial of degree t
    ZZ_pX randpoly = random_ZZ_pX(t+1);
    // Set the constant term to the secret
    SetCoeff(randpoly, 0, secret);
    // std::cerr << "Poly is (" << randpoly << ")\n";

    // Evaluate the polynomial at each of the indices
    for (unsigned int i=0; i<l; ++i) {
	eval(values[i], randpoly, indices[i]);
	// std::cerr << "(" << indices[i] << ", " << values[i] << ")\n\n";
    }
}

// Generate t-private (t+1)-of-l shares of a given secret value.
static void genshares_GF28(unsigned short t, unsigned short l,
    const GF28_Element *indices, GF28_Element *values, GF28_Element secret)
{
    int i;

    // Pick a random polynomial of degree t with the right constant term
    GF28_Element coeffs[t+1];
    coeffs[0] = secret;
    for (i=1;i<=t;++i) {
	coeffs[i] = RandomBits_ulong(8);
    }

    // Evaluate the polynomial at each of the indices
    for (unsigned int i=0; i<l; ++i) {
	values[i] = evalpoly_GF28(coeffs, t, indices[i]);
	// std::cerr << "(" << indices[i] << ", " << values[i] << ")\n\n";
    }
}

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
int PercyClient::send_request_GF28(vector<unsigned int> block_numbers,
	std::vector<iostream*> &iosvec)
{
    if (num_servers != iosvec.size()) {
	std::cerr << "Incorrect iostream vector size passed to "
	    "send_request_GF28.\n";
	std::cerr << "Was " << iosvec.size() << ", should be " << num_servers
	    << ".\n";
	return -1;
    }

    // Construct the vector of server indices
    num_queries = block_numbers.size();
    int q;
    indices_gf28 = new GF28_Element[num_queries*num_servers];
    for (unsigned int j=0; j<num_servers; ++j) {
	if (params.tau() > 0) {
	    // Use the constant indices 1, 2, ..., num_servers.
	    for (q=0; q<num_queries; ++q) {
		indices_gf28[q*num_servers+j] = j+1;
	    }
	} else {
	    // Use random indices
	    GF28_Element r;
	    bool ok = false;
	    do {
		r = RandomLen_long(8);
		if (r != 0) {
		    ok = true;
		    for (unsigned int k=0;k<j;++k) {
			if (indices_gf28[k] == r) {
			    ok = false;
			}
		    }
		}
	    } while (!ok);

	    for (q=0; q<num_queries; ++q) {
		indices_gf28[q*num_servers+j] = r;
	    }
	}
    }

    // Construct the shares of the e_{index} vector
    unsigned int num_blocks = params.num_blocks();
    GF28_Element shares[num_queries][num_blocks][num_servers];

    for (q=0; q<num_queries; ++q) {
	for (unsigned int i = 0; i < num_blocks; ++i) {
	    genshares_GF28(t, num_servers, indices_gf28 + q*num_servers,
		    shares[q][i], i == block_numbers[q]);
	}
    }

    // Send the params and query to each server
    for (unsigned int j = 0; j < num_servers; ++j) {
	*(iosvec[j]) << params;
	unsigned char nq[2];
	nq[0] = (q >> 8) & 0xff;
	nq[1] = (q) & 0xff;
	iosvec[j]->write((char *)nq, 2);

	for (q=0; q<num_queries; ++q) {
	    for (unsigned int i = 0; i < num_blocks; ++i) {
		char *sharebyte = (char *)&(shares[q][i][j]);
		iosvec[j]->write(sharebyte, 1);
	    }
	}
	iosvec[j]->flush();
    }

    return 0;
}

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
int PercyClient::send_request(vector<unsigned int> block_numbers,
	std::vector<iostream*> &iosvec)
{
    if (params.is_gf28()) {
	return send_request_GF28(block_numbers, iosvec);
    }

    if (num_servers != iosvec.size()) {
	std::cerr << "Incorrect iostream vector size passed to "
	    "send_request.\n";
	std::cerr << "Was " << iosvec.size() << ", should be " << num_servers
	    << ".\n";
	return -1;
    }

    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();

    // Construct the vector of server indices
    params.mod_modulus();
    num_queries = block_numbers.size();
    int q;
    indicesp = new vec_ZZ_p[num_queries];
    for (q=0; q<num_queries; ++q) {
	indicesp[q].SetLength(num_servers);
    }
    for (unsigned int j=0; j<num_servers; ++j) {
	if (params.tau() > 0) {
	    // Use the constant indices 1, 2, ..., num_servers.
	    for (q=0; q<num_queries; ++q) {
		indicesp[q][j] = j+1;
	    }
	} else {
	    // Use random indices
	    ZZ_p r = random_ZZ_p();
	    for (q=0; q<num_queries; ++q) {
		indicesp[q][j] = r;
	    }
	}
    }

    // Construct the shares of the e_{index} vector
    vec_vec_ZZ_p shares[num_queries];
    unsigned int num_blocks = params.num_blocks();

    for (q=0; q<num_queries; ++q) {
	shares[q].SetLength(num_blocks);
	for (unsigned int i = 0; i < num_blocks; ++i) {
	    params.mod_modulussq();
	    shares[q][i].SetLength(num_servers);
	    params.mod_modulus();
	    genshares(t, num_servers, indicesp[q], shares[q][i],
		    i == block_numbers[q]);
	}
    }

    // Optionally encrypt the shares
    if (params.hybrid()) {
	for (unsigned int i = 0; i < num_blocks; ++i) {
	    for (unsigned int j = 0; j < num_servers; ++j) {
		for (q=0; q<num_queries; ++q) {
		    params.mod_modulus();
		    ZZ share = rep(shares[q][i][j]);
		    params.mod_modulussq();
		    shares[q][i][j] = params.encrypt(share);
		}
		std::cerr << i << " " << j << "\n";
	    }
	}
	// If we're encrypting, leave modulus^2 as the active modulus
	params.mod_modulussq();
    }

    // Send the params and query to each server
    unsigned int modulus_bytes = params.hybrid() ? params.modulussq_bytes() :
	params.modulus_bytes();
    for (unsigned int j = 0; j < num_servers; ++j) {
	*(iosvec[j]) << params;
	unsigned char nq[2];
	nq[0] = (q >> 8) & 0xff;
	nq[1] = (q) & 0xff;
	iosvec[j]->write((char *)nq, 2);

	for (q=0; q<num_queries; ++q) {
	    for (unsigned int i = 0; i < num_blocks; ++i) {
		unsigned char bytes[modulus_bytes];
		BytesFromZZ(bytes, rep(shares[q][i][j]), modulus_bytes);
		iosvec[j]->write((char *)bytes, modulus_bytes);
	    }
	}
	iosvec[j]->flush();
    }

    savectx.restore();
    return 0;
}

// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
unsigned short PercyClient::receive_replies_GF28(std::vector<iostream*> &iosvec)
{
    unsigned int words_per_block = params.words_per_block();
    int q;

    // The vector of servers that have responded properly
    goodservers.clear();

    // The responses from the servers
    answers_gf28 = new GF28_Element[num_queries*words_per_block*num_servers];

    // Read the replies
    for (unsigned int j = 0; j < num_servers; ++j) {
	bool isgood = true;
	for (unsigned int i = 0; isgood && i < words_per_block; ++i) {
	    for (q=0; q<num_queries; ++q) {
		iosvec[j]->read((char *)(answers_gf28 + 
			    (q*words_per_block + i)*num_servers + j), 1);
		if ((unsigned int)(iosvec[j]->gcount()) < 1) {
		    // Mark this server as bad
		    std::cerr << "Marking server " << j+1 << " as bad.\n";
		    isgood = false;
		    break;
		}
	    }
	}
	if (isgood) {
	    goodservers.push_back(j);
	}
    }

    return goodservers.size();
}

// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
unsigned short PercyClient::receive_replies(std::vector<iostream*> &iosvec)
{
    if (params.is_gf28()) {
	return receive_replies_GF28(iosvec);
    }

    unsigned int words_per_block = params.words_per_block();
    int q;

    // Choose the right modulus
    unsigned int modulus_bytes;
    if (params.hybrid()) {
	modulus_bytes = params.modulussq_bytes();
	params.mod_modulussq();
    } else {
	modulus_bytes = params.modulus_bytes();
	params.mod_modulus();
    }

    // The vector of servers that have responded properly
    goodservers.clear();

    // The responses from the servers
    answersp = new vec_vec_ZZ_p[num_queries];
    for (q=0; q<num_queries; ++q) {
	answersp[q].SetLength(words_per_block);
	for (unsigned long c = 0; c < words_per_block; ++c) {
	    answersp[q][c].SetLength(num_servers);
	}
    }

    // Read the replies
    for (unsigned int j = 0; j < num_servers; ++j) {
	bool isgood = true;
	for (unsigned int i = 0; isgood && i < words_per_block; ++i) {
	    for (q=0; q<num_queries; ++q) {
		unsigned char bytes[modulus_bytes];
		iosvec[j]->read((char *)bytes, modulus_bytes);
		if ((unsigned int)(iosvec[j]->gcount()) < modulus_bytes) {
		    // Mark this server as bad
		    std::cerr << "Marking server " << j+1 << " as bad.\n";
		    isgood = false;
		    break;
		}
		ZZ ans;
		ZZFromBytes(ans, bytes, modulus_bytes);
		answersp[q][i][j] = to_ZZ_p(ans);
	    }
	}
	if (isgood) {
	    goodservers.push_back(j);
	}
    }

    // Optionally decrypt the answers
    if (params.hybrid()) {
	params.mod_modulussq();
	for (unsigned long i = 0; i < words_per_block; ++i) {
	    for (unsigned long j = 0; j < num_servers; ++j) {
		for (q=0; q<num_queries; ++q) {
		    ZZ_p dec = params.decrypt(answersp[q][i][j]);
		    clear(answersp[q][i][j]);
		    answersp[q][i][j] = dec;
		}
	    }
	}
    }

    // Now we're mod modulus for sure
    params.mod_modulus();

    return goodservers.size();
}

vector< vector<PercyResult> > PercyClient::process_replies_GF28(
	unsigned short h)
{
    unsigned int words_per_block = params.words_per_block();
    unsigned short tau = params.tau();

    std::cerr << goodservers.size() << " of " << num_servers << " servers responded.\n";

    int q;
    vector< vector<PercyResult> > allH;

    for (q=0; q<num_queries; ++q) {
	vector<PercyResult> H, Hprime;
	string empty;
	PercyResult n(goodservers, empty);
	H.push_back(n);

	for (unsigned long i = 0; i < words_per_block; ++i) {
	    if (goodservers.size() <= (t+tau)) {
		std::cerr << "Too few honest servers to recover data!\n";
		allH.clear();
		return allH;
	    }

	    // std::cerr << i << " " << t+tau << " " << numgood << " " << h << "\n";
	    Hprime = EasyRecover_GF28(t+tau, h, H,
		    answers_gf28 + (q*words_per_block + i)*num_servers,
		    indices_gf28 + q*num_servers);
	    // std::cerr << "easyrecover returned " << res << "\n";
	    if (Hprime.empty()) {
		// We seem to have some Byzantine servers.  Try to identify
		// them with the expensive algorithm
		std::cerr << "Switching to HardRecover...\n";
		GSDecoder_GF2E decoder;

		// Convert the raw bytes to vec_GF2E
		vec_GF2E indices_vec;
		indices_vec.SetLength(num_servers);
		for (int ix=0; ix<num_servers; ++ix) {
		    conv(indices_vec[ix],
			    GF2XFromBytes(indices_gf28 + q*num_servers + ix,
				1));
		}
		vec_GF2E answers_vec;
		answers_vec.SetLength(num_servers);
		for (int ix=0; ix<num_servers; ++ix) {
		    conv(answers_vec[ix],
			    GF2XFromBytes(answers_gf28 +
				(q*words_per_block + i)*num_servers + ix,
				1));
		}
		Hprime = decoder.HardRecover(1, t+tau, h, H,
			goodservers, answers_vec, indices_vec);
	    }
	    if (Hprime.empty()) {
		// Unable to recover data
		std::cerr << "Too few honest servers to recover data!\n";
		allH.clear();
		return allH;
	    }
	    H = Hprime;
	    // std::cerr << i << " " << t << " " << numgood << " " << x << "\n";
	}
	allH.push_back(H);
    }
    return allH;
}

vector< vector<PercyResult> > PercyClient::process_replies(unsigned short h,
	ZZ p1, ZZ p2)
{
    if (params.is_gf28()) {
	return process_replies_GF28(h);
    }

    unsigned int words_per_block = params.words_per_block();
    unsigned int bytes_per_word = params.bytes_per_word();
    unsigned short tau = params.tau();

    std::cerr << goodservers.size() << " of " << num_servers << " servers responded.\n";

    int q;
    vector< vector<PercyResult> > allH;

    for (q=0; q<num_queries; ++q) {
	vector<PercyResult> H, Hprime;
	string empty;
	PercyResult n(goodservers, empty);
	H.push_back(n);

	for (unsigned long i = 0; i < words_per_block; ++i) {
	    if (goodservers.size() <= (t+tau)) {
		std::cerr << "Too few honest servers to recover data!\n";
		allH.clear();
		return allH;
	    }

	    // std::cerr << i << " " << t+tau << " " << numgood << " " << h << "\n";
	    Hprime = EasyRecover(bytes_per_word, t+tau, h, H, answersp[q][i], indicesp[q]);
	    // std::cerr << "easyrecover returned " << res << "\n";
	    if (Hprime.empty()) {
		// We seem to have some Byzantine servers.  Try to identify
		// them with the expensive algorithm
		std::cerr << "Switching to HardRecover...\n";
		GSDecoder_ZZ_p decoder(p1, p2);
		Hprime = decoder.HardRecover(bytes_per_word, t+tau, h, H,
			goodservers, answersp[q][i], indicesp[q]);
	    }
	    if (Hprime.empty()) {
		// Still unable to recover data
		std::cerr << "Too few honest servers to recover data!\n";
		allH.clear();
		return allH;
	    }
	    H = Hprime;
	    // std::cerr << i << " " << t << " " << numgood << " " << x << "\n";
	}
	allH.push_back(H);
    }
    return allH;
}

// Do all of the above in one shot
vector< vector<PercyResult> > PercyClient::fetch_blocks(
	vector<unsigned int> block_numbers, std::vector<iostream*> &iosvec,
	ZZ p1, ZZ p2)
{
    int res = send_request(block_numbers, iosvec);
    if (res < 0) {
	vector< vector<PercyResult> > empty;
	return empty;
    }
    unsigned short k = receive_replies(iosvec);
    // Calculate what h to use
    unsigned short h = (unsigned short)(floor(sqrt((t+params.tau())*k)))+1;
    return process_replies(h, p1, p2);
}
