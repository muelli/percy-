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

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
int PercyClient::send_request(unsigned int block_number,
	std::vector<iostream*> &iosvec)
{
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
    indices.SetLength(num_servers);
    for (unsigned int j=0; j<num_servers; ++j) {
	if (params.tau() > 0) {
	    // Use the constant indices 1, 2, ..., num_servers.
	    indices[j] = j+1;
	} else {
	    // Use random indices
	    indices[j] = random_ZZ_p();
	}
    }

    // Construct the shares of the e_{index} vector
    vec_vec_ZZ_p shares;
    unsigned int num_blocks = params.num_blocks();
    shares.SetLength(num_blocks);
    for (unsigned int i = 0; i < num_blocks; ++i) {
	params.mod_modulussq();
	shares[i].SetLength(num_servers);
	params.mod_modulus();
	genshares(t, num_servers, indices, shares[i], i == block_number);
    }

    // Optionally encrypt the shares
    if (params.hybrid()) {
	for (unsigned int i = 0; i < num_blocks; ++i) {
	    for (unsigned int j = 0; j < num_servers; ++j) {
		params.mod_modulus();
		ZZ share = rep(shares[i][j]);
		params.mod_modulussq();
		shares[i][j] = params.encrypt(share);
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
	for (unsigned int i = 0; i < num_blocks; ++i) {
	    unsigned char bytes[modulus_bytes];
	    BytesFromZZ(bytes, rep(shares[i][j]), modulus_bytes);
	    iosvec[j]->write((char *)bytes, modulus_bytes);
	}
	iosvec[j]->flush();
    }

    savectx.restore();
    return 0;
}

// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
unsigned short PercyClient::receive_replies(std::vector<iostream*> &iosvec)
{
    unsigned int words_per_block = params.words_per_block();

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
    answers.SetLength(words_per_block);
    for (unsigned long c = 0; c < words_per_block; ++c) {
	answers[c].SetLength(num_servers);
    }

    // Read the replies
    for (unsigned int j = 0; j < num_servers; ++j) {
	bool isgood = true;
	for (unsigned int i = 0; i < words_per_block; ++i) {
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
	    answers[i][j] = to_ZZ_p(ans);
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
		ZZ_p dec = params.decrypt(answers[i][j]);
		clear(answers[i][j]);
		answers[i][j] = dec;
	    }
	}
    }

    // Now we're mod modulus for sure
    params.mod_modulus();

    return goodservers.size();
}

vector<PercyResult> PercyClient::process_replies(unsigned short h,
	ZZ p1, ZZ p2)
{
    unsigned int words_per_block = params.words_per_block();
    unsigned int bytes_per_word = params.bytes_per_word();
    unsigned short tau = params.tau();

    std::cerr << goodservers.size() << " of " << num_servers << " servers responded.\n";

    vector<PercyResult> H, Hprime;
    string empty;
    PercyResult n(goodservers, empty);
    H.push_back(n);

    for (unsigned long i = 0; i < words_per_block; ++i) {
	if (goodservers.size() <= (t+tau)) {
	    std::cerr << "Too few honest servers to recover data!\n";
	    H.clear();
	    return H;
	}

	// std::cerr << i << " " << t+tau << " " << numgood << " " << h << "\n";
	Hprime = EasyRecover(bytes_per_word, t+tau, h, H, answers[i], indices);
	// std::cerr << "easyrecover returned " << res << "\n";
	if (Hprime.empty()) {
	    // We seem to have some Byzantine servers.  Try to identify
	    // them with the expensive algorithm
	    Hprime = HardRecover(bytes_per_word, t+tau, h, H, goodservers,
		    answers[i], indices, p1, p2);
	}
	if (Hprime.empty()) {
	    // Still unable to recover data
	    std::cerr << "Too few honest servers to recover data!\n";
	    H.clear();
	    return H;
	}
	H = Hprime;
	// std::cerr << i << " " << t << " " << numgood << " " << x << "\n";
    }
    return H;
}

// Do all of the above in one shot
vector<PercyResult> PercyClient::fetch_block(unsigned int block_number,
	std::vector<iostream*> &iosvec, ZZ p1, ZZ p2)
{
    int res = send_request(block_number, iosvec);
    if (res < 0) {
	vector<PercyResult> empty;
	return empty;
    }
    unsigned short k = receive_replies(iosvec);
    // Calculate what h to use
    unsigned short h = (unsigned short)(floor(sqrt((t+params.tau())*k)))+1;
    return process_replies(h, p1, p2);
}
