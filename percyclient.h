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

#ifndef __PERCYCLIENT_H__
#define __PERCYCLIENT_H__

#include <vector>
#include <iostream>
#include <vec_vec_ZZ_p.h>
#include "percyresult.h"
#include "percyparams.h"

NTL_CLIENT

// The class for Percy++ clients
class PercyClient {
public:
    // Construct a Percy client
    PercyClient(PercyParams &params, unsigned short num_servers,
	    unsigned short t) : params(params), num_servers(num_servers),
	    t(t) {}

    ~PercyClient() {}

    // Send a request for the given block number (0-based) to the
    // servers connected with the ostreams in the given vector.
    int send_request(unsigned int block_number,
	    std::vector<iostream*> &iosvec);

    // Receive the servers' replies.  Return k, the number of servers
    // that returned anything.
    unsigned short receive_replies(std::vector<iostream*> &iosvec);

    // Process the server's replies, and return a new string
    // (of size params.bytes_per_block()) with the result.
    vector<PercyResult> process_replies(unsigned short h, ZZ p1, ZZ p2);

    // Do all of the above in one shot
    vector<PercyResult> fetch_block(unsigned int block_number,
	    std::vector<iostream*> &iosvec, ZZ p1, ZZ p2);

private:
    PercyParams params;
    unsigned short num_servers, t;
    vec_ZZ_p indices;
    vector<unsigned short> goodservers;
    vec_vec_ZZ_p answers;
};

#endif
