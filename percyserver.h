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

#ifndef __PERCYSERVER_H__
#define __PERCYSERVER_H__

#include <iostream>
#include <vec_ZZ_p.h>
#include "datastore.h"
#include "percyparams.h"

NTL_CLIENT

// The class for Percy++ servers
class PercyServer {
public:
    // Initialize a server with the given DataStore.
    PercyServer(DataStore &datastore) : byzantine(false),
					datastore(datastore) {}
    ~PercyServer() {}

    // Tell the server to be Byzantine
    void be_byzantine() { byzantine = true; }

    // Handle a request.
    void handle_request(PercyParams &params, istream &is, ostream &os);

private:
    bool byzantine;
    DataStore &datastore;
    ZZ_p compute_one(bool hybrid_protection, unsigned int num_blocks,
	    const vec_ZZ_p &inputvector, unsigned int c);
};

#endif
