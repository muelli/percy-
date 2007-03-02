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
#include "datastore.h"
#include "percyserver.h"

int main(int argc, char **argv)
{
    // Initialize NTL and the random number stream
    ZZ modinit;
    modinit = 257;
    ZZ_p::init(modinit);
    unsigned char randbuf[128];
    ifstream urand("/dev/urandom");
    urand.read((char *)randbuf, sizeof(randbuf));
    urand.close();
    ZZ randzz = ZZFromBytes(randbuf, sizeof(randbuf));
    SetSeed(randzz);

    unsigned int servernum = 0;
    const char *progname = argv[0];

    if (argc > 1 && argv[1][0] == '-') {
	servernum = strtoul(argv[1]+1, NULL, 10);
	--argc;
	++argv;
    }

    if (!(argc == 2 || (argc == 3 && !strcmp(argv[2], "tau")))) {
	std::cerr << "Usage: " << progname <<
	    " [-servernum] database [\"tau\"]\n";
	exit(1);
    }
    const char *dbname = argv[1];
    bool tau_independent = (argc == 3);

    // Read the params for the database from the client
    PercyParams params;
    std::cin >> params;

    // Load the datastore
    FileDataStore datastore(dbname, params, tau_independent);

    // Create the PIR server
    PercyServer server(datastore);

    // With probability $PIRS_FAIL/100, fail completely
    // With probability $PIRS_BYZ/100, be Byzantine
    unsigned long rndval = RandomBnd(100);
    unsigned long failat = 0, byzat = 0;
    const char *failenv = getenv("PIRS_FAIL");
    if (failenv) failat = atoi(failenv);
    const char *byzenv = getenv("PIRS_BYZ");
    if (byzenv) byzat = atoi(byzenv);

    if (rndval < failat) {
	std::cerr << "Server " << servernum << " failing.\n";
	exit(0);
    }
    if (rndval < failat + byzat) {
	std::cerr << "Server " << servernum << " Byzantine.\n";
	server.be_byzantine();
    }

    // Handle the request
    server.handle_request(params, cin, cout);
}
