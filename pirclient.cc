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

#include "pstream.h"
#include <signal.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ZZ_p.h>
#include "percyclient.h"
#include "config.h"

NTL_CLIENT

int main(int argc, char **argv)
{
    int c;
    for (c = 1; c < argc; c++) {
	if (strcmp("--version", argv[c]) == 0) {
	    std::cerr << "Percy++ pirclient version " << VERSION << std::endl;
	    std::cerr << AUTHOR << std::endl;
	    exit(0);
	}
    }

    // Ignore SIGPIPE
    signal(SIGPIPE, SIG_IGN);

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

    if (argc < 9) {
	std::cerr << "Usage: " << argv[0] << " n b w l k t i tau [hybrid|gf28]\n";
	std::cerr << "   n = length of database (bytes)\n";
	std::cerr << "   b = block size (bytes)\n";
	std::cerr << "   w = word size (bits)\n";
	std::cerr << "   l = number of servers\n";
	std::cerr << "   k = number of servers that need to respond\n";
	std::cerr << "   t = number of servers that can collude\n";
	std::cerr << "   i = index of block to fetch (0-based)\n";
	std::cerr << "   tau = tau-independence value of database shares\n";
	exit(1);
    }

    unsigned long long n = strtoull(argv[1], NULL, 10) * 8;
    unsigned int b = strtoul(argv[2], NULL, 10) * 8;
    unsigned int w = strtoul(argv[3], NULL, 10);
    unsigned int num_servers = strtoul(argv[4], NULL, 10);
    unsigned int k = strtoul(argv[5], NULL, 10);
    unsigned int t = strtoul(argv[6], NULL, 10);
    vector<unsigned int> indices;
    istringstream iss(argv[7]);
    while(1) {
	unsigned int i;
	iss >> i;
	indices.push_back(i);
	if (iss.eof()) break;
    }
    unsigned int tau = strtoul(argv[8], NULL, 10);
    int do_hybrid = (argc > 9 && !strcmp(argv[9], "hybrid"));
    int do_gf28 = (argc > 9 && !strcmp(argv[9], "gf28"));

    if (do_gf28 && w != 8) {
	std::cerr << "Error: w must be 8 for gf28.\n";
	exit(1);
    }

    if (n % b != 0 || b % w != 0) {
	std::cerr << "Error: b must divide n and w must divide b.\n";
	exit(1);
    }

    if (w % 8 != 0) {
	std::cerr << "Error: 8 must divide w.\n";
	exit(1);
    }

    unsigned int num_blocks = n/b;
    unsigned int words_per_block = b/w;
    std::cerr << "Number of blocks: " << num_blocks << "\n";
    std::cerr << "Words per block:  " << words_per_block << "\n";

    if (num_blocks != words_per_block) {
	std::cerr << "Warning: non-optimal choice of blocksize detected.\n";
    }

    stringstream ss (stringstream::in | stringstream::out);
    ZZ p1, p2;
    const char *p1s, *p2s;

    if (w == 2048) {
	p1s = "208647130951457402363969335056365957472826150618980217460328400485971950387185944410889077723063406198415802830757517777351462262669194793047360775411639408116452523756687066355086195124187048682420529316060567502352699557841412039275095485224490337148164650010000499984813523719988826268799665657866626493329 ";
	p2s = "245210205383950153265232956846271987008710436579074459102383753214859717124121302267932009072054546430711727811323033561244148876933172687995163379778095734152594201215411509169035373484564340604271927100344464582888777887746436564737355045100633587336239754449508771770564607896955672999950235015535154415867 ";
    } else if (w == 1536) {
	p1s = "1762848592595080314705600925431624874456855439794595868418019480189213868063348394981842423875338178991362893712297567682392276281463789141688306484765105096429658863055172316227409205756175078509101834587188923103831602929062176351 ";
	p2s = "2306072568237159640249655953989533876736033293267891813492402870069702343561490811306173449455816459207943593196801405361355605814646339972518285709494570145269396000374210678514250118174550977925517522811232946347459478425104006037 ";
    } else if (w == 1024) {
	p1s = "14710132128541592475387440366744304824352604767753216777226640368050037133836174845369895150342922969891066267019166301546403100960464521216972792406229873 ";
	p2s = "23338263930359653850870152235447790059566299230626918909126959284529524161146399225204807633841208114717867386116272471253667601589249734409576687328158817 ";
    } else if (do_gf28) {
	p1s = "1 ";
	p2s = "256 ";
    } else if (w == 8 && !do_hybrid) {
	p1s = "1 ";
	p2s = "257 ";
    } else if (w == 16 && !do_hybrid) {
	p1s = "1 ";
	p2s = "65537 ";
    } else if (w == 32 && !do_hybrid) {
	p1s = "1 ";
	p2s = "4294967311 ";
    } else if (w == 96 && !do_hybrid) {
	p1s = "1 ";
	p2s = "79228162514264337593543950397 ";
    } else if (w == 128 && !do_hybrid) {
	p1s = "1 ";
	p2s = "340282366920938463463374607431768211507 ";
    } else if (w == 192 && !do_hybrid) {
	p1s = "1 ";
	p2s = "6277101735386680763835789423207666416102355444464034513029 ";
    } else if (w == 256 && !do_hybrid) {
	p1s = "1 ";
	p2s = "115792089237316195423570985008687907853269984665640564039457584007913129640233 ";
    } else {
	std::cerr << "No modulus available for w = " << w << "\n";
	exit(1);
    }

    ss << p1s;
    ss >> p1;

    ss << p2s;
    ss >> p2;

    ZZ modulus;
    modulus = p1 * p2;

    unsigned int ybytes = NumBytes(modulus);

    if (ybytes <= w/8) {
	std::cerr << "Error: w must be at most " << (ybytes-1)*8 <<
	    " (was " << w << ").\n";
	exit(1);
    }

    if (k < 1) {
	std::cerr << "Error: k must be at least 1.\n";
	exit(1);
    }

    if (k > num_servers) {
	std::cerr << "Error: k must be at most l.\n";
	exit(1);
    }

    if (t+tau >= k) {
	std::cerr << "Error: t+tau must be less than k.\n";
	exit(1);
    }

    for (size_t q = 0; q < indices.size(); ++q) {
	if (indices[q] >= num_blocks) {
	    std::cerr << "Error: i must be less than " << num_blocks << ".\n";
	    exit(1);
	}
    }

    // Set up the iostreams to the servers.  In reality, you'd probably
    // use an SSL connection to a remote server here.  As a simple demo,
    // we'll just spawn a local process.
    vector<iostream*> serverstreams;
    for (unsigned int j = 0; j < num_servers; ++j) {
	cerr << "Spawning server " << j+1 << " of " << num_servers << "...\n";
	char command[200];
	char dbname[30];
	strcpy(dbname, "database");
	if (tau > 0) {
	    sprintf(dbname, "database.%u tau", j+1);
	}
	sprintf(command, "time ./pirserver -%u %s", j+1, dbname);
	serverstreams.push_back(new redi::pstream(command));
    }

    // Do the PIR query
    PercyParams *params;
    if (do_hybrid) {
	params = new PercyParams(words_per_block, num_blocks, tau, p1, p2);
    } else {
	params = new PercyParams(words_per_block, num_blocks, tau, modulus);
    }
    PercyClient client(*params, num_servers, t);
    vector< vector<PercyResult> > results = client.fetch_blocks(indices,
	    serverstreams, p1, p2);
    int ret = 0;
    int num_res = results.size();
    for (int r=0; r < num_res; ++r) {
	if (results[r].empty()) {
	    std::cerr << "PIR query failed.\n";
	    ret = 1;
	} else if (results[r].size() > 1) {
	    std::cerr << results[r].size() << " possible blocks returned.\n";
	}
	// Output the retrieved block(s)
	vector<PercyResult>::const_iterator resiter;
	for (resiter = results[r].begin(); resiter != results[r].end();
		++resiter) {
	    std::cout << resiter->sigma;
	}
    }

    // Clean up
    for (unsigned int j = 0; j < num_servers; ++j) {
	delete serverstreams[j];
	serverstreams[j] = NULL;
    }
    serverstreams.clear();
    delete params;

    return ret;
}
