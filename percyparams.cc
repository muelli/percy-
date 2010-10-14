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

#include <string.h>
#include "percyparams.h"
#include "percyio.h"

#define PERCY_VERSION 1

PercyParams::PercyParams()
{
    this->version = PERCY_VERSION;
    this->hybrid_protection = false;
    this->_tau = 0;
    this->_words_per_block = 1;
    this->_num_blocks = 1;
    this->modulus = 257;
    create_ZZ_pContexts();
}

PercyParams::PercyParams(unsigned int words_per_block, unsigned int num_blocks,
	unsigned short tau, ZZ modulus)
{
    this->version = PERCY_VERSION;
    this->hybrid_protection = false;
    this->_tau = tau;
    this->_words_per_block = words_per_block;
    this->_num_blocks = num_blocks;
    this->modulus = modulus;
    create_ZZ_pContexts();
}

ostream& operator<<(ostream& os, const PercyParams &params)
{
    unsigned char minibuf[4];

    // Output the magic header
    os.write("PIRC", 4);

    // Output the version number and flags
    minibuf[0] = params.version;
    minibuf[1] = (params.hybrid_protection ? 1 : 0) | (params._tau ? 2 : 0);
    os.write((char *)minibuf, 2);

    // Output the words_per_block and num_blocks values
    percy_write_le_int(os, params._words_per_block);
    percy_write_le_int(os, params._num_blocks);

    // Output the modulus
    percy_write_ZZ(os, params.modulus);

    // Output g, if appropriate
    if (params.hybrid_protection) percy_write_ZZ(os, rep(params.g));

    return os;
}

istream& operator>>(istream& is, PercyParams &params)
{
    unsigned char minibuf[4];

    // Input the magic header
    is.read((char *)minibuf, 4);
    if (memcmp(minibuf, "PIRC", 4)) {
	std::cerr << "Did not find expected PercyParams header.\n";
	return is;
    }

    // Input the version number and flags
    is.read((char *)minibuf, 2);
    params.version = minibuf[0];
    params.hybrid_protection = ( (minibuf[1] & 1) == 1 );
    params._tau = ( (minibuf[1] & 2) == 2 );
    if (params.version != PERCY_VERSION) {
	std::cerr << "Did not find expected PercyParams version number " <<
	    PERCY_VERSION << ".\n";
	return is;
    }

    // Input the words_per_block and num_blocks values
    percy_read_le_int(is, params._words_per_block);
    percy_read_le_int(is, params._num_blocks);

    // Input the modulus
    percy_read_ZZ(is, params.modulus);
    params.create_ZZ_pContexts();

    // Input g, if appropriate
    if (params.hybrid_protection) {
	ZZ_pContext savectx;
	savectx.save();
	params.modsqctx.restore();
	ZZ gz;
	percy_read_ZZ(is, gz);
	params.g = to_ZZ_p(gz);
	savectx.restore();
    }

    return is;
}

PercyParams::PercyParams(unsigned int words_per_block,
	unsigned int num_blocks, unsigned short tau, ZZ p, ZZ q)
{
    init_hybrid(words_per_block, num_blocks, tau, p, q);
}

void PercyParams::init_hybrid(unsigned int words_per_block,
	unsigned int num_blocks, unsigned short tau, ZZ p, ZZ q)
{
    this->version = PERCY_VERSION;
    this->hybrid_protection = true;
    this->_tau = tau;
    this->_words_per_block = words_per_block;
    this->_num_blocks = num_blocks;
    ZZ modulus;
    modulus = p * q;
    this->modulus = modulus;
    create_ZZ_pContexts();

    // Generate the Paillier public and private parts
    ZZ pm1, qm1;
    pm1 = p - 1;
    qm1 = q - 1;
    this->lambda = pm1 * qm1 / GCD(pm1, qm1);

    ZZ_pContext savectx;
    savectx.save();
    mod_modulussq();
    random(this->g);
    ZZ muinv = rep(power(this->g, this->lambda) - 1) / modulus;
    mod_modulus();
    this->mu = inv(to_ZZ_p(muinv));
    savectx.restore();
}

void PercyParams::create_ZZ_pContexts()
{
    // Create the ZZ_pContexts
    ZZ_pContext modctx(modulus);
    ZZ_pContext modsqctx(modulus * modulus);
    this->modctx = modctx;
    this->modsqctx = modsqctx;
}

PercyParams::PercyParams(unsigned int words_per_block,
	unsigned int num_blocks, unsigned short tau, int modulus_bits)
{
    // Pick the sizes for the primes
    int qsize = modulus_bits / 2;
    int psize = modulus_bits - qsize;

    // Generate random primes of the appropriate size.  We ensure the
    // top two bits are set so that their product is of the right
    // bitlength.
    ZZ pbase, qbase, p, q;
    RandomBits(pbase, psize);
    if (psize >= 1) SetBit(pbase, psize-1);
    if (psize >= 2) SetBit(pbase, psize-2);
    NextPrime(p, pbase);
    RandomBits(qbase, qsize);
    if (qsize >= 1) SetBit(qbase, qsize-1);
    if (qsize >= 2) SetBit(qbase, qsize-2);
    NextPrime(q, qbase);

    init_hybrid(words_per_block, num_blocks, tau, p, q);
}
