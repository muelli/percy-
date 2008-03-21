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

#ifndef __PERCYPARAMS_H__
#define __PERCYPARAMS_H__

#include <iostream>
#include <NTL/ZZ_p.h>

NTL_CLIENT

// The base class for specifying (usually on the client side) the
// parameters of the database.  This is the class to use except when you
// want hybrid protection.
class PercyParams {
    friend ostream& operator<<(ostream& os, const PercyParams &params);
    friend istream& operator>>(istream& is, PercyParams &params);

public:
    PercyParams();

    // Use the given modulus directly; only encryption will be possible
    PercyParams(unsigned int words_per_block, unsigned int num_blocks,
	    unsigned short tau, ZZ modulus);

    // Generate a new public/private key pair of the given keysize
    PercyParams(unsigned int words_per_block, unsigned int num_blocks,
	    unsigned short tau, int modulus_bits);

    // Use the given factors to generate a public/private key pair
    PercyParams(unsigned int words_per_block, unsigned int num_blocks,
	    unsigned short tau, ZZ p, ZZ q);

    unsigned int modulus_bytes() const {
	return NumBytes(modulus);
    }
    unsigned int modulussq_bytes() const {
	return NumBytes(modulus*modulus);
    }
    unsigned int bytes_per_word() const {
	unsigned int r = NumBytes(modulus);
	if (r > 1) return r - 1;
	return 1;
    }
    unsigned int words_per_block() const {
	return _words_per_block;
    }
    unsigned int num_blocks() const {
	return _num_blocks;
    }
    unsigned int bytes_per_block() const {
	return bytes_per_word() * _words_per_block;
    }
    unsigned short tau() const {
	return _tau;
    }
    void mod_modulus() const {
	modctx.restore();
    }
    void mod_modulussq() const {
	modsqctx.restore();
    }
    bool hybrid() const {
	return hybrid_protection;
    }
    bool modulus_match(ZZ testmod) const {
	return modulus == testmod;
    }
    bool is_gf28() const {
	return modulus == 256;
    }
    // Encrypt the given plaintext.  The current ZZ_p context must be
    // modsqctx.
    ZZ_p encrypt(ZZ plaintext) const {
	ZZ_p r;
	random(r);
	return power(g, plaintext) * power(r, modulus);
    }
    // Decrypt the given ciphertext.  This routine will change the
    // current ZZ_p context to modctx.
    ZZ_p decrypt(ZZ_p ciphertext) const {
	modsqctx.restore();
	ZZ Lval = rep(power(ciphertext, lambda) - 1) / modulus;
	modctx.restore();
	ZZ_p ret = to_ZZ_p(Lval) * mu;
	return ret;
    }

protected:
    void create_ZZ_pContexts();
    void init_hybrid(unsigned int words_per_block, unsigned int num_blocks,
	    unsigned short tau, ZZ p, ZZ q);
    unsigned char version;
    bool hybrid_protection;
    unsigned short _tau;
    unsigned int _words_per_block, _num_blocks;
    ZZ_pContext modctx, modsqctx;
    // Paillier public key
    ZZ modulus;
    ZZ_p g;   // mod modulus^2
    // Paillier private key
    ZZ lambda;
    ZZ_p mu;  // mod modulus
};

#endif
