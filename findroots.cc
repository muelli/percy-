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

#include <ZZ_pX.h>
#include "pstream.h"
#include "findroots.h"

NTL_CLIENT

// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  The global ZZ_pContext should already be set to
// p1 * p2, where p2 is prime, and p1 is either prime or 1.
// This routine may also return some spurious values.
vec_ZZ_pX findroots(const vec_ZZ_pX /* really ZZ_pXY */ &P,
	int degreebound, ZZ p1, ZZ p2)
{
    vec_ZZ_pX roots;

    // Until we implement ZZ_pXY and a root-finding algorithm for it in
    // NTL, we'll have to call out to MuPAD.  But MuPAD can only do this
    // for prime moduli, so we'll have to do it one prime at a time, and
    // combine the results with the CRT if p1 > 1.

    vec_ZZ_pX roots_p1, roots_p2;

    redi::pstream mupadclient("./mupadclient");
    mupadclient << p1 << "\n" << p2 << "\n";
    bool printed = false;
    long ydegree = P.length();
    for (long ydeg = 0; ydeg < ydegree; ++ydeg) {
	long xdegree = deg(P[ydeg]);
	for (long xdeg = 0; xdeg <= xdegree; ++xdeg) {
	    ZZ_p c = coeff(P[ydeg], xdeg);
	    if (c != 0) {
		if (printed) {
		    mupadclient << " + ";
		}
		mupadclient << c << "*x^" << xdeg << "*y^" << ydeg;
		printed = true;
	    }
	}
    }
    mupadclient << "\n";
    mupadclient.flush();

    // Now read back the results.  First for p1
    while (1) {
	ZZ_p ycoeff;
	mupadclient >> ycoeff;
	if (ycoeff == 0) break;
	ZZ_pX xcoeffs;
	mupadclient >> xcoeffs;
	xcoeffs /= (-ycoeff);
	if (degreebound < 0 || deg(xcoeffs) <= degreebound) {
	    append(roots_p1, xcoeffs);
	}
    }

    // Now for p2
    while (1) {
	ZZ_p ycoeff;
	mupadclient >> ycoeff;
	if (ycoeff == 0) break;
	ZZ_pX xcoeffs;
	mupadclient >> xcoeffs;
	xcoeffs /= (-ycoeff);
	if (degreebound < 0 || deg(xcoeffs) <= degreebound) {
	    append(roots_p2, xcoeffs);
	}
    }

    if (p1 > 1) {
	// Calcuate a1 and a2 s.t. a1 = (0, 1) mod (p1,p2) and
	// a2 = (1, 0) mod (p1, p2).
	ZZ_p a1 = to_ZZ_p(p1);
	ZZ_p a2 = to_ZZ_p(p2);
	// cerr << "p1 = " << p1 << "\np2 = " << p2 << "\np1 * p2 = " << p1 * p2 << "\n";
	a1 *= to_ZZ_p(InvMod(AddMod(p1, 0, p2), p2));
	// cerr << "a1 = " << a1 << "\n";
	a2 *= to_ZZ_p(InvMod(AddMod(p2, 0, p1), p1));
	// cerr << "a2 = " << a2 << "\n";
	
	// For each pair, use the CRT to combine them
	unsigned short num_p1 = roots_p1.length();
	unsigned short num_p2 = roots_p2.length();
	for (unsigned short i=0; i<num_p1; ++i) {
	    for (unsigned short j=0; j<num_p2; ++j) {
		ZZ_pX comb = roots_p1[i] * a2 + roots_p2[j] * a1;
		append(roots, comb);
	    }
	}
    } else {
	roots = roots_p2;
    }

    return roots;
}
