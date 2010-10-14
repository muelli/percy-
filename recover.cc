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

#include <vector>
#include "subset.h"
#include "subset_iter.h"
#include "recover.h"

NTL_CLIENT

vector<PercyResult> EasyRecover_GF28(unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const GF28_Element *values, const GF28_Element *indices)
{
    vector<PercyResult> Hprime;
    vector<PercyResult>::const_iterator Hiter;

    for (Hiter = H.begin(); Hiter != H.end(); ++Hiter) {
	// Pick a random subset I of G, of size t+1
	vector<unsigned short> I;
	vector<unsigned short>::const_iterator Iiter;
	random_subset(Hiter->G, I, t+1);

	// Use Lagrange interpolation to find the unique polynomial phi
	// of degree t which matches the points indexed by I
	GF28_Element I_indices[t+1], I_values[t+1];
	unsigned short i = 0;
	for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
	    I_indices[i] = indices[*Iiter];
	    I_values[i] = values[*Iiter];
	}
	GF28_Element wz;
	wz = interpolate_GF28(I_indices, I_values, t+1, 0);

	// Count the number of points in G that agree, and that
	// disagree, with phi
	unsigned short numagree = 0;
	unsigned short numdisagree = 0;
	vector<unsigned short> vecagree;
	vector<unsigned short>::const_iterator Giter;
	for (Giter = Hiter->G.begin(); Giter != Hiter->G.end(); ++Giter) {
	    if (interpolate_GF28(I_indices, I_values, t+1, indices[*Giter])
			== values[*Giter]) {
		++numagree;
		vecagree.push_back(*Giter);
	    } else {
		++numdisagree;
	    }
	}

	// If at least h agreed, and less than h-t disagreed, then phi
	// can be the *only* polynomial of degree t matching at least
	// h points.
	if (numagree >= h && numdisagree < h-t) {
	    string news(Hiter->sigma);
	    news.append((char*)&wz, 1);
	    PercyResult n(vecagree, news);
	    Hprime.push_back(n);
	} else {
	    // This either isn't the right polynomial, or there may be
	    // more than one.  Abort.
	    Hprime.clear();
	    return Hprime;
	}
    }
    return Hprime;
}

vector<PercyResult> EasyRecover(unsigned int bytes_per_word, unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const vec_ZZ_p &values, const vec_ZZ_p &indices)
{
    vector<PercyResult> Hprime;
    vector<PercyResult>::const_iterator Hiter;

    for (Hiter = H.begin(); Hiter != H.end(); ++Hiter) {
	// Pick a random subset I of G, of size t+1
	vector<unsigned short> I;
	random_subset(Hiter->G, I, t+1);

	unsigned short numagree, numdisagree;
	vector<unsigned short> vecagree;

	ZZ_pX phi;
	GSDecoder_ZZ_p::test_interpolate(t, values, indices, I, Hiter->G,
		numagree, numdisagree, vecagree, phi);

	// If at least h agreed, and less than h-t disagreed, then phi
	// can be the *only* polynomial of degree t matching at least
	// h points.
	if (numagree >= h && numdisagree < h-t) {
	    // Find the secret determined by phi
	    ZZ_p wz;
	    eval(wz, phi, ZZ_p::zero());

	    PercyResult n(vecagree, GSDecoder_ZZ_p::append(Hiter->sigma,
			wz, bytes_per_word));
	    Hprime.push_back(n);
	} else {
	    // This either isn't the right polynomial, or there may be
	    // more than one.  Abort, and we'll use HardRecover.
	    Hprime.clear();
	    return Hprime;
	}
    }
    return Hprime;
}

