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
#include <sstream>
#include "gsdecoder.h"

NTL_CLIENT

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
template<>
string GSDecoder_ZZ_p::append(const string &s, const ZZ_p &wz,
	unsigned int bytes_per_word)
{
    unsigned char w[bytes_per_word];
    BytesFromZZ(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    return r;
}

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
template<>
string GSDecoder_GF2E::append(const string &s, const GF2E &wz,
	unsigned int bytes_per_word)
{
    unsigned char w[bytes_per_word];
    BytesFromGF2X(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    return r;
}

// Project a polynomial created with a ZZ_pContext of k*p down to the
// current ZZ_pContext of p.
static ZZ_pXY project_down(const ZZ_pXY &P)
{
    ZZ_pXY newP;
    stringstream ss (stringstream::in | stringstream::out);
    ss << P;
    ss >> newP;
    return newP;
}

// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  The global ZZ_pContext should already be set to
// p1 * p2, where p2 is prime, and p1 is either prime or 1.
// This routine may also return some spurious values.
template<>
vector<ZZ_pX> GSDecoder_ZZ_p::findroots(
	const ZZ_pXY &P, int degreebound)
{
    // If we're already working mod a prime, just go ahead
    if (p1 == 1) {
	return rr_findroots(P, degreebound);
    }

    // We have to find the roots mod each prime separately, and combine
    // the results with the CRT if p1 > 1.
    ZZ_pBak pbak;
    pbak.save();

    vector<ZZ_pX> roots_p1, roots_p2, roots;

    ZZ_p::init(p1);
    ZZ_pXY P1 = project_down(P);
    roots_p1 = rr_findroots(P1, degreebound);

    ZZ_p::init(p2);
    ZZ_pXY P2 = project_down(P);
    roots_p2 = rr_findroots(P2, degreebound);

    pbak.restore();

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
    unsigned short num_p1 = roots_p1.size();
    unsigned short num_p2 = roots_p2.size();
    for (unsigned short i=0; i<num_p1; ++i) {
	for (unsigned short j=0; j<num_p2; ++j) {
	    ZZ_pX comb = roots_p1[i] * a2 + roots_p2[j] * a1;
	    roots.push_back(comb);
	}
    }

    return roots;
}

#ifdef TEST_FINDPOLYS

#include <sstream>
int main()
{
    ZZ modulus, one;
    modulus = 10007;
    one = 1;
    ZZ_p::init(modulus);
    GSDecoder_ZZ_p decoder(one, modulus);

    stringstream ss(stringstream::in | stringstream::out);

    vector<unsigned short> goodservers;
    goodservers.push_back(1);
    goodservers.push_back(2);
    goodservers.push_back(3);
    goodservers.push_back(4);
    goodservers.push_back(5);
    goodservers.push_back(6);
    vec_ZZ_p indices, shares;
    ss << "[ 1 2 3 4 5 6 ]";
    ss >> indices;
    ss << "[ 4 3 7507 1 2504 7 ]";
    ss >> shares;
    vector<RecoveryPoly<ZZ_pX> > ret = decoder.findpolys(3, 5, goodservers, indices, shares);
    for (vector<RecoveryPoly<ZZ_pX> >::const_iterator iter = ret.begin(); iter != ret.end(); ++iter) {
	cout << "{ " << iter->phi << ", [ ";
	for (vector<unsigned short>::const_iterator gi = iter->G.begin(); gi != iter->G.end(); ++gi) {
	    cout << *gi << " ";
	}
	cout << "] }\n";
    }
}

#endif

#ifdef TEST_RR

// Test the algorithm on the Example 15 from the McElice paper.
int main()
{
    ZZ modulus;
    modulus = 19;
    ZZ p1;
    p1 = 1;
    ZZ_p::init(modulus);
    GSDecoder_ZZ_p decoder(p1, modulus);

    stringstream ss(stringstream::in | stringstream::out);
    ZZ_pXY P;

    ss << "[[4 12 5 11 8 13] [14 14 9 16 8] [14 13 1] [2 11 1] [17]] ";
    ss >> P;

    vector<ZZ_pX> roots = decoder.findroots(P, 1);
    cout << "\nRoots found:\n";
    for (unsigned int i=0; i<roots.size(); ++i) {
	cout << roots[i] << endl;
    }

    cout << "\nExpected output (in some order):\n";
    cout << "[18 14]\n[14 16]\n[8 8]\n";
    return 0;
}
#endif
