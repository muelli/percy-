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

// This file is an implementation of the Roth-Ruckenstein algorithm
// for finding roots of bivariate polynomials, as explained in section
// IX of:  R. J. McEliece. The Guruswami-Sudan Decoding Algorithm for
// Reed-Solomon Codes. IPN Progress Report 42-153, May 15, 2003.
// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html

#include <sstream>
#include <vector>
#include "ZZ_pXY.h"

NTL_CLIENT

// Depth-first search on the tree of coefficients
static void dfs(vec_ZZ_pX &res, 
		int u, 
		vector<int> &pi, 
		vector<int> &Deg,
		vector<ZZ_p> &Coeff, 
		vector<ZZ_pXY> &Q, 
		int &t, 
		int degreebound)
{
#ifdef TEST_RR
    cout << "\nVertex " << u << ": pi[" << u << "] = " << pi[u] <<
	", deg[" << u << "] = " << Deg[u] << ", Coeff[" << u <<
	"] = " << Coeff[u] << "\n";
    cout << "Q[" << u << "] = " << Q[u] << "\n";
#endif
    // if Q_u(x,0) == 0 then output y-root
    if ( IsZero(Q[u]) || IsZero(Q[u].rep[0])) { 
	// Output f_u[x]

	ZZ_pX fux;
	while (Deg[u] >= 0) {
	    SetCoeff(fux, Deg[u], Coeff[u]);
	    u = pi[u];
	}
#ifdef TEST_RR
	cout << "Outputting " << fux << "\n";
#endif
	append(res, fux);
    } else if (degreebound < 0 || Deg[u] < degreebound) {
	// Construct Q_u(0,y)
	ZZ_pX Qu0y;
	int degqu = deg(Q[u]);
	for (int d = 0; d <= degqu; ++d) {
	    SetCoeff(Qu0y, d, coeff(Q[u].rep[d], 0));
	}
      
#ifdef TEST_RR
	cout << "Q[" << u << "](0,y) = " << Qu0y << "\n";
#endif
	// Find its roots
	vec_ZZ_p rootlist = findroots(Qu0y);
#ifdef TEST_RR
	cout << "Rootlist = " << rootlist << "\n";
#endif
	int numroots = rootlist.length();
	for (int r=0; r<numroots; ++r) {
	    // For each root a, recurse on <<Q[u](x, x*y + a)>>
	    // where <<Q>> is Q/x^k for the maximal k such that the
	    // division goes evenly.
	    int v = t;
	    ++t;
	    pi.push_back(u);
	    Deg.push_back(Deg[u]+1);
	    Coeff.push_back(rootlist[r]);
	    // mapit(Q(x,y), a) computes <<Q(x, x*y + a)>>
	    Q.push_back(mapit(Q[u], rootlist[r]));
	    dfs(res, v, pi, Deg, Coeff, Q, t, degreebound);
	}
    }
}

// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  The global ZZ_pContext should already be set to
// the appropriate prime.
vec_ZZ_pX rr_findroots(const ZZ_pXY &P, int degreebound)
{
    vector<int> pi;
    vector<int> Deg;
    vector<ZZ_p> Coeff;
    vector<ZZ_pXY> Q;
    vec_ZZ_pX res;
    int t = 1;

    ZZ_pXY P0 = backShiftX(P, minX(P));
    pi.push_back(-1);
    Deg.push_back(-1);
    Coeff.push_back(ZZ_p::zero());
    Q.push_back(P0);
    int u = 0; 

    if (degreebound < 0) {
	degreebound = degX(P0);
    }
    dfs(res, u, pi, Deg, Coeff, Q, t, degreebound);

    return res;
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
vec_ZZ_pX findroots(const ZZ_pXY &P, int degreebound, ZZ p1, ZZ p2)
{
    // If we're already working mod a prime, just go ahead
    if (p1 == 1) {
	return rr_findroots(P, degreebound);
    }

    // We have to find the roots mod each prime separately, and combine
    // the results with the CRT if p1 > 1.
    ZZ_pBak pbak;
    pbak.save();

    vec_ZZ_pX roots_p1, roots_p2, roots;

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
    unsigned short num_p1 = roots_p1.length();
    unsigned short num_p2 = roots_p2.length();
    for (unsigned short i=0; i<num_p1; ++i) {
	for (unsigned short j=0; j<num_p2; ++j) {
	    ZZ_pX comb = roots_p1[i] * a2 + roots_p2[j] * a1;
	    append(roots, comb);
	}
    }

    return roots;
}

#ifdef TEST_RR

// Test the algorithm on the Example 15 from the McElice paper.
int main()
{
    ZZ modulus;
    modulus = 19;
    ZZ p1;
    p1 = 1;
    ZZ_p::init(modulus);

    stringstream ss(stringstream::in | stringstream::out);
    ZZ_pXY P;

    ss << "[[4 12 5 11 8 13] [14 14 9 16 8] [14 13 1] [2 11 1] [17]] ";
    ss >> P;

    vec_ZZ_pX roots = rr_findroots(P, 1);
    cout << "\nRoots found = " << roots << "\n";

    return 0;
}
#endif
