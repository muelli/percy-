//  Percy++ Copyright 2007 Ian Goldberg <iang@cs.uwaterloo.ca>
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

#ifndef __GSDECODER_IMPL_H__
#define __GSDECODER_IMPL_H__

#include <vector>
#include <map>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <vec_ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include "subset.h"
#include "subset_iter.h"
#include "percyresult.h"
#include "FXY.h"

// An implementation of the Guruswami-Sudan decoder.  This decoder will
// produce all possible codewords of degree at most t that match the
// given codeword with k non-erasures in more than sqrt(k*t) places.
//
// It runs in polynomial time, but in the worst case, it's O(k^12).
// Luckily, the worst case is when t=k-2, which is easily solved by
// brute force.  (You will need to do this brute force outside this
// class.)
//
// This implementation is based on the K\"otter and Roth-Ruckenstein
// algorithms, as described in:  R. J. McEliece. The Guruswami-Sudan
// Decoding Algorithm for Reed-Solomon Codes. IPN Progress Report
// 42-153, May 15, 2003.
// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html

// This file contains the (templated) implementations of the class
// methods.

// Specializations for the derived types.  Only vec_FX and
// vec_pair_FX_long are in here, because they're only used once.  The
// rest of the derived types (vec_F, FX, FXY) are used many times, so
// it's more consice to make them explicit instead of having to write
// DT::vec_F every time.
template<>
struct GSDecoder_ZZ_p::DT {
    typedef vec_ZZ_pX vec_FX;
    typedef vec_pair_ZZ_pX_long vec_pair_FX_long;
};

template<>
struct GSDecoder_GF2E::DT {
    typedef vec_GF2EX vec_FX;
    typedef vec_pair_GF2EX_long vec_pair_FX_long;
};

// Set phi to the degree-t polynomial which interpolates the
// (indices,values) pairs indexed by I.  Check which
// (indices,values) pairs indexed by G agree and disagree with
// phi, and set numagree, numdisagree, and vecagree accordingly.
template<class F, class vec_F, class FX, class FXY>
void GSDecoder<F,vec_F,FX,FXY>::test_interpolate(unsigned short t,
	const vec_F &values, const vec_F &indices,
	const vector<unsigned short> &I, const vector<unsigned short> &G,
	unsigned short &numagree, unsigned short &numdisagree,
	vector<unsigned short> &vecagree, FX &phi)
{
    numagree = numdisagree = 0;

    // Use Lagrange interpolation to find the unique polynomial phi
    // of degree t which matches the points indexed by I
    vec_F I_indices, I_values;
    I_indices.SetLength(t+1);
    I_values.SetLength(t+1);
    vector<unsigned short>::const_iterator Iiter;
    unsigned short i = 0;
    for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
	I_indices[i] = indices[*Iiter];
	I_values[i] = values[*Iiter];
    }
    interpolate(phi, I_indices, I_values);

    // Count the number of points in G that agree, and that
    // disagree, with phi
    vector<unsigned short>::const_iterator Giter;
    for (Giter = G.begin(); Giter != G.end(); ++Giter) {
	F phival;
	eval(phival, phi, indices[*Giter]);
	if (phival == values[*Giter]) {
	    ++numagree;
	    vecagree.push_back(*Giter);
	} else {
	    ++numdisagree;
	}
    }
}

extern unsigned long long hasseop;

// Evaluate the (r,s)th Hasse mixed partial derivative of g at the point
// (alpha, beta), which is:
//   \sum_{i,j} C(i,r) C(j,s) a_{i,j} alpha^{i-r} beta^{j-s}
// = \sum_j C(j,s) beta^{j-s} \sum_i C(i,r) a_{i,j} alpha^{i-r}
// where g(x,y) = \sum_{i,j} a_{i,j} x^i y^j
// See page 14 of R. J. McEliece. The Guruswami-Sudan Decoding
// Algorithm for Reed-Solomon Codes. IPN Progress Report 42-153, May 15,
// 2003.  http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
template<class F, class vec_F, class FX, class FXY>
F GSDecoder<F,vec_F,FX,FXY>::evalhasse(const FXY &g,
	unsigned int r, unsigned int s,
	F alpha, F beta, Ccache_t &Ccache)
{
    F res;
    res = 0;
    int ydeg = deg(g);
    for (int j = ydeg; j >= (int)s; --j) {
	const FX &gj = coeff(g,j);
	int xdeg = deg(gj);
	// Use Horner's method to evaluate the inner sum (which is the
	// coefficient of beta^{j-s})
	F resx;
	resx = 0;
	for (int i = xdeg; i >= (int) r; --i) {
	    resx *= alpha;
	    resx += C(Ccache, i,r) * coeff(gj, i);
	    ++hasseop;
	}
	// Use Horner's method to accumulate the results into the outer
	// sum
	res *= beta;
	res += C(Ccache, j,s) * resx;
	++hasseop;
    }

    return res;
}

extern unsigned long long kotter_usec;

// Construct a bivariate polynomial P(x,y) such that, for any polynomial
// f(x) of degree at most v that agrees with at least t of the given
// points, (y-f(x)) is a factor of P(x,y).  This version is K"otter's
// Interpolation Algorithm, as described in Section VII of R. J.
// McEliece. The Guruswami-Sudan Decoding Algorithm for Reed-Solomon
// Codes. IPN Progress Report 42-153, May 15, 2003.
// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
template<class F, class vec_F, class FX, class FXY>
FXY GSDecoder<F,vec_F,FX,FXY>::interpolate_kotter(
	unsigned int v, unsigned int t,
	const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
    struct timeval st, et;
    gettimeofday(&st, NULL);

    unsigned int n = goodservers.size();

    // Compute the m and L parameters
    unsigned int m = 1 + (unsigned int)(floor( v*n / (t*t-v*n)));
    unsigned int L = (m*t - 1)/v;

    if (getenv("PIRC_L")) L = atoi(getenv("PIRC_L"));
    if (getenv("PIRC_m")) m = atoi(getenv("PIRC_m"));

    std::cerr << "Constructing (1," << v << ")-degree " << L*v << " polynomial...\n";
    std::cerr << "Estimated work: " << n * m * (m+1) / 2 * (L+1) << "\n";
    std::cerr << "L = " << L << "\n";
    std::cerr << "m = " << m << "\n";
    std::cerr << "Min matches: " << t << "\n";
    std::cerr << "Max degree: " << v << "\n";
#if 0
    double Km = v * n * (m+1);
    Km /= (double) m;
    Km = floor(sqrt(Km));
    std::cerr << "Km ~= " << Km << "\n";
    unsigned int C = n * m * (m+1) / 2;
    for (int K=0;;++K) {
	cerr << (K*(K+v)+(K%v)*(v-(K%v)))/(2*v) << " " << C << "\n";
	if ( ((K*(K+v)+(K%v)*(v-(K%v)))/(2*v)) > C ) {
	    std::cerr << "Km: " << (K-1)/m + 1 << "\n";
	    break;
	}
    }
#endif

    // Initialize the g vector
    typedef pair<FXY, unsigned int> polydeg;
    polydeg g[L+1];
    for (unsigned int j = 0; j <= L; ++j) {
	SetCoeff(g[j].first, j);
	g[j].second = j * v;
    }

    Ccache_t Ccache;

    for (unsigned int i = 0; i < n; ++i) {
	F alpha = indices[goodservers[i]];
	F beta = shares[goodservers[i]];
	for (unsigned int r = 0; r < m; ++r) {
	    for (unsigned int s = 0; s < m - r; ++s) {
		int seennonzero = 0;
		unsigned int seendeg = 0, jstar = 0;
		F Delta[L+1];
		for (unsigned int j = 0; j <= L; ++j) {
		    Delta[j] = evalhasse(g[j].first, r, s, alpha, beta,
			    Ccache);
		    // cerr << i << " " << r << " " << s << " " << j;
		    if (Delta[j] != 0) {
			// cerr << " nonzero";
			seennonzero = 1;
			seendeg = g[j].second;
			jstar = j;
		    }
		    // cerr << "\n";
		}
		if (seennonzero) {
		    for (unsigned int j = 0; j <= L; ++j) {
			if (Delta[j] != 0 && g[j].second <= seendeg) {
			    seendeg = g[j].second;
			    jstar = j;
			}
		    }
		    F Deltajstar = Delta[jstar];
		    FXY f = g[jstar].first;
		    // cerr << "Deltajstar = " << Deltajstar << "\n";
		    // cerr << "f = " << f << "\n";
		    for (unsigned int j = 0; j <= L; ++j) {
			if (Delta[j] != 0) {
			    if (j != jstar) {
				// cerr << "g["<<j<<"] = " << Deltajstar << " * " << g[j].first << " - " << Delta[j] << " * " << f << " = ";

				g[j].first = Deltajstar * g[j].first -
				    Delta[j] * f;
				// cerr << g[j].first << "\n";
			    } else {
				FX xminusalpha;
				SetCoeff(xminusalpha, 1, 1);
				SetCoeff(xminusalpha, 0, -alpha);
				// cerr << "g["<<j<<"] = " << Deltajstar << " * " << xminusalpha << " * " << f << " = ";
				g[j].first = Deltajstar * xminusalpha * f;
				// cerr << g[j].first << "\n";
				g[j].second += 1;
			    }
			}
			// cerr << "Now -> " << evalhasse(g[j].first, r, s, alpha, beta, Ccache) << "\n";
		    }
		}
	    }
	}
    }
    // Return the poly of least weighted degree from g
    unsigned int minweight = g[0].second;
    unsigned int minindex = 0;
    for (unsigned int i=1; i<=L; ++i) {
	if (g[i].second <= minweight) {
	    minweight = g[i].second;
	    minindex = i;
	}
    }

    gettimeofday(&et, NULL);
    kotter_usec = ((unsigned long long)(et.tv_sec - st.tv_sec))*1000000
	+ (et.tv_usec - st.tv_usec);

    return g[minindex].first;
}

// Find all polynomials of degree at most k that agree with the given
// (index,share) pairs at at least t of the n points.  The notation is
// from Venkatesan Guruswami and Madhu Sudan, "Improved Decoding of
// Reed-Solomon and Algebraic-Geometry Codes".
template<class F, class vec_F, class FX, class FXY>
vector< RecoveryPoly<FX> > GSDecoder<F,vec_F,FX,FXY>::findpolys(unsigned int k,
	unsigned int t, const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
    char *bruteenv = getenv("PIRC_BRUTE");
    if (bruteenv && atoi(bruteenv)) {
	return findpolys_brute(k, t, goodservers, indices, shares);
    }

    FXY P;
    //char *naiveenv = getenv("PIRC_NAIVE");
    //if (naiveenv && atoi(naiveenv)) {
	//P = interpolate_naive(k, t, goodservers, indices, shares);
    //} else {
	P = interpolate_kotter(k, t, goodservers, indices, shares);
    //}

#if 0
    // cerr << "factor(poly(0";
    for(int j=0;j<=deg(P);++j) {
	FX x = coeff(P,j);
	for(int i=0; i<=deg(x); ++i) {
	    // cerr << " + " << coeff(x,i) << "*x^" << i << "*y^" << j;
	}
    }
#endif
    // cerr << ", [x,y], IntMod(" << p1*p2 << ")));\n\n";
    std::cerr << "Finding roots of resulting polynomial...\n";

    // It turns out that any polynomial phi(x) that we should be
    // returning (since phi is of degree at most k and agrees with the
    // input data on at least t points) is such that (y - phi(x)) is a
    // factor of P(x,y).  So we first find all roots for y of P(x,y)
    // which are polynomials of degree at most k.
    vector<FX> roots = findroots(P, k);
    // cerr << "roots = " << roots << "\n";

    // For each of these roots, check how many input points it agrees
    // with.  If it's at least t, add it to the list of polys to return.
    vector< RecoveryPoly<FX> > polys;
    unsigned int numroots = roots.size();
    for (unsigned int i=0; i<numroots; ++i) {
	if (deg(roots[i]) > (long)k) continue;
	vector<unsigned short>::const_iterator gooditer;
	vector<unsigned short> vecagree;
	unsigned short numagree = 0;
	for (gooditer = goodservers.begin(); gooditer != goodservers.end();
		++gooditer) {
	    F phival;
	    eval(phival, roots[i], indices[*gooditer]);
	    if (phival == shares[*gooditer]) {
		++numagree;
		vecagree.push_back(*gooditer);
	    }
	}
	if (numagree >= t) {
	    RecoveryPoly<FX> n(vecagree, roots[i]);
	    polys.push_back(n);
	}
    }

    return polys;
}

// Find all polynomials of degree at most k that agree with the given
// (index,share) pairs at at least t of the n points.  This version
// simply uses brute force and interpolates all subsets of size k+1 of
// the n points.  Note that in general, this version does *not* run in
// polynomial time, but for some cases (with n-k small, for example) it
// is faster than Guruswami-Sudan.
//
// Runtime is about C(k,t+1) *
// [2.55712398139188+1.17564033117597*k+4.20999565858028*t+1.21270558067138*t^2]
// microseconds, according to observations and regression analysis.
template<class F, class vec_F, class FX, class FXY>
vector<RecoveryPoly<FX> > GSDecoder<F,vec_F,FX,FXY>::findpolys_brute(
	unsigned int k, unsigned int t,
	const vector<unsigned short> &goodservers,
	const vec_F &indices, const vec_F &shares)
{
    vector<RecoveryPoly<FX> > polys;

unsigned int numcombs = 0;
    // Iterate over all subsets of goodservers of size k+1
    subset_iterator iter(goodservers, k+1);

    while(!iter.atend()) {
++numcombs;

	unsigned short numagree, numdisagree;
	vector<unsigned short> vecagree;
	FX phi;

	// We should probably avoid calling this if polys[i].G is a
	// superset of *iter for some i, but "is a superset" is a
	// non-trivial test with the current data structure, so we'll
	// just run it anyway for now, and check for uniqueness of phi
	// on the way out.
	test_interpolate(k, shares, indices, *iter, goodservers,
		numagree, numdisagree, vecagree, phi);

	if (numagree >= t) {
	    // As above: check to see if we've seen this phi before
	    typename vector<RecoveryPoly<FX> >::iterator polysiter;
	    bool matched = false;
	    for (polysiter = polys.begin(); polysiter != polys.end();
		    ++polysiter) {
		if (polysiter->phi == phi) {
		    matched = true;
		    break;
		}
	    }
	    if (!matched) {
		RecoveryPoly<FX> n(vecagree, phi);
		polys.push_back(n);
	    }
	}

	++iter;
    }

    return polys;
}

template<class F, class vec_F, class FX, class FXY>
vector<PercyResult> GSDecoder<F,vec_F,FX,FXY>::HardRecover(
	unsigned int bytes_per_word, unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const vector<unsigned short> &goodservers,
	const vec_F &values, const vec_F &indices)
{
    vector<PercyResult> Hprime;

    // Find all polynomials of degree at most t that match at least h of
    // the (index, value) pairs, mod p1*p2  (p1 might be 1).
    vector< RecoveryPoly<FX> > polys = findpolys(t, h, goodservers,
	    indices, values);
    typename vector< RecoveryPoly<FX> >::const_iterator Piter;

    for (Piter = polys.begin(); Piter != polys.end(); ++Piter) {
	// Find the secret determined by this poly
	F wz;
	eval(wz, Piter->phi, F::zero());

	// Find any elements (G,sigma) of H with |Piter->G \cap G| \ge h
	vector<PercyResult>::const_iterator Hiter;
	for (Hiter = H.begin(); Hiter != H.end(); ++Hiter) {
	    vector<unsigned short> intersection =
		intersect(Piter->G, Hiter->G);
	    if (intersection.size() >= h) {
		PercyResult n(intersection,
			append(Hiter->sigma, wz, bytes_per_word));
		Hprime.push_back(n);
	    }
	}

    }

    return Hprime;
}

// Depth-first search on the tree of coefficients; used by the
// Roth-Ruckenstein algorithm.
template<class F, class vec_F, class FX, class FXY>
void GSDecoder<F,vec_F,FX,FXY>::dfs(vector<FX> &res,
		int u, 
		vector<int> &pi, 
		vector<int> &Deg,
		vector<F> &Coeff, 
		vector<FXY> &Q, 
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

	FX fux;
	while (Deg[u] >= 0) {
	    SetCoeff(fux, Deg[u], Coeff[u]);
	    u = pi[u];
	}
#ifdef TEST_RR
	cout << "Outputting " << fux << "\n";
#endif
	res.push_back(fux);
    } else if (degreebound < 0 || Deg[u] < degreebound) {
	// Construct Q_u(0,y)
	FX Qu0y;
	int degqu = deg(Q[u]);
	for (int d = 0; d <= degqu; ++d) {
	    SetCoeff(Qu0y, d, coeff(Q[u].rep[d], 0));
	}
      
#ifdef TEST_RR
	cout << "Q[" << u << "](0,y) = " << Qu0y << "\n";
#endif
	// Find its roots
	vec_F rootlist =
	    findroots_FX<F,vec_F,FX,
		typename DT::vec_FX,typename DT::vec_pair_FX_long>(Qu0y);
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
// degreebound.  This routine only works over fields F.
template<class F, class vec_F, class FX, class FXY>
vector<FX> GSDecoder<F,vec_F,FX,FXY>::rr_findroots(
	const FXY &P, int degreebound)
{
    vector<int> pi;
    vector<int> Deg;
    vector<F> Coeff;
    vector<FXY> Q;
    vector<FX> res;
    int t = 1;

    FXY P0 = backShiftX(P, minX(P));
    pi.push_back(-1);
    Deg.push_back(-1);
    Coeff.push_back(F::zero());
    Q.push_back(P0);
    int u = 0; 

    if (degreebound < 0) {
	degreebound = degX(P0);
    }
    dfs(res, u, pi, Deg, Coeff, Q, t, degreebound);

    return res;
}


// Return a list of roots for y of the bivariate polynomial P(x,y).
// If degreebound >= 0, only return those roots with degree <=
// degreebound.  This routine handles the case where F is the
// integers mod p1*p2 as well as fields.  (But it handles the ring case
// in a specialization of this templated method.) This routine may also
// return some spurious values.
template<class F, class vec_F, class FX, class FXY>
vector<FX> GSDecoder<F,vec_F,FX,FXY>::findroots(const FXY &P, int degreebound)
{
    return rr_findroots(P, degreebound);
}

#endif
