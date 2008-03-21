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

#ifndef __RECOVER_H__
#define __RECOVER_H__

#include <vector>
#include <map>
#include <string>
#include <vec_ZZ_p.h>
#include <GF2EX.h>
#include "percyresult.h"
#include "FXY.h"
#include "rr_roots.h"
#include "gf28.h"

NTL_CLIENT

vector<PercyResult> EasyRecover(unsigned int bytes_per_word, unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const vec_ZZ_p &values, const vec_ZZ_p &indices);


template <class FX>
struct RecoveryPoly {
    RecoveryPoly(vector<unsigned short> G, FX phi) : G(G), phi(phi) {}
    vector<unsigned short> G;
    FX phi;
};

typedef pair<unsigned int, unsigned int> pairint;
#define Ccache_t map<pairint, F>

template <class F>
static F C(Ccache_t &Ccache, unsigned int n, unsigned int k)
{
    pairint nk;
    nk.first = n;
    nk.second = k;
    typename Ccache_t::const_iterator Ci = Ccache.find(nk);
    if (Ci != Ccache.end()) {
	return Ci->second;
    }
    
    F num, dem;
    num = 1;
    dem = 1;
    unsigned int i;

    for (i=0;i<k;++i) {
	num *= (n-i);
	dem *= (k-i);
    }

    num /= dem;
    Ccache[nk] = num;
    return num;
}

template <class F, class vec_F, class FX, class FXY>
class GSDecoder {
    public:
	GSDecoder() {}
	GSDecoder(const ZZ &p1, const ZZ &p2): p1(p1), p2(p2) {}

	vector<PercyResult> HardRecover(unsigned int bytes_per_word,
		unsigned short t, unsigned short h,
		const vector<PercyResult> &H,
		const vector<unsigned short> &goodservers,
		const vec_F &values, const vec_F &indices);

	static string append(const string &s, const F &wz,
		unsigned int bytes_per_word);

    private:
	static vector<unsigned short> intersect(
		const vector<unsigned short> &v1,
		const vector<unsigned short> &v2);

	vector< RecoveryPoly<FX> > findpolys(unsigned int k, unsigned int t,
		const vector<unsigned short> &goodservers,
		const vec_F& indices, const vec_F& shares);

	// Construct a bivariate polynomial P(x,y) such that, for any polynomial
	// f(x) of degree at most v that agrees with at least t of the given
	// points, (y-f(x)) is a factor of P(x,y).  This version is K"otter's
	// Interpolation Algorithm, as described in Section VII of R. J.
	// McEliece. The Guruswami-Sudan Decoding Algorithm for Reed-Solomon
	// Codes. IPN Progress Report 42-153, May 15, 2003.
	// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
	FXY interpolate_kotter(unsigned int v, unsigned int t,
		const vector<unsigned short> &goodservers,
		const vec_F &indices, const vec_F &shares);

	// Evaluate the (r,s)th Hasse mixed partial derivative of g at the point
	// (alpha, beta), which is:
	// \sum_{i,j} C(i,r) C(j,s) a_{i,j} alpha^{i-r} beta^{j-s}
	// = \sum_j C(j,s) beta^{j-s} \sum_i C(i,r) a_{i,j} alpha^{i-r}
	// where g(x,y) = \sum_{i,j} a_{i,j} x^i y^j
	// See page 14 of R. J. McEliece. The Guruswami-Sudan Decoding
	// Algorithm for Reed-Solomon Codes. IPN Progress Report 42-153, May 15,
	// 2003.  http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html
	F evalhasse(const FXY &g, unsigned int r, unsigned int s,
		F alpha, F beta, Ccache_t &Ccache);

	ZZ p1, p2;
};

typedef GSDecoder<ZZ_p, vec_ZZ_p, ZZ_pX, ZZ_pXY> GSDecoder_ZZ_p;
typedef GSDecoder<GF2E, vec_GF2E, GF2EX, GF2EXY> GSDecoder_GF2E;

// Compute the intersection of the two *strictly sorted* vectors v1 and v2
template<class F, class vec_F, class FX, class FXY>
vector<unsigned short> GSDecoder<F,vec_F,FX,FXY>::intersect(
	const vector<unsigned short> &v1,
	const vector<unsigned short> &v2)
{
    vector<unsigned short> intersection;
    vector<unsigned short>::const_iterator v1p, v2p;
    v1p = v1.begin();
    v2p = v2.begin();
    while (v1p != v1.end() && v2p != v2.end()) {
	if (*v1p == *v2p) {
	    // This is an element of the intersection
	    intersection.push_back(*v1p);
	    ++v1p;
	    ++v2p;
	} else if (*v1p > *v2p) {
	    ++v2p;
	} else {  // *v1p < *v2p
	    ++v1p;
	}
    }
    return intersection;
}

// Evaluate the (r,s)th Hasse mixed partial derivative of g at the point
// (alpha, beta), which is:
// \sum_{i,j} C(i,r) C(j,s) a_{i,j} alpha^{i-r} beta^{j-s}
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
    ZZ_p res;
    res = 0;
    int ydeg = deg(g);
    for (int j = ydeg; j >= (int)s; --j) {
	const ZZ_pX &gj = coeff(g,j);
	int xdeg = deg(gj);
	// Use Horner's method to evaluate the inner sum (which is the
	// coefficient of beta^{j-s})
	ZZ_p resx;
	resx = 0;
	for (int i = xdeg; i >= (int) r; --i) {
	    resx *= alpha;
	    resx += C(Ccache, i,r) * coeff(gj, i);
	}
	// Use Horner's method to accumulate the results into the outer
	// sum
	res *= beta;
	res += C(Ccache, j,s) * resx;
    }

    return res;
}

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

    unsigned int n = goodservers.size();

    // Compute the m and L parameters
    unsigned int m = 1 + (unsigned int)(floor( v*n / (t*t-v*n)));
    unsigned int L = (m*t - 1)/v;

    std::cerr << "Constructing (1," << v << ")-degree " << L*v << " polynomial...\n";
    std::cerr << "Estimated work: " << n * m * (m+1) / 2 * (L+1) << "\n";
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
    FXY P;
    //char *naiveenv = getenv("PIRC_NAIVE");
    //if (naiveenv && atoi(naiveenv)) {
	//P = interpolate_naive(k, t, goodservers, indices, shares);
    //} else {
	P = interpolate_kotter(k, t, goodservers, indices, shares);
    //}

    // cerr << "factor(poly(0";
    for(int j=0;j<=deg(P);++j) {
	FX x = coeff(P,j);
	for(int i=0; i<=deg(x); ++i) {
	    // cerr << " + " << coeff(x,i) << "*x^" << i << "*y^" << j;
	}
    }
    // cerr << ", [x,y], IntMod(" << p1*p2 << ")));\n\n";
    std::cerr << "Finding roots of resulting polynomial...\n";

    // It turns out that any polynomial phi(x) that we should be
    // returning (since phi is of degree at most k and agrees with the
    // input data on at least t points) is such that (y - phi(x)) is a
    // factor of P(x,y).  So we first find all roots for y of P(x,y)
    // which are polynomials of degree at most k.
    vec_ZZ_pX roots = findroots(P, k, p1, p2);
    // cerr << "roots = " << roots << "\n";

    // For each of these roots, check how many input points it agrees
    // with.  If it's at least t, add it to the list of polys to return.
    vector< RecoveryPoly<FX> > polys;
    unsigned int numroots = roots.length();
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
    vector< RecoveryPoly<FX> > polys = findpolys(t, h, goodservers, indices, values);
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

vector<PercyResult> EasyRecover_GF28(unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const GF28_Element *values, const GF28_Element *indices);

//template class GSDecoder<GF2E, vec_GF2E, GF2EX>;

#endif
