#ifndef __GF28_H__
#define __GF28_H__

typedef unsigned char GF28_Element;

/* Use the same representation as AES */
static const unsigned char GF28_Generator = 0x1b;

extern const GF28_Element GF28_mult_table[256][256];

extern const GF28_Element GF28_inv_table[256];

GF28_Element evalpoly_GF28(GF28_Element *coeffs, unsigned short degree,
	GF28_Element index);

// return f(alpha), where f is the polynomial of degree numpoints-1 such
// that f(indices[i]) = values[i] for i=0..(numpoints-1)
GF28_Element interpolate_GF28(const GF28_Element *indices,
	const GF28_Element *values, unsigned short numpoints,
	const GF28_Element alpha);

#endif
