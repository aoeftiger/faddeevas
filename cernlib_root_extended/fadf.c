#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "fadf.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029  /* pi */
#endif

#define SQR(a) ((a)*(a))

static inline cnum_t
fexp(cnum_t z, int tauMi) // the Fourier expansion approximation
{
  assert(tauMi == 0 || tauMi == 1); // sanity check

  enum { maxN = 23 }; // number of summation terms
  static const num_t tauM[2] = { 12.0, 12.1 };
  static num_t aN[2][maxN]; // 0 -> 12.0, 1 -> 12.1
  static int done = 0;

  if (!done) {
    done = 1; 
    for (int n = 1; n <= maxN; n++) { // Fourier coefficients
      aN[0][n] = 2 * sqrt(M_PI) / tauM[0] * exp(-SQR(n) * SQR(M_PI)/SQR(tauM[0])); // tauM=12.0
      aN[1][n] = 2 * sqrt(M_PI) / tauM[1] * exp(-SQR(n) * SQR(M_PI)/SQR(tauM[1])); // tauM=12.1
    }
  }

  const cnum_t z1 = cexp(I * tauM[tauMi] * z); // define first value
  const cnum_t z2 = SQR(tauM[tauMi]) * SQR(z); // define second value
  
  int s = 1;
  cnum_t FE = sqrt(M_PI) / tauM[tauMi] * (1 - z1)/z2; // initiate FE

  for (int n = 1; n <= maxN; n++) {
    s = -s;
    FE += aN[tauMi][n] * (s * z1 - 1) / (SQR(n) * SQR(M_PI) - z2);
  }
  FE *= I * (SQR(tauM[tauMi])/sqrt(M_PI)) * z;

  return FE;
}

static inline cnum_t
contfr(cnum_t z) // the Laplace continued fraction approximation
{
  enum { aN = 8 }; // initial integer

  cnum_t CF = aN/2/z; // start computing from the last aN

  for (int n = aN-1; n >= 1; n--) {
    CF = n / (2*(z - CF));
  }
  CF = I/sqrt(M_PI) / (z - CF);

  return CF;
}

static inline cnum_t
smallim(cnum_t z) // approximation at small imag(z)
{
  num_t z_are = fabs(creal(z));
  cnum_t SIm;

  if (!(z_are <= 1e-3) && !(z_are > 4)) {
    SIm = fexp(z,1);
  }
  else if (z_are <= 1e-3) {
    SIm = (1 - SQR(z)) * (1 + 2*I/sqrt(M_PI) * z); 
  }
  else if (z_are > 4) {
    const num_t del = 1e-5; // assign delta
    cnum_t z_del = creal(z) + I*del; // incorporate delta
    cnum_t zr = cimag(z)/del; // assign ratio
    cnum_t ff_ref = fexp(z_del,0); // assign reference of the Faddeeva function
  
    // Finally, for |Re[z]| > 4 apply
    SIm = zr * creal(ff_ref) + (1 - zr) * exp(-SQR(creal(z_del))) + I*cimag(ff_ref);
  }
  else assert(NULL && "smallim: invalid region"); // sanity check

  return SIm;
}

cnum_t
fadf(cnum_t z)
{
  int is_neg = 0;
  if (cimag(z) < 0) { is_neg = 1 ; z = conj(z); }

  num_t z_abs = cabs(z);
  cnum_t FF = 0;

  if (!(z_abs > 15) && !(z_abs <= 15 && cimag(z) < 1e-5)) {
//    printf("fexp: ");
    FF = fexp(z,0); // internal area. This area is the most difficult for accurate computation
  }
  else if (z_abs > 15) {
//    printf("cont: ");
    FF = contfr(z); // external area
  }
  else if (z_abs <= 15) {
//    printf("smal: ");
    FF = smallim(z); // narrow band
  }
  else assert(NULL && "fadf: invalid region"); // sanity check

  // Convert for negative imag(z) values
  if (is_neg) FF = conj(2 * cexp(-SQR(z)) - FF);

  return FF;
}
