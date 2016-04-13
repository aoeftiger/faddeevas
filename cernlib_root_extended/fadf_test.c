// gcc -std=c99 -W -Wall -Wextra -pedantic -O3 -ffast-math -ftree-vectorize -flto fadf.c fadf_test.c -o fadf
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "fadf.h"

#define CNUM(re,im) (* (cnum_t*) & (num_t[2]) { (re), (im) })

int main(int argc, char *argv[])
{
  num_t step = argc >= 2 ? strtod(argv[1],0) : 0.05;
  num_t stop = argc >= 3 ? strtod(argv[2],0) : 20;

  cnum_t z = 0;
  for (num_t x = 0.0; x <= stop; x += step) {
    for (num_t y = 0.0; y <= stop; y += step) {
      z = fadf(CNUM(x,y));
      printf("fadf(%.2f, %.2f) = %.17e, %.17e\n", x, y, creal(z), cimag(z));    
    }  
    printf("\n");
  }

  return 0;
}