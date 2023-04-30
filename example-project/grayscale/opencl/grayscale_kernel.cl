#include "grayscale.h"
__kernel void grayscale(__global PPMPixel *a, __global PPMPixel *c, const unsigned int n) {
     int id = get_global_id(0);
     if (id < n) {
          double f = 1;
          double l = 0.3 * a[id].red + 0.6 * a[id].green + 0.1 * a[id].blue;
          c[id].red = a[id].red + f * (l - a[id].red);
          c[id].green = a[id].green + f * (l - a[id].green);
          c[id].blue = a[id].blue + f * (l - a[id].blue);
     }
}   