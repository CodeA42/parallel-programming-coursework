#include "grayscale.h"

// Чрез указателя а ще имаме достъп до входните данни, а в чрез указателя c ще запазваме трансформираните данните
// така, че след това тези трансформирани данни да бъдат копирани в паметта на хоста, за да може да се запазят в ново изображение
__kernel void grayscale(__global PPMPixel *a, __global PPMPixel *c, const unsigned int n) {                                                             
     // Идентификатор на нишката
     int id = get_global_id(0);
     if (id < n) {
          double f = 1;
          double l = 0.3 * a[id].red + 0.6 * a[id].green + 0.1 * a[id].blue;
          c[id].red = a[id].red + f * (l - a[id].red);
          c[id].green = a[id].green + f * (l - a[id].green);
          c[id].blue = a[id].blue + f * (l - a[id].blue);
     }
}   