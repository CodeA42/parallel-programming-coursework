#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 smallpt.cpp -o smallpt
#include <stdio.h>  // Usage: time ./smallpt 5000 && xv image.ppm
#pragma omp parallel
struct Vec
{
  double x, y, z; // position, also color (r,g,b)
  Vec(double x_ = 0, double y_ = 0, double z_ = 0)
  {
    x = x_;
    y = y_;
    z = z_;
  }
  Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
  Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
  Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
  Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
  Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
  double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
  Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

struct Ray
{
  Vec o, d;
  Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
  DIFF,
  SPEC,
  REFR
}; // material types, used in radiance()
struct Sphere
{
  double rad;  // radius
  Vec p, e, c; // position, emission, color
  Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
  Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
  double intersect(const Ray &r) const
  {                   // returns distance, 0 if nohit
    Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
    if (det < 0)
      return 0;
    else
      det = sqrt(det);
    return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
  }
};

// #include "extraScenes.h"
Sphere spheres[] = {
    // Scene: radius, position, emission, color, material
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),   // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),         // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),               // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),         // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), // Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),        // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),        // Glas
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF)     // Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray &r, double &t, int &id)
{
  double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
  for (int i = int(n); i--;)
    if ((d = spheres[i].intersect(r)) && d < t)
    {
      t = d;
      id = i;
    }
  return t < inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
  double t;   // distance to intersection
  int id = 0; // id of intersected object
  if (!intersect(r, t, id))
    return Vec();                  // if miss, return black
  const Sphere &obj = spheres[id]; // the hit object
  Vec x = r.o + r.d * t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
  double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                      : f.z; // max refl
  if (++depth > 5)
  {
    if (erand48(Xi) < p)
      f = f * (1 / p);
    else
      return obj.e;
  } // R.R.
  if (obj.refl == DIFF)
  { // Ideal DIFFUSE reflection
    double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
    Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
    Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
    return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
  }
  else if (obj.refl == SPEC) // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
  Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
  bool into = n.dot(nl) > 0;                // Ray from outside going in?
  double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
    return obj.e + f.mult(radiance(reflRay, depth, Xi));
  Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
  double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
  double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
  return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
                                         radiance(reflRay, depth, Xi) * RP
                                                     : radiance(Ray(x, tdir), depth, Xi) * TP)
                                  : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}

int main(int argc, char *argv[])
{
  int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
  Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
  Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
  for (int y = 0; y < h; y++)
  { // Loop over image rows
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
    for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < w; x++) // Loop cols
      for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)         // 2x2 subpixel rows
        for (int sx = 0; sx < 2; sx++, r = Vec())
        { // 2x2 subpixel cols
          for (int s = 0; s < samps; s++)
          {
            double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                    cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
            r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
          } // Camera rays are pushed ^^^^^ forward to start in interior
          c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
        }
  }
  FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}

// This code appears to be a small path tracer, implementing the Monte Carlo method for rendering realistic images of 3D scenes with lighting and reflection effects. It defines a few structs, classes, and functions that work together to calculate the color of each pixel of the image, given the scene and the observer's viewpoint. Here is a summary of what each part does:

// Vec struct: represents a vector in 3D space, with x, y, and z components. Also used to store colors, with each component representing the intensity of one of the primary colors (red, green, and blue).
// Ray struct: represents a ray of light, with an origin (o) and a direction (d).
// Refl_t enum: used to specify the type of material that a sphere is made of, which determines how it reflects light.
// Sphere struct: represents a sphere in the scene, with a radius, a position, an emission color (light emitted by the sphere), a color (reflective properties of the sphere), and a material type.
// clamp() function: clamps a value between 0 and 1.
// toInt() function: converts a floating-point value in the range 0-1 to an integer value in the range 0-255, using a gamma correction factor of 2.2.
// intersect() function: tests if a given ray intersects any of the spheres in the scene, and returns the distance to the closest intersection (if any), and the index of the closest sphere.
// radiance() function: calculates the color of a given ray by tracing it through the scene, bouncing it off of surfaces, and sampling the incoming light at each intersection point. This is the heart of the path tracer algorithm. The function takes a Ray object, a depth counter (to limit the recursion depth), and a random number generator seed (Xi). The algorithm works by first finding the closest intersection point of the ray with any sphere in the scene. If there is no intersection, the function returns black (no light). Otherwise, it calculates the color of the sphere at that point by computing the contribution of each light source in the scene, the reflection and transmission of the light by the sphere's material, and the indirect illumination caused by other objects in the scene. This is done by recursively calling radiance() with new rays that are randomly generated based on the material properties of the current sphere (i.e., its reflection and transmission coefficients). The recursion terminates when the maximum depth is reached or the incoming light becomes too dim. The final color of the ray is the sum of all the contributions calculated along its path. The function also implements Russian roulette, a technique to randomly terminate the recursion based on the probability of continuing (determined by the maximum reflection coefficient encountered along the path). This helps reduce noise and speed up the rendering process.

// Optimizing code for parallel performance requires identifying opportunities for parallelism and applying parallel programming techniques to exploit them. Here are some suggestions for optimizing the given code for parallel performance:

// Parallelize ray-tracing computation: The core computation in this code is the ray-sphere intersection test, which is done sequentially for each sphere. This computation can be parallelized by dividing the image into smaller regions and assigning each region to a separate thread for computation. This approach is known as domain decomposition, and it can significantly speed up the computation for large images.

// Use SIMD instructions: SIMD (Single Instruction Multiple Data) instructions can be used to perform the same operation on multiple data elements in parallel. The code can be optimized using SIMD instructions such as SSE (Streaming SIMD Extensions) or AVX (Advanced Vector Extensions), which can be used to perform vector operations in parallel.

// Use multi-threading for sphere intersection: Since each sphere intersection test is independent of the others, it can be executed in parallel using multiple threads. This can be achieved by using a thread pool and assigning each intersection test to a separate thread.

// Use cache optimization techniques: Cache optimization techniques such as cache blocking can be used to improve memory access patterns and reduce cache misses. This can improve performance by reducing the time spent waiting for memory accesses.

// Use loop unrolling: Loop unrolling can be used to reduce the overhead of loop control statements and improve the efficiency of the computation. This can be achieved by manually unrolling loops and performing multiple iterations of the loop in a single iteration.

// Use compiler optimizations: Compiler optimizations such as loop unrolling, vectorization, and inlining can be used to improve the performance of the code. These optimizations can be enabled using compiler flags such as -O3 or -Ofast.

// Use data parallelism: Data parallelism can be used to process data in parallel by dividing the data into smaller chunks and assigning each chunk to a separate thread for computation. This approach can be used to parallelize operations such as image filtering or color correction.

// Use task parallelism: Task parallelism can be used to parallelize tasks that have different execution times. This approach involves dividing the tasks into smaller sub-tasks and assigning each sub-task to a separate thread for computation. Task parallelism can be used to parallelize operations such as scene traversal or shading.

// Use OpenMP: OpenMP is a popular API for parallel programming in C and C++. It provides a simple and portable way to parallelize code using compiler directives. The code can be parallelized using OpenMP directives such as #pragma omp parallel for, which can be used to parallelize loops.

// Use CUDA: If the code is running on a GPU, it can be parallelized using CUDA, which is a parallel computing platform and programming model developed by NVIDIA. CUDA provides a way to write parallel programs that can be executed on NVIDIA GPUs, which can significantly speed up the computation for certain types of computations.

// Note that optimizing code for parallel performance can be challenging and requires careful consideration of the trade-offs between performance, scalability, and complexity. It is important to measure the performance of the optimized code and validate that it meets the performance requirements before deploying it in a production environment.
