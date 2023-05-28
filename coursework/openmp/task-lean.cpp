#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 smallpt.cpp -o smallpt
#include <stdio.h>  // Usage: time ./smallpt 5000 && xv image.ppm
#include <omp.h>

/**
 * @brief Represents a vector in 3D space, with x, y, and z components.
 *
 * Also used to store colors, with each component representing the
 * intensity of one of the primary colors (red, green, and blue).
 */
struct Vec
{
  double x, y, z; // position, also color (r,g,b)
  Vec(double x_ = 0, double y_ = 0, double z_ = 0)
  {
    x = x_;
    y = y_;
    z = z_;
  }
  Vec operator+(const Vec &b) const
  {
    return Vec(x + b.x, y + b.y, z + b.z);
  }
  Vec operator-(const Vec &b) const
  {
    return Vec(x - b.x, y - b.y, z - b.z);
  }
  Vec operator*(double b) const
  {
    return Vec(x * b, y * b, z * b);
  }
  Vec mult(const Vec &b) const
  {
    return Vec(x * b.x, y * b.y, z * b.z);
  }
  Vec &norm()
  {
    return *this = *this * (1 / sqrt(x * x + y * y + z * z));
  }
  // cross:
  double dot(const Vec &b) const
  {
    return x * b.x + y * b.y + z * b.z;
  }
  Vec operator%(Vec &b)
  {
    return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
  }
};

/**
 * @brief represents a ray of light, with an origin (o) and a direction (d).
 */
struct Ray
{
  Vec origin;
  Vec direction;
  Ray(Vec origin_, Vec direction_) : origin(origin_), direction(direction_) {}
};

Vec radiance(const Ray &ray, int depth, unsigned short *Xi);

/**
 * @brief Used to specify the type of material that a sphere
 * is made of, which determines how it reflects light.
 */
enum Refl_t
{
  DIFF,
  SPEC,
  REFR
}; // material types, used in radiance()

/**
 * @brief Represents a sphere in the scene, with a radius, a position,
 * an emission color (light emitted by the sphere), a color
 * (reflective properties of the sphere), and a material type.
 */
struct Sphere
{
  double rad;            // radius
  Vec position;          // position
  Vec emission;          // emission
  Vec color;             // color
  Refl_t reflectionType; // reflection type (DIFFuse, SPECular, REFRactive)
  Sphere(double rad_, Vec position_, Vec emission_, Vec color_, Refl_t refl_) : rad(rad_), position(position_), emission(emission_), color(color_), reflectionType(refl_) {}
  double intersect(const Ray &ray) const // returns distance, 0 if nohit
  {
    Vec op = position - ray.origin; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t;
    double eps = 1e-4;
    double b = op.dot(ray.direction);
    double det = b * b - op.dot(op) + rad * rad;
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

/**
 * @brief Clamps a value between 0 and 1.
 *
 * @param x Value to be clamped
 * @return double
 */
inline double clamp(double x)
{
  if (x < 0)
  {
    return 0;
  }
  else if (x > 1)
  {
    return 1;
  }
  return x;
}

/**
 * @brief Converts a floating-point value in the range 0-1
 * to an integer value in the range 0-255, using a
 * gamma correction factor of 2.2.
 *
 * @param x Value in the range 0-1
 * @return int
 */
inline int toInt(double x)
{
  return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}
/**
 *  @brief tests if a given ray intersects any of the spheres in the scene,
 *  and returns the distance to the closest intersection (if any),
 *  and the index of the closest sphere.
 **/
inline bool intersect(const Ray &ray, double &closestIntersection, int &closestSphereId)
{
  double numberOfSpheres = sizeof(spheres) / sizeof(Sphere);
  double distance;
  double infinity = 1e20;
  closestIntersection = 1e20;
  for (int id = int(numberOfSpheres); id--;)
    if ((distance = spheres[id].intersect(ray)) && distance < closestIntersection)
    {
      closestIntersection = distance;
      closestSphereId = id;
    }
  return closestIntersection < infinity;
}

inline Vec calculateIdealDiffuseReflection(unsigned short *Xi, Vec nl, const Sphere &obj, Vec f, Vec x, int depth)
{
  double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
  Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
  Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
  return obj.emission + f.mult(radiance(Ray(x, d), depth, Xi));
}

inline Vec getRayColor(int depth, unsigned short *seed, double P, Ray reflRay, double RP, Vec x, Vec tdir, double TP, double Re, double Tr)
{
  if (depth > 2)
  {
    if (erand48(seed) < P) // Russian roulette
    {
      return radiance(reflRay, depth, seed) * RP;
    }
    else
    {
      return radiance(Ray(x, tdir), depth, seed) * TP;
    }
  }
  else
  {
    return radiance(reflRay, depth, seed) * Re + radiance(Ray(x, tdir), depth, seed) * Tr;
  }
}

/**
 * @brief Calculates the color of a given ray by tracing it through the scene,
 * bouncing it off of surfaces, and sampling the incoming light at each
 * intersection point. This is the heart of the path tracer algorithm.
 * The function takes a Ray object, a depth counter (to limit the recursion depth),
 * and a random number generator seed (Xi). The algorithm works by first
 * finding the closest intersection point of the ray with any sphere in the scene.
 * If there is no intersection, the function returns black (no light).
 * Otherwise, it calculates the color of the sphere at that point by computing
 * the contribution of each light source in the scene, the reflection and
 * transmission of the light by the sphere's material, and the indirect
 * illumination caused by other objects in the scene. This is done by recursively
 * calling radiance() with new rays that are randomly generated based on the material
 * properties of the current sphere (i.e., its reflection and transmission
 * coefficients). The recursion terminates when the maximum depth is reached or the
 * incoming light becomes too dim. The final color of the ray is the sum of all the
 * contributions calculated along its path. The function also implements Russian
 * roulette, a technique to randomly terminate the recursion based on the probability
 * of continuing (determined by the maximum reflection coefficient encountered along
 * the path). This helps reduce noise and speed up the rendering process.
 *
 * @param ray Ray to base the calculations for
 * @param depth Depth Counter
 * @param seed Seed used for the calculations
 * @return Vec
 */
Vec radiance(const Ray &ray, int depth, unsigned short *seed)
{
  double distToIntersection;
  int idOfIntersectedObj;
  if (!intersect(ray, distToIntersection, idOfIntersectedObj))
    return Vec();                                                // if miss, return black
  const Sphere &intersectedSphere = spheres[idOfIntersectedObj]; // the hit object
  Vec x = ray.origin + ray.direction * distToIntersection;
  Vec n = (x - intersectedSphere.position).norm();
  Vec nl = n.dot(ray.direction) < 0 ? n : n * -1, f = intersectedSphere.color;
  double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                      : f.z; // max refl
  if (++depth > 5)
  {
    if (erand48(seed) < p)
    {
      f = f * (1 / p);
    }
    else
    {
      return intersectedSphere.emission;
    }
  }                                             // R.R.
  if (intersectedSphere.reflectionType == DIFF) // Ideal DIFFUSE reflection
  {
    return calculateIdealDiffuseReflection(seed, nl, intersectedSphere, f, x, depth);
  }
  else if (intersectedSphere.reflectionType == SPEC) // Ideal SPECULAR reflection
  {
    return intersectedSphere.emission + f.mult(radiance(Ray(x, ray.direction - n * 2 * n.dot(ray.direction)), depth, seed));
  }
  Ray reflRay(x, ray.direction - n * 2 * n.dot(ray.direction)); // Ideal dielectric REFRACTION
  bool into = n.dot(nl) > 0;                                    // Ray from outside going in?
  double nc = 1;
  double nt = 1.5;
  double nnt = into ? nc / nt : nt / nc;
  double ddn = ray.direction.dot(nl);
  double cos2t;
  if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
  {
    return intersectedSphere.emission + f.mult(radiance(reflRay, depth, seed));
  }
  Vec tdir = (ray.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
  double a = nt - nc;
  double b = nt + nc;
  double R0 = a * a / (b * b);
  double c = 1 - (into ? -ddn : tdir.dot(n));
  double Re = R0 + (1 - R0) * c * c * c * c * c;
  double Tr = 1 - Re;
  double P = .25 + .5 * Re;
  double RP = Re / P;
  double TP = Tr / (1 - P);
  Vec rayColor = getRayColor(depth, seed, P, reflRay, RP, x, tdir, TP, Re, Tr);
  return intersectedSphere.emission + f.mult(rayColor);
}

int main(int argc, char *argv[])
{
  int width = 1024;
  int height = 768;
  int samples = argc == 2 ? atoi(argv[1]) / 4 : 1;
  Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
  Vec cx = Vec(width * .5135 / height);
  Vec cy = (cx % cam.direction).norm() * .5135;
  Vec r;
  Vec *c = new Vec[width * height];
  for (int y = 0; y < height; y++) // Loop over image rows
  {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples * 4, 100. * y / (height - 1));
    for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < width; x++) // Loop cols
    {
      // Subpixel y
      int sy;
      int i;
      for (sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) // 2x2 subpixel rows
      {
        // Subpixel x
        int sx;
        for (sx = 0; sx < 2; sx++, r = Vec()) // 2x2 subpixel cols
        {
          for (int s = 0; s < samples; s++)
          {
            double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
            double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
            Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                    cy * (((sy + .5 + dy) / 2 + y) / height - .5) + cam.direction;
            r = r + radiance(Ray(cam.origin + d * 140, d.norm()), 0, Xi) * (1. / samples);
          } // Camera rays are pushed ^^^^^ forward to start in interior
          c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
        }
      }
    }
  }
  FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
  for (int i = 0; i < width * height; i++)
  {
    fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
  }
}
