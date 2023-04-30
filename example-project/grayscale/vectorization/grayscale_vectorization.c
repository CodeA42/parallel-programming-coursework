#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <sys/time.h>
typedef struct
{
     int x, y;
     float *red;
     float *green;
     float *blue;
} PPMImage;
#define CREATOR "ParallelProgrammer"
#define RGB_COMPONENT_COLOR 255
static PPMImage *ReadPPM(const char *filename)
{
     char buff[16];
     PPMImage *img;
     FILE *fp;
     int c, rgb_comp_color;
     fp = fopen(filename, "rb");
     if (!fp)
     {
          fprintf(stderr, "Unable to open file '%s'\n", filename);
          exit(1);
     }
     if (!fgets(buff, sizeof(buff), fp))
     {
          perror(filename);
          exit(1);
     }
     if (buff[0] != 'P' || buff[1] != '6')
     {
          fprintf(stderr, "Invalid image format (must be 'P6')\n");
          exit(1);
     }
     img = (PPMImage *)malloc(sizeof(PPMImage));
     if (!img)
     {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }
     c = getc(fp);
     while (c == '#')
     {
          while (getc(fp) != '\n')
               ;
          c = getc(fp);
     }
     ungetc(c, fp);
     if (fscanf(fp, "%d %d", &img->x, &img->y) != 2)
     {
          fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
          exit(1);
     }
     if (fscanf(fp, "%d", &rgb_comp_color) != 1)
     {
          fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
          exit(1);
     }
     if (rgb_comp_color != RGB_COMPONENT_COLOR)
     {
          fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
          exit(1);
     }
     while (fgetc(fp) != '\n')
          ;
     img->red = (float *)malloc(img->x * img->y * sizeof(float));
     img->green = (float *)malloc(img->x * img->y * sizeof(float));
     img->blue = (float *)malloc(img->x * img->y * sizeof(float));
     if (!img)
     {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }
     for (int i = 0; i < img->x * img->y; i++)
     {
          unsigned char color[3];
          fread(color, 1, 3, fp);
          img->red[i] = color[0];
          img->green[i] = color[1];
          img->blue[i] = color[2];
     }
     fclose(fp);
     return img;
}
void WritePPM(const char *filename, PPMImage *img)
{
     FILE *fp;
     fp = fopen(filename, "wb");
     if (!fp)
     {
          fprintf(stderr, "Unable to open file '%s'\n", filename);
          exit(1);
     }
     fprintf(fp, "P6\n");
     fprintf(fp, "# Created by %s\n", CREATOR);
     fprintf(fp, "%d %d\n", img->x, img->y);
     fprintf(fp, "%d\n", RGB_COMPONENT_COLOR);
     for (int i = 0; i < img->x * img->y; i++)
     {
          unsigned char color[3];
          color[0] = img->red[i];
          color[1] = img->green[i];
          color[2] = img->blue[i];
          fwrite(color, 1, 3, fp);
     }
     fclose(fp);
}
void ChangeColorPPM(PPMImage *img)
{
     int i;
     if (img)
     {
          const __m128 f = _mm_set1_ps(1.0);
          __m128 l, r, g, b;
          for (i = 0; i < img->x * img->y; i += 4)
          {
               r = _mm_loadu_ps(&img->red[i]);
               g = _mm_loadu_ps(&img->green[i]);
               b = _mm_loadu_ps(&img->blue[i]);
               const __m128 red_comp = _mm_set1_ps(0.3);
               const __m128 green_comp = _mm_set1_ps(0.6);
               const __m128 blue_comp = _mm_set1_ps(0.1);
               const __m128 red_gray = _mm_mul_ps(red_comp, r);
               const __m128 green_gray = _mm_mul_ps(green_comp, g);
               const __m128 blue_gray = _mm_mul_ps(blue_comp, b);
               l = _mm_add_ps(_mm_add_ps(_mm_mul_ps(red_comp, red_gray), _mm_mul_ps(green_comp, green_gray)), _mm_mul_ps(blue_comp, blue_gray));
               r = _mm_add_ps(_mm_mul_ps(_mm_sub_ps(l, r), f), r);
               g = _mm_add_ps(_mm_mul_ps(_mm_sub_ps(l, g), f), g);
               b = _mm_add_ps(_mm_mul_ps(_mm_sub_ps(l, b), f), b);
               _mm_store_ps(&img->red[i], r);
               _mm_store_ps(&img->green[i], g);
               _mm_store_ps(&img->blue[i], b);
          }
     }
}
int main()
{
     PPMImage *image;
     struct timeval tval_before, tval_after, tval_result;
     gettimeofday(&tval_before, NULL);
     image = ReadPPM("image.ppm");
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld секунди за четене на данните от изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     gettimeofday(&tval_before, NULL);
     ChangeColorPPM(image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld секунди за обработка на данните от изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     gettimeofday(&tval_before, NULL);
     WritePPM("grayscale_vectorized_result.ppm", image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld секунди за запис на данните в изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}