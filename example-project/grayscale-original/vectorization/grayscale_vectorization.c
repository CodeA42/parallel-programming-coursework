// gcc -x c -g -mavx grayscale_vectorization.c -o grayscale_vectorization

#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <sys/time.h>

// Структурата е трансформирана така, че лесно да се работи със SIMD функции
// SIMD функциите, които се използват, работят с параметър - указател към масив от стойности
typedef struct {
     int x, y;
     float *red;
     float *green;
     float *blue;
} PPMImage;

#define CREATOR "ParallelProgrammer"
#define RGB_COMPONENT_COLOR 255

static PPMImage *ReadPPM(const char *filename) {
     char buff[16];
     PPMImage *img;
     FILE *fp;
     int c, rgb_comp_color;

     // зареждане на PPM файла за четене
     fp = fopen(filename, "rb");
     if (!fp) {
          fprintf(stderr, "Unable to open file '%s'\n", filename);
          exit(1);
     }

     // четене на метаданните от него 
     if (!fgets(buff, sizeof(buff), fp)) {
          perror(filename);
          exit(1);
     }

     // проверка на метаданните
     if (buff[0] != 'P' || buff[1] != '6') {
          fprintf(stderr, "Invalid image format (must be 'P6')\n");
          exit(1);
     }

     // заделяне на памет
     img = (PPMImage *)malloc(sizeof(PPMImage));
     if (!img) {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }

     // проверка за коментари вътре в самото изображение
     c = getc(fp);
     while (c == '#') {
          while (getc(fp) != '\n') ;
          c = getc(fp);
     }

     ungetc(c, fp);
     // проверка на данните за размера на изображението
     if (fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
          fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
          exit(1);
     }

     // проверка на RGB компонента
     if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
          fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
          exit(1);
     }

     // проверка на размерността на RGB компонента
     if (rgb_comp_color!= RGB_COMPONENT_COLOR) {
          fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
          exit(1);
     }

     while (fgetc(fp) != '\n') ;
     // заделяне на памет за информацията във всеки пиксел
     img->red = (float*)malloc(img->x * img->y * sizeof(float));
     img->green = (float*)malloc(img->x * img->y * sizeof(float));
     img->blue = (float*)malloc(img->x * img->y * sizeof(float));

     if (!img) {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }

     // зареждане на данните за всеки пиксел
     for (int i = 0; i < img->x * img->y; i++) {
          unsigned char color[3];

          fread(color, 1, 3, fp);
          img->red[i] = color[0];
          img->green[i] = color[1];
          img->blue[i] = color[2];
    }

     fclose(fp);
     return img;
}
void WritePPM(const char *filename, PPMImage *img) {
     FILE *fp;
     // отваряне на файл в режим за писане
     fp = fopen(filename, "wb");
     if (!fp) {
          fprintf(stderr, "Unable to open file '%s'\n", filename);
          exit(1);
     }

     // записване на метаданни за типа на PPM изображението
     fprintf(fp, "P6\n");

     // запис на коментари
     fprintf(fp, "# Created by %s\n", CREATOR);

     // запис на размера на изображението
     fprintf(fp, "%d %d\n",img->x,img->y);

     // запис на размерността на RGB компонента
     fprintf(fp, "%d\n", RGB_COMPONENT_COLOR);

     // запис на данните за пикселите от изображението
     for (int i = 0; i < img->x * img->y; i++) {
           unsigned char color[3];
               
          color[0] = img->red[i];
          color[1] = img->green[i];
          color[2] = img->blue[i];

          fwrite(color, 1, 3, fp);
    }
     fclose(fp);
}

void ChangeColorPPM(PPMImage *img) {
     int i;
     if (img) {
          // 128 битова променлива за 4 32-битови float променливи
          const __m128 f = _mm_set1_ps(1.0);

          // 128 битови променливи за 4 32-битови float стойности
          __m128 l, r, g, b;
          for (i = 0; i < img->x * img->y; i += 4) {
               r = _mm_loadu_ps(&img->red[i]);
               g = _mm_loadu_ps(&img->green[i]);
               b = _mm_loadu_ps(&img->blue[i]);

               // printf("scalar: %f %f %f %f %f %f %f %f \n", &img->red[0], img->red[1], img->red[2], img->red[3], img->red[4], img->red[5], img->red[6], img->red[7]);
               // printf("vector: %f %f %f %f %f %f %f %f \n", r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);

               const __m128 red_comp = _mm_set1_ps(0.3);
               const __m128 green_comp = _mm_set1_ps(0.6);
               const __m128 blue_comp = _mm_set1_ps(0.1);

               const __m128 red_gray = _mm_mul_ps(red_comp, r);
               const __m128 green_gray = _mm_mul_ps(green_comp, g);
               const __m128 blue_gray = _mm_mul_ps(blue_comp, b);

               // l = 0.3 * img->data[i].red + 0.6 * img->data[i].green + 0.1 * img->data[i].blue;
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

int main() {
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