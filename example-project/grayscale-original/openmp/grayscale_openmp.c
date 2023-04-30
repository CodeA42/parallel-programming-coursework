// gcc -x c -g -fopenmp grayscale_openmp.c -o grayscale_openmp
// /usr/local/bin/gcc-10 -x c -g -fopenmp grayscale_openmp.c -o grayscale_openmp

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

typedef struct {
     unsigned char red, green, blue;
} PPMPixel;

typedef struct {
     int x, y;
     PPMPixel *data;
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
     img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));

     if (!img) {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }

     // зареждане на данните за всеки пиксел
     if (fread(img->data, 3 * img->x, img->y, fp) != img->y) {
          fprintf(stderr, "Error loading image '%s'\n", filename);
          exit(1);
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
     fwrite(img->data, 3 * img->x, img->y, fp);
     fclose(fp);
}

void ChangeColorPPM(PPMImage *img) {
     // Искане на 4 нишки от системата
     omp_set_num_threads(4);

     if (img) {
          // Задаване на паралелен блок за for цикъл като размера на блоковете,
          // които ще бъдат създадени се определят по време на компилация
          #pragma omp parallel for schedule(static)
          for (int i = 0; i < img->x * img->y; i++) {
               // Внасяне на локалните променливи вътре в паралелния блок за да станат частни за отделните нишки
               int r, g, b, f = 1;
               double l = 0.3 * img->data[i].red + 0.6 * img->data[i].green + 0.1 * img->data[i].blue;

               img->data[i].red = img->data[i].red + f * (l - img->data[i].red);
               img->data[i].green = img->data[i].green + f * (l - img->data[i].green);
               img->data[i].blue = img->data[i].blue + f * (l - img->data[i].blue);
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
     WritePPM("grayscale_openmp_result.ppm", image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld секунди за запис на данните в изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}
