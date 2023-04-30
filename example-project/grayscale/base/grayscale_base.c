// gcc -x c -g grayscale_base.c -o grayscale_base

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// Структура с данни за цветовата стойност на пиксел
typedef struct {
     unsigned char red, green, blue;
} PPMPixel;

// Структура, която съдържа данните за ширина и височина на изображението,
// както и указател към данните за пикселите
typedef struct {
     int x, y;
     PPMPixel *data;
} PPMImage;

// Записване на коментар в изображението
#define CREATOR "ParallelProgrammer"

// Максималната стойност на цвят, която може да се съдържа в изображението
#define RGB_COMPONENT_COLOR 255

// Четене на данните от изображението
static PPMImage *ReadPPM(const char *filename) {
     // Формат на PPM P6 изображение:
     // P6
     // 3 2
     // 255
     // # Частта отгоре е хедър
     // # "P6" е типа на PPM изображението
     // # "3 2" е широчината и височината на изображението в пиксели
     // # "255" е максималната стойност за всеки цвят
     // # Частта отдолу е данните за пикселите чрез RGB тройки
     // 255   0   0  \ # червен
     // 0    255  0  \ # зелен
     // 0     0  255 \ # син
     // 255  255  0  \ # жълт
     // 255  255 255 \ # бял
     //  0    0   0  \ # черен

     // Буфер, който съдържа метаданните на изображението
     char buff[16];

     // Указател към изображението, чиито данни ще бъдат модифицирани
     PPMImage *img;

     // Указател към физическия файл, който ще бъде модифициран
     FILE *fp;

     // Помощни променливи при четенето на данните
     int c, rgb_comp_color;

     // Зареждане на PPM файла за четене
     fp = fopen(filename, "rb");
     if (!fp) {
          fprintf(stderr, "Unable to open file '%s'\n", filename);
          exit(1);
     }

     // Четене на метаданните от него 
     if (!fgets(buff, sizeof(buff), fp)) {
          perror(filename);
          exit(1);
     }

     // Проверка на метаданните - дали форматът е P6, тъй като четенето и записа е съобразено с него
     if (buff[0] != 'P' || buff[1] != '6') {
          fprintf(stderr, "Invalid image format (must be 'P6')\n");
          exit(1);
     }

     // Заделяне на памет за данните от изображението
     img = (PPMImage *)malloc(sizeof(PPMImage));
     if (!img) {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }

     // Проверка за коментари вътре в самото изображение
     c = getc(fp);
     while (c == '#') {
          while (getc(fp) != '\n') ;
          c = getc(fp);
     }

     ungetc(c, fp);

     // Проверка на данните за размера на изображението
     if (fscanf(fp, "%d %d", &img->x, &img->y) != 2) {
          fprintf(stderr, "Invalid image size (error loading '%s')\n", filename);
          exit(1);
     }

     // Проверка на RGB компонента
     if (fscanf(fp, "%d", &rgb_comp_color) != 1) {
          fprintf(stderr, "Invalid rgb component (error loading '%s')\n", filename);
          exit(1);
     }

     // Проверка на размерността на RGB компонента
     if (rgb_comp_color!= RGB_COMPONENT_COLOR) {
          fprintf(stderr, "'%s' does not have 8-bits components\n", filename);
          exit(1);
     }

     while (fgetc(fp) != '\n') ;

     // Заделяне на памет за информацията за всеки пиксел
     img->data = (PPMPixel*)malloc(img->x * img->y * sizeof(PPMPixel));

     if (!img) {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }

     // Зареждане на данните за всеки пиксел
     if (fread(img->data, 3 * img->x, img->y, fp) != img->y) {
          fprintf(stderr, "Error loading image '%s'\n", filename);
          exit(1);
     }

     // Затваряне на изображението
     fclose(fp);
     return img;
}

// Запис на трансформираните данни във файл
void WritePPM(const char *filename, PPMImage *img) {
     FILE *fp;

     // Отваряне на файл в режим за писане
     fp = fopen(filename, "wb");
     if (!fp) {
          fprintf(stderr, "Unable to open file '%s'\n", filename);
          exit(1);
     }

     // Записване на метаданни за типа на PPM изображението
     fprintf(fp, "P6\n");

     // Запис на коментари
     fprintf(fp, "# Created by %s\n", CREATOR);

     // Запис на размера на изображението
     fprintf(fp, "%d %d\n",img->x,img->y);

     // Запис на размерността на RGB компонента
     fprintf(fp, "%d\n", RGB_COMPONENT_COLOR);

     // Запис на (трансформираните) данни за пикселите от изображението
     fwrite(img->data, 3 * img->x, img->y, fp);
     fclose(fp);
}

// Функция за трансформация
void ChangeColorPPM(PPMImage *img) {
     int i;
     if (img) {
          int r, g, b;
          // Стойност на намаляване на наситеността (saturation) - стойност 1 за 100%
          int f = 1;
          double l;
          for (i = 0; i < img->x * img->y; i++) {
               l = 0.3 * img->data[i].red + 0.6 * img->data[i].green + 0.1 * img->data[i].blue;

               // Намаляване на наситеността 
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
     WritePPM("grayscale_base_result.ppm", image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld секунди за запис на данните в изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}
