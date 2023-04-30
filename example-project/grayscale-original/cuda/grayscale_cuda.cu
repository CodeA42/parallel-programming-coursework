// nvcc grayscale_cuda.cu -o grayscale_cuda

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

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

// kernel функция
__global__ void grayscale(int n, PPMPixel *source, PPMPixel *target) {
     int id = blockIdx.x*blockDim.x + threadIdx.x;
     if (id < n) {
          double f = 1;
          double l = 0.3 * source[id].red + 0.6 * source[id].green + 0.1 * source[id].blue;
          target[id].red = source[id].red + f * (l - source[id].red);
          target[id].green = source[id].green + f * (l - source[id].green);
          target[id].blue = source[id].blue + f * (l - source[id].blue);
     }
}

void ChangeColorPPM(PPMImage *img) {
     struct timeval tval_before, tval_after, tval_result;

     int N = img->x * img->y;
     PPMPixel *source, *d_source, *target, *d_target;

     // Заделяне на памет за входните и трансформираните данни в хоста
     source = (PPMPixel*)malloc(N*sizeof(PPMPixel));
     target = (PPMPixel*)malloc(N*sizeof(PPMPixel));

     // Заделяне на памет за входните и трансфомираните данни в cuda устройството
     cudaMalloc(&d_source, N*sizeof(PPMPixel));
     cudaMalloc(&d_target, N*sizeof(PPMPixel));

     // Измерване колко време отнема копирането на данните в устройството
     gettimeofday(&tval_before, NULL);

     // Копиране на данните от хоста в cuda устройството
     cudaMemcpy(d_source, img->data, N*sizeof(PPMPixel), cudaMemcpyHostToDevice);
     cudaMemcpy(d_target, img->data, N*sizeof(PPMPixel), cudaMemcpyHostToDevice);

     gettimeofday(&tval_after, NULL);
     
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за копиране на масивите в устройството\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     
     // Измерване колко време отнема изпънението на kernel функцията
     gettimeofday(&tval_before, NULL);

     // Изпълнение на kernel функцияа
     grayscale<<<(N+383)/384, 384>>>(N, d_source, d_target);
     gettimeofday(&tval_after, NULL);
     
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за изпълнението на kernel функцията\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     // printf("error: %s\n", cudaGetErrorString(cudaGetLastError()));

     gettimeofday(&tval_before, NULL);

     // Копиране на данните от устройството в хоста, за да можем да ги прочетем и запишем в изображението
     cudaMemcpy(target, d_target, N*sizeof(PPMPixel), cudaMemcpyDeviceToHost);
     gettimeofday(&tval_after, NULL);
     // Измерване колко време отнема копирането на данните в хоста
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за копиране на данните обратно в хоста\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     // printf("error: %s\n", cudaGetErrorString(cudaGetLastError()));

     gettimeofday(&tval_before, NULL);

     // Записване на получение трансформирани данни в структурата от данни, която се използва от записващата функция
     for (int i = 0; i < N; i++) {
          img->data[i] = target[i];
     }
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за копиране на данните обратно в масива с пикселите\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     // printf("error: %s\n", cudaGetErrorString(cudaGetLastError()));

     // Освобождаване на ресурси
     cudaFree(d_source);
     cudaFree(d_target);
     free(source);
     free(target);
}

int main(void) {
     PPMImage *image;
     struct timeval tval_before, tval_after, tval_result;

     gettimeofday(&tval_before, NULL);
     image = ReadPPM("image.ppm");
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld     секунди за четене на данните от изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

     gettimeofday(&tval_before, NULL);
     ChangeColorPPM(image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  --------\n");
     printf("  %ld.%06ld   секунди за обработка на данните от изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

     gettimeofday(&tval_before, NULL);
     WritePPM("grayscale_cuda_result.ppm", image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld     секунди за запис на данните в изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}
