#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
typedef struct
{
     unsigned char red, green, blue;
} PPMPixel;
typedef struct
{
     int x, y;
     PPMPixel *data;
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
     img->data = (PPMPixel *)malloc(img->x * img->y * sizeof(PPMPixel));

     if (!img)
     {
          fprintf(stderr, "Unable to allocate memory\n");
          exit(1);
     }
     if (fread(img->data, 3 * img->x, img->y, fp) != img->y)
     {
          fprintf(stderr, "Error loading image '%s'\n", filename);
          exit(1);
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
     fwrite(img->data, 3 * img->x, img->y, fp);
     fclose(fp);
}
__global__ void grayscale(int n, PPMPixel *source, PPMPixel *target)
{
     int id = blockIdx.x * blockDim.x + threadIdx.x;
     if (id < n)
     {
          double f = 1;
          double l = 0.3 * source[id].red + 0.6 * source[id].green + 0.1 * source[id].blue;
          target[id].red = source[id].red + f * (l - source[id].red);
          target[id].green = source[id].green + f * (l - source[id].green);
          target[id].blue = source[id].blue + f * (l - source[id].blue);
     }
}
void ChangeColorPPM(PPMImage *img)
{
     struct timeval tval_before, tval_after, tval_result;
     int N = img->x * img->y;
     PPMPixel *source, *d_source, *target, *d_target;
     source = (PPMPixel *)malloc(N * sizeof(PPMPixel));
     target = (PPMPixel *)malloc(N * sizeof(PPMPixel));
     cudaMalloc(&d_source, N * sizeof(PPMPixel));
     cudaMalloc(&d_target, N * sizeof(PPMPixel));
     gettimeofday(&tval_before, NULL);
     cudaMemcpy(d_source, img->data, N * sizeof(PPMPixel), cudaMemcpyHostToDevice);
     cudaMemcpy(d_target, img->data, N * sizeof(PPMPixel), cudaMemcpyHostToDevice);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за копиране на масивите в устройството\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     gettimeofday(&tval_before, NULL);
     grayscale<<<(N + 383) / 384, 384>>>(N, d_source, d_target);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за изпълнението на kernel функцията\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     gettimeofday(&tval_before, NULL);
     cudaMemcpy(target, d_target, N * sizeof(PPMPixel), cudaMemcpyDeviceToHost);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за копиране на данните обратно в хоста\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     gettimeofday(&tval_before, NULL);
     for (int i = 0; i < N; i++)
     {
          img->data[i] = target[i];
     }
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("  %ld.%06ld   секунди за копиране на данните обратно в масива с пикселите\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     cudaFree(d_source);
     cudaFree(d_target);
     free(source);
     free(target);
}
int main(void)
{
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
