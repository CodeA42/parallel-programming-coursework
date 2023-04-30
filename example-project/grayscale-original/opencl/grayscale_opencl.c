// gcc -x c -g -framework OpenCL grayscale_opencl.c -o grayscale_opencl
// gcc -x c -g grayscale_opencl.c -o grayscale_opencl -lOpenCL -lm
// clinfo

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#ifdef __APPLE__
     #include <OpenCL/cl.h>
#else
     #include <CL/cl.h>
#endif
#include "grayscale.h"

#define CREATOR "ParallelProgrammer"
#define RGB_COMPONENT_COLOR 255

// Проверка за грешки от OpenCL
#define CL_CHECK(_expr)                                                         \
  do {                                                                          \
    cl_int _err = _expr;                                                        \
    if (_err == CL_SUCCESS)                                                     \
      break;                                                                    \
    fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);    \
    exit(1); /* abort(); */                                                     \
  } while (0)

#define CL_CHECK_ERR(_expr)                                                     \
  ({	                                                                           \
     cl_int _err = CL_INVALID_VALUE;                                            \
     typeof(_expr) _ret = _expr;                                                \
     if (_err != CL_SUCCESS) {                                                  \
       fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err); \
       exit(1); /* abort(); */                                                  \
     }                                                                          \
     _ret;                                                                      \
  })

// Четене на изображение от външен файл
static PPMImage *ReadPPM(const char *filename) {
     char buff[16];
     PPMImage *img;
     FILE *fp;
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

     // Проверка на метаданните
     if (buff[0] != 'P' || buff[1] != '6') {
          fprintf(stderr, "Invalid image format (must be 'P6')\n");
          exit(1);
     }

     // Заделяне на памет
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
     // Заделяне на памет за информацията във всеки пиксел
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

     fclose(fp);
     return img;
}

// Запис на изображение във външен файл
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

     // запис на размера на изображението
     fprintf(fp, "%d %d\n",img->x,img->y);

     // Запис на размерността на RGB компонента
     fprintf(fp, "%d\n", RGB_COMPONENT_COLOR);

     // Запис на данните за пикселите от изображението
     fwrite(img->data, 3 * img->x, img->y, fp);
     fclose(fp);
}

// Четене на сорс кода на kernel функцията от външен файл
char *ReadKernelProgram() {
     FILE *program_handle;
     char *program_buffer, *program_log;
     size_t program_size, log_size;
     program_handle = fopen("grayscale_kernel.cl", "r");
     if (program_handle == NULL) {
          perror("Couldn't find the program file");
          exit(1);
     }

     fseek(program_handle, 0, SEEK_END);
     program_size = ftell(program_handle);
     rewind(program_handle);
     program_buffer = (char*)malloc(program_size + 1);
     program_buffer[program_size] = '\0';
     fread(program_buffer, sizeof(char), program_size, program_handle);
     fclose(program_handle);

     return program_buffer;
}

// Информация за платформата и устройствата
void PrintSystemInfo(cl_platform_id *platforms, cl_uint platforms_n) {
     char buffer[10240];

     cl_device_id devices[100];

	printf("=== %d OpenCL platform(s) found:\n", platforms_n);
     for (int i = 0; i < platforms_n; i++) {
          cl_uint devices_n = 0;
	     CL_CHECK(clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 100, devices, &devices_n));
          clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(buffer), buffer, NULL);
          printf("  Platform ID: %d \n", i);
          printf("    CL_PLATFORM_NAME = %s\n", buffer);

          for (int j = 0; j < devices_n; j++){
               
               cl_uint buf_uint;
               cl_ulong buf_ulong;
               printf("    Device ID: %d \n", j);
               clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
               printf("      DEVICE_NAME = %s\n", buffer);
               clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL);
               printf("      DEVICE_VENDOR = %s\n", buffer);
               clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL);
               printf("      DEVICE_VERSION = %s\n", buffer);
               clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, sizeof(buffer), buffer, NULL);
               printf("      DRIVER_VERSION = %s\n", buffer);
               clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL);
               printf("      DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
               clGetDeviceInfo(devices[j], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL);
               printf("      DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
               clGetDeviceInfo(devices[j], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
               printf("      DEVICE_GLOBAL_MEM_SIZE = %llu\n", (unsigned long long)buf_ulong);
          }
          printf("\n");
     }
}

// Информация за конкретно устройствато
void PrintDeviceInfo(cl_device_id device) {
     char buffer[10240];

     cl_uint buf_uint;
     cl_ulong buf_ulong;
     printf("Selecting device for execution: \n");
     clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL);
     printf("  DEVICE_NAME = %s\n", buffer);
     clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL);
     printf("  DEVICE_VENDOR = %s\n", buffer);
     clGetDeviceInfo(device, CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL);
     printf("  DEVICE_VERSION = %s\n", buffer);
     clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(buffer), buffer, NULL);
     printf("  DRIVER_VERSION = %s\n", buffer);
     clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL);
     printf("  DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
     clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL);
     printf("  DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
     clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
     printf("  DEVICE_GLOBAL_MEM_SIZE = %llu\n", (unsigned long long)buf_ulong);
}

void ChangeColorPPM(PPMImage *img) {
     struct timeval tval_before, tval_after, tval_result;

     // Размер на изображението
     unsigned int n = img->x * img->y;
     
     // Брой на OpenCL платформите
     cl_uint platforms_n = 0;

     // Масив на наличните платформи
     cl_platform_id platforms[100];
     clGetPlatformIDs(100, platforms, &platforms_n);

     // Отпечатване на наличните платформи и устройства на екрана
     // В зависимост какво е налично в системата се модифицира clGetDeviceIDs извикването
     PrintSystemInfo(platforms, platforms_n);

     // Масив на наличните устройства в избрана платформа
     cl_device_id devices[100];

     // Контекст, чрез който ще се извърши изпълнението
     cl_context context;

     // Опашката за изпънителят (в общия случай графичен ускорител или процесор)
     cl_command_queue queue;

     // Декларации за OpenCL програма и kernel функция
     cl_program program;
     cl_kernel kernel;
 
     // Декларации за брой локални групи и размер на локална група
     size_t globalSize, localSize;
     cl_int err;
     
     // Размерност на локална група
     localSize = 64;
     
     // Брой локални групи
     globalSize = ceil(n/(float)localSize)*localSize;
     
     // Брой налични устройства
     cl_uint num_devices_returned;

     // Извличане на наличните устройства от зададената платформа в devices
     err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 100, devices, &num_devices_returned);

     // Задаване на конкретно устройсто за изпълнител
     cl_device_id target_device = devices[0];

     // Отпечатване на данните за избраното устройство за изпълнител
     PrintDeviceInfo(target_device);
     
     // Създаване на контекст от избраното устройство в devices
     context = clCreateContext(0, 1, &target_device, NULL, NULL, &err);
     
     // Създаване на опашка за изпълнение
     queue = clCreateCommandQueue(context, target_device, 0, &err);
     
     // Четене на kernel функцията от външен файл
     char *program_buffer = ReadKernelProgram();

     gettimeofday(&tval_before, NULL);
     // Създаване на програма на базата на прочетената kernel функция
     program = clCreateProgramWithSource(context, 1, (const char **) & program_buffer, NULL, &err);
     
     // Компилиране на kernel функцията според типа на устройството
     err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
     if (err) {
          char log[10240] = "";
          err = clGetProgramBuildInfo(program, target_device, CL_PROGRAM_BUILD_LOG, sizeof(log), log, NULL);
          printf("Program build log:\n%s\n", log);
     }

     // Създаване на изчисляващата kernel функция в програмата, която ще се изпълни
     kernel = clCreateKernel(program, "grayscale", &err);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld     секунди за създаване на kernel функция\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

     // Определяне на размера на паметта
     size_t bytes = n*sizeof(PPMPixel);

     // Променлива, която ще съдържа резултата от изпълнението на kernel функцията
     PPMPixel *host_result;

     // Входни данни за устройството
     cl_mem device_pixel_data;

     // Трансформирани данни от устройството
     cl_mem device_result;

     // Заделяне на памет за резултата в хоста
     host_result = (PPMPixel*)malloc(bytes);

     // Заделяне на памет за съответните входни и трансформирани данни използвани в устройството
     device_pixel_data = clCreateBuffer(context, CL_MEM_READ_ONLY, bytes, NULL, NULL);
     device_result = clCreateBuffer(context, CL_MEM_WRITE_ONLY, bytes, NULL, NULL);
     
     // Копиране на входните данни от хоста във входния масив в изчисляващото устройство
     err = clEnqueueWriteBuffer(queue, device_pixel_data, CL_TRUE, 0, bytes, img->data, 0, NULL, NULL);

     // Задаване на аргументите на kernel функцията
     err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &device_pixel_data);
     err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &device_result);
     err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &n);
     
     gettimeofday(&tval_before, NULL);

     // Изпълнение на kernel функцията върху входните данните 
     err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalSize, &localSize, 0, NULL, NULL);
     // Изчакаване за завършване на обработката на опашката
     clFinish(queue);
     
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld   секунди за изпълнението на kernel функцията\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     
     gettimeofday(&tval_before, NULL);

     // Четене на резултата от устройството
     clEnqueueReadBuffer(queue, device_result, CL_TRUE, 0, bytes, host_result, 0, NULL, NULL );

     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld   секунди за четене на резултата от устройството\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     
     gettimeofday(&tval_before, NULL);

     for (int i = 0; i < n; i++)
          img->data[i] = host_result[i];

     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld   секунди за копиране на данните обратно в масива с пикселите\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
     
     // Освобождаване на OpenCL resources
     clReleaseMemObject(device_pixel_data);
     clReleaseMemObject(device_result);
     clReleaseProgram(program);
     clReleaseKernel(kernel);
     clReleaseCommandQueue(queue);
     clReleaseContext(context);

     // Освобождаване на хост ресурси
     free(host_result);
}

int main() {
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
     printf("%ld.%06ld     секунди за обработка на данните от изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);


     gettimeofday(&tval_before, NULL);
     WritePPM("grayscale_opencl_result.ppm", image);
     gettimeofday(&tval_after, NULL);
     timersub(&tval_after, &tval_before, &tval_result);
     printf("%ld.%06ld     секунди за запис на данните в изображението\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
}

