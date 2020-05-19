#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;
int n_threads;

int image_size;

int i_x_max;
int i_y_max;

void init(int argc, char *argv[]){
    if(argc < 7){
        printf("usage: ./mandelbrot_omp_noio c_x_min c_x_max c_y_min c_y_max image_size n_threads\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_omp_noio -2.5 1.5 -2.0 2.0 11500 10\n");
        printf("    Seahorse Valley:      ./mandelbrot_omp_noio -0.8 -0.7 0.05 0.15 11500 10\n");
        printf("    Elephant Valley:      ./mandelbrot_omp_noio 0.175 0.375 -0.1 0.1 11500 10\n");
        printf("    Triple Spiral Valley: ./mandelbrot_omp_noio -0.188 -0.012 0.554 0.754 11500 10\n");
        exit(0);
    }
    else{
        sscanf(argv[1], "%lf", &c_x_min);
        sscanf(argv[2], "%lf", &c_x_max);
        sscanf(argv[3], "%lf", &c_y_min);
        sscanf(argv[4], "%lf", &c_y_max);
        sscanf(argv[5], "%d", &image_size);
        sscanf(argv[6], "%d", &n_threads);

        i_x_max           = image_size;
        i_y_max           = image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;
    };
};

void compute_mandelbrot(){
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;

    int iteration;
    int i_x;
    int i_y;

    double c_x;
    double c_y;

    #pragma omp parallel num_threads(n_threads) \
                            private(i_x, i_y, iteration, c_x, c_y, z_x, z_y, z_x_squared, z_y_squared) \
                            shared(i_y_max, pixel_height, pixel_width, c_x_min, i_x_max, iteration_max, escape_radius_squared)
    #pragma omp for \
                    schedule(dynamic)
    for(i_y = 0; i_y < i_y_max; i_y++){
        c_y = c_y_min + i_y * pixel_height;

        if(fabs(c_y) < pixel_height / 2){
            c_y = 0.0;
        };

        for(i_x = 0; i_x < i_x_max; i_x++){
            c_x         = c_x_min + i_x * pixel_width;

            z_x         = 0.0;
            z_y         = 0.0;

            z_x_squared = 0.0;
            z_y_squared = 0.0;

            for(iteration = 0;
                iteration < iteration_max && \
                ((z_x_squared + z_y_squared) < escape_radius_squared);
                iteration++){
                z_y         = 2 * z_x * z_y + c_y;
                z_x         = z_x_squared - z_y_squared + c_x;

                z_x_squared = z_x * z_x;
                z_y_squared = z_y * z_y;
            };
        };
    };
};

int main(int argc, char *argv[]){
    init(argc, argv);

    compute_mandelbrot();

    return 0;
};
