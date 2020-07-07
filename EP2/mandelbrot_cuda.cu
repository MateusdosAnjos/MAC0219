#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"


double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
unsigned char *image_buffer_host;
unsigned char *image_buffer_device;

int i_x_max;
int i_y_max;
int image_buffer_size;

int rgb_size = 3;

int gradient_size = 16;
int colors[17][3] = {
                        {66, 30, 15},
                        {25, 7, 26},
                        {9, 1, 47},
                        {4, 4, 73},
                        {0, 7, 100},
                        {12, 44, 138},
                        {24, 82, 177},
                        {57, 125, 209},
                        {134, 181, 229},
                        {211, 236, 248},
                        {241, 233, 191},
                        {248, 201, 95},
                        {255, 170, 0},
                        {204, 128, 0},
                        {153, 87, 0},
                        {106, 52, 3},
                        {16, 16, 16},
                    };



void allocate_image_buffer(){
    image_buffer_host = (unsigned char *) malloc(sizeof(unsigned char) * image_buffer_size * rgb_size);

    // for(int i = 0; i < image_buffer_size; i++){
    //     image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
    // };
};

void init(int argc, char *argv[]){
    if(argc < 6){
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    }
    else{
        sscanf(argv[1], "%lf", &c_x_min);
        sscanf(argv[2], "%lf", &c_x_max);
        sscanf(argv[3], "%lf", &c_y_min);
        sscanf(argv[4], "%lf", &c_y_max);
        sscanf(argv[5], "%d", &image_size);

        i_x_max           = image_size;
        i_y_max           = image_size;
        image_buffer_size = image_size * image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;
    };
};

// void update_rgb_buffer(int iteration, int x, int y){
//     int color;

//     if(iteration == iteration_max){
//         image_buffer[(i_y_max * y) + x][0] = colors[gradient_size][0];
//         image_buffer[(i_y_max * y) + x][1] = colors[gradient_size][1];
//         image_buffer[(i_y_max * y) + x][2] = colors[gradient_size][2];
//     }
//     else{
//         color = iteration % gradient_size;

//         image_buffer[(i_y_max * y) + x][0] = colors[color][0];
//         image_buffer[(i_y_max * y) + x][1] = colors[color][1];
//         image_buffer[(i_y_max * y) + x][2] = colors[color][2];
//     };
// };

void write_to_file(){
    FILE * file;
    const char * filename               = "output.ppm";
    const char * comment                = "# ";

    int max_color_component_value = 255;

    file = fopen(filename,"wb");

    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            i_x_max, i_y_max, max_color_component_value);

    for(int i = 0; i < image_buffer_size * rgb_size; i++){
        fwrite(image_buffer_host + i, 1, 1, file);
    };

    fclose(file);
};


__global__ void compute_mandelbrot_gpu(double pixel_height, double pixel_width, double c_x_min, double c_y_min, \
                                       int image_size, int iteration_max, unsigned char* image_buffer_device){
    int i_x = threadIdx.x + blockDim.x * blockIdx.x;
    int i_y = threadIdx.y + blockDim.y * blockIdx.y;
    
    // printf("i_x=%d | i_y=%d", i_x, i_y);

    // declaração variáveis para a função update_rgb_buffer
    int color;
    int rgb_size = 3;
    int gradient_size = 16;
    int colors[17][3] = {
                        {66, 30, 15},
                        {25, 7, 26},
                        {9, 1, 47},
                        {4, 4, 73},
                        {0, 7, 100},
                        {12, 44, 138},
                        {24, 82, 177},
                        {57, 125, 209},
                        {134, 181, 229},
                        {211, 236, 248},
                        {241, 233, 191},
                        {248, 201, 95},
                        {255, 170, 0},
                        {204, 128, 0},
                        {153, 87, 0},
                        {106, 52, 3},
                        {16, 16, 16},
                    };
    
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;

    int iteration;
    // int i_x;

    double c_x;
    double c_y;

    // int i_x_max = image_size;
    int i_y_max = image_size;

    if(i_x < image_size && i_y < image_size){

        c_y = c_y_min + i_y * pixel_height;

        if(fabs(c_y) < pixel_height / 2){
            c_y = 0.0;
        };
        
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

            
            if(iteration == iteration_max){
                image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 0] = colors[gradient_size][0];
                image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 1] = colors[gradient_size][1];
                image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 2] = colors[gradient_size][2];
            }
            else{
                color = iteration % gradient_size;

                image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 0] = colors[color][0];
                image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 1] = colors[color][1];
                image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 2] = colors[color][2];
            };

            // printf("color= %d | i_y= %d | i_x= %d | i_y_max= %d | iteration= %d | image_buffer_device(0,1,2)= (%u, %u, %u)\n", color, i_y, i_x, i_y_max, iteration, \
            //         image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 0], image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 1], \
            //         image_buffer_device[((i_y_max * i_y) + i_x) * rgb_size + 2]);
    };
}


int main(int argc, char *argv[]){
    init(argc, argv);
    
    allocate_image_buffer();
       
    int dimBlock, dimGrid;
    
    // define estrategia do grid e block: quanto mais proximo de 32 a dimensao do block melhor
    // devido ao warp size
    if(image_size > 32){
        dimBlock = 32;
        dimGrid = (int) (image_size / dimBlock) + 1;
    }
    else{
        dimBlock = image_size;
        dimGrid = 1;
    };

    // printf("dimBlock = %d | dimGrid = %d\n", dimBlock, dimGrid);

    // alocando espaço no device
    cudaMalloc((void **)&image_buffer_device, sizeof(unsigned char) * image_buffer_size * rgb_size);
    
    // transferir dados do device para o host com cudaMemcpy
    cudaMemcpy(image_buffer_device, image_buffer_host, sizeof(unsigned char) * image_buffer_size * rgb_size, cudaMemcpyHostToDevice);
    
    // dimensionamento do grid e do block
    dim3 block(dimBlock, dimBlock);
    dim3 grid(dimGrid, dimGrid);

    // chama função para executar no device (GPU)
    compute_mandelbrot_gpu<<<grid, block>>>(pixel_height, pixel_width, c_x_min, c_y_min, image_size, iteration_max, image_buffer_device);
    cudaDeviceSynchronize();

    // passando os dados do array image_buffer do device para o host
    cudaMemcpy(image_buffer_host, image_buffer_device, sizeof(unsigned char) * image_buffer_size * rgb_size, cudaMemcpyDeviceToHost);

    cudaFree(image_buffer_device);

    cudaDeviceReset();

    write_to_file();

    return 0;
};
