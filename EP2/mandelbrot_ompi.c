#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define MASTER 0

int numtasks, taskid, len, dest, offset, i, j, tag1,
    tag2, source, chunksize, leftover;
char hostname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;


double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
unsigned char **image_buffer;

int i_x_max;
int i_y_max;
int image_buffer_size;

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
    int rgb_size = 3;
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);
    printf("eu aloquei image buffer no id %d", taskid);
    for(int i = 0; i < image_buffer_size; i++){
        image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
    };
};

void init(int argc, char *argv[]){
    if(argc < 7){
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    }
    else{
        sscanf(argv[2], "%lf", &c_x_min);
        sscanf(argv[3], "%lf", &c_x_max);
        sscanf(argv[4], "%lf", &c_y_min);
        sscanf(argv[5], "%lf", &c_y_max);
        sscanf(argv[6], "%d", &image_size);

        i_x_max           = image_size;
        i_y_max           = image_size;
        
        image_buffer_size = image_size * image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;
    };
};

void update_rgb_buffer(int iteration, int x, int y){
    int color;
    //printf("entrei na update rgb buffer com o taskid %d\n", taskid);
    if(iteration == iteration_max){
        image_buffer[(i_y_max * y) + x][0] = colors[gradient_size][0];
        image_buffer[(i_y_max * y) + x][1] = colors[gradient_size][1];
        image_buffer[(i_y_max * y) + x][2] = colors[gradient_size][2];
        // printf("%d %d %d", colors[gradient_size][0], colors[gradient_size][1], colors[gradient_size][2]);
    }
    else{
        color = iteration % gradient_size;

        image_buffer[(i_y_max * y) + x][0] = colors[color][0];
        image_buffer[(i_y_max * y) + x][1] = colors[color][1];
        image_buffer[(i_y_max * y) + x][2] = colors[color][2];
        // printf("%d %d %d", colors[color][0], colors[color][1], colors[color][2]);

    };
};

void write_to_file(){
    FILE * file;
    char * filename               = "output.ppm";
    char * comment                = "# ";

    int max_color_component_value = 255;

    file = fopen(filename,"wb");

    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            i_x_max, i_y_max, max_color_component_value);

    for(int i = 0; i < image_buffer_size; i++){
        fwrite(image_buffer[i], 1 , 3, file);
    };

    fclose(file);
};

void update (int myoffset, int chunk, int myid) {
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

    printf("ENTREI EM update ID: %d\n", myid);
    
    for(i_y = myoffset; i_y < myoffset + chunk; i_y++){
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

            update_rgb_buffer(iteration, i_x, i_y);
        };
    };

}

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

    chunksize = (i_y_max / numtasks);
    leftover = (i_y_max % numtasks);
    tag1 = 1;
    tag2 = 2;

    printf("ENTREI EM compute_mandelbrot ID: %d\n", taskid);

    if (taskid == MASTER) {
        printf("numtasks= %d  chunksize= %d  leftover= %d\n", numtasks, chunksize, leftover);

        /* Send each task its portion of the array - master keeps 1st part plus leftover elements */
        offset = chunksize + leftover;
        for (dest = 1; dest < numtasks; dest++) {
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            MPI_Send(&image_buffer[offset], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
            printf("Sent %d elements to task %d offset= %d\n", chunksize, dest,offset);
            offset = offset + chunksize;
        }


        /* Master does its part of the work */
        offset = 0;
        update (offset,  chunksize + leftover, taskid);

        /* Wait to receive results from each task */
        for (i=1; i < numtasks; i++) {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            printf("CHEGAMOS AQUI PADRIN!\n");
            MPI_Recv(&image_buffer[offset], chunksize, MPI_DOUBLE, source, tag2,
                     MPI_COMM_WORLD, &status);
            printf("EXECUTEI UMA VEZ");
        }

     }
    if (taskid > MASTER) {
        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&image_buffer[offset], chunksize, MPI_DOUBLE, source, tag2,
                 MPI_COMM_WORLD, &status);
        printf("Received %d elements to task %d offset= %d\n", chunksize, taskid, offset);
        
        /* Do my part of the work */
        update(offset, chunksize, taskid);
        
        /* Send my results back to the master task */
        dest = MASTER;
        MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        MPI_Send(&image_buffer[offset], chunksize, MPI_DOUBLE, MASTER, tag2, MPI_COMM_WORLD);

        /* Use sum reduction operation to obtain final sum */
        //MPI_Reduce(&image_buffer, offset, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    } /* end of non-master */

};

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    init(argc, argv);
    allocate_image_buffer();
    
    compute_mandelbrot();

    MPI_Finalize();

    //write_to_file();

    return 0;
};
