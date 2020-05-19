#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

// thread
struct thread_data{
    int chunk_size;
    int current_index;
    pthread_mutex_t mutex_index;
    pthread_mutex_t mutex_atomic;
};

typedef struct thread_data thread_data_array;
#define THREAD_DATA_ARRAY_INITIALIZER() { 32, 0, PTHREAD_MUTEX_INITIALIZER, PTHREAD_MUTEX_INITIALIZER }

// mandelbrot
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
int image_buffer_size;

void init(int argc, char *argv[]){
    if(argc < 7){
        printf("usage: ./mandelbrot_pth_noio c_x_min c_x_max c_y_min c_y_max image_size n_threads\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_pth_noio -2.5 1.5 -2.0 2.0 11500 10\n");
        printf("    Seahorse Valley:      ./mandelbrot_pth_noio -0.8 -0.7 0.05 0.15 11500 10\n");
        printf("    Elephant Valley:      ./mandelbrot_pth_noio 0.175 0.375 -0.1 0.1 11500 10\n");
        printf("    Triple Spiral Valley: ./mandelbrot_pth_noio -0.188 -0.012 0.554 0.754 11500 10\n");
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
        image_buffer_size = image_size * image_size;

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

void *compute_mandelbrot_thread(void *args){
    // thread data
    thread_data_array *data = (thread_data_array *) args;
    int from_index;
    int to_index;
    int chunk_start_x;
    int chunk_start_y;

    // mandelbrot data
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

    while(1)
    {
       // atomic chunk retrieval
        pthread_mutex_lock(&(data->mutex_index));
        from_index = data->current_index;
        to_index = from_index + data->chunk_size;
        to_index = to_index > image_buffer_size ? image_buffer_size : to_index;
        data->current_index = to_index;
        pthread_mutex_unlock(&(data->mutex_index));

        // job done
        if(from_index >= image_buffer_size)
            break;

        chunk_start_y = from_index/i_x_max;
        chunk_start_x = from_index%i_x_max;

        for(i_y = chunk_start_y; i_y < i_y_max; i_y++){
            c_y = c_y_min + i_y * pixel_height;

            if(fabs(c_y) < pixel_height / 2){
                c_y = 0.0;
            };

            int start_x = i_y == chunk_start_y ? chunk_start_x : 0;
            for(i_x = start_x; i_x < i_x_max; i_x++){
                // chunk end?
                if((i_y_max * i_y) + i_x >= to_index)
                    break;

                c_x         = c_x_min + i_x * pixel_width;

                z_x         = 0.0;
                z_y         = 0.0;

                z_x_squared = 0.0;
                z_y_squared = 0.0;

                for(iteration = 0;
                    iteration < iteration_max &&                        \
                    ((z_x_squared + z_y_squared) < escape_radius_squared);
                    iteration++){
                    z_y         = 2 * z_x * z_y + c_y;
                    z_x         = z_x_squared - z_y_squared + c_x;

                    z_x_squared = z_x * z_x;
                    z_y_squared = z_y * z_y;
                };
            };
            // chunk end?
            if((i_y_max * i_y) + i_x >= to_index)
                break;
        };
    }
    pthread_exit(NULL);
}

int main(int argc, char *argv[]){
    init(argc, argv);

    if(n_threads == 1){
        compute_mandelbrot();
    } else {
        thread_data_array data_array = THREAD_DATA_ARRAY_INITIALIZER();

        int i;
        pthread_t thread_pool[n_threads];
        for (i = 0; i < sizeof(thread_pool) / sizeof(thread_pool[0]); ++i)
        {
            pthread_create(&thread_pool[i], NULL, compute_mandelbrot_thread, &data_array);
        }
        for (i = 0; i < sizeof(thread_pool) / sizeof(thread_pool[0]); ++i)
        {
            pthread_join(thread_pool[i], NULL);
        }
    }

    return 0;
};
