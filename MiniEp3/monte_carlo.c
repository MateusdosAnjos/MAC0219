#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif

#define FUNCTIONS 1

struct timer_info {
    clock_t c_start;
    clock_t c_end;
    struct timespec t_start;
    struct timespec t_end;
    struct timeval v_start;
    struct timeval v_end;
};

struct timer_info timer;

char *usage_message = "usage: ./monte_carlo SAMPLES FUNCTION_ID N_THREADS\n";

struct function {
    long double (*f)(long double);
    long double interval[2];
};

long double rand_interval[] = {0.0, (long double) RAND_MAX};

long double f1(long double x){
    return 2 / (sqrt(1 - (x * x)));
}

struct function functions[] = {
                               {&f1, {0.0, 1.0}}
};

// Your thread data structures go here

struct thread_data{
    long double (*f)(long double);
    long double *samples;
    int size;
    int chunk_size;
    int current_index;
    long double sum;
    pthread_mutex_t mutex_index;
    pthread_mutex_t mutex_sum;
};

typedef struct thread_data thread_data_array;
#define THREAD_DATA_ARRAY_INITIALIZER(func, samples, size) { func, samples, size, 1000, 0, 0, PTHREAD_MUTEX_INITIALIZER, PTHREAD_MUTEX_INITIALIZER }

// End of data structures

long double *samples;
long double *results;

long double map_intervals(long double x, long double *interval_from, long double *interval_to){
    x -= interval_from[0];
    x /= (interval_from[1] - interval_from[0]);
    x *= (interval_to[1] - interval_to[0]);
    x += interval_to[0];
    return x;
}

long double *uniform_sample(long double *interval, long double *samples, int size){
    for(int i = 0; i < size; i++){
        samples[i] = map_intervals((long double) rand(),
                                   rand_interval,
                                   interval);
    }
    return samples;
}

void print_array(long double *sample, int size){
    printf("array of size [%d]: [", size);

    for(int i = 0; i < size; i++){
        printf("%Lf", sample[i]);

        if(i != size - 1){
            printf(", ");
        }
    }

    printf("]\n");
}

long double monte_carlo_integrate(long double (*f)(long double), long double *samples, int size){
    // Your sequential code goes here
    long double accum = 0.0;
    for(int i = 0 ; i < size; i++)
    {
        accum += (*f)(samples[i]);
    }
    return accum / (long double) size;
}

void *monte_carlo_integrate_thread(void *args){
    // Your pthreads code goes here
    thread_data_array *data = (thread_data_array *) args;
    int from_index;
    int to_index;
    long double accum_chunk;
    while(1)
    {
        // atomic sample retrieve
        pthread_mutex_lock(&(data->mutex_index));
        from_index = data->current_index;
        to_index = from_index + data->chunk_size;
        to_index = to_index > data->size ? data->size : to_index;
        data->current_index = to_index;
        pthread_mutex_unlock(&(data->mutex_index));
        
        // job done
        if(from_index >= data->size)
            break;
        
        // calculate chunk result
        accum_chunk = 0.0;
        for(int i = from_index; i < to_index; i++)
            accum_chunk += (*data->f)(data->samples[i]);
  
        // atomic sum aggregation
        pthread_mutex_lock(&(data->mutex_sum));
        data->sum += accum_chunk;
        pthread_mutex_unlock(&(data->mutex_sum));
    }
    pthread_exit(NULL);
}

int main(int argc, char **argv){
    if(argc != 4){
        printf("%s", usage_message); // erro
        exit(-1);
    } else if(atoi(argv[2]) >= FUNCTIONS || atoi(argv[2]) < 0){
        printf("Error: FUNCTION_ID must in [0,%d]\n", FUNCTIONS - 1);
        printf("%s", usage_message); // erro
        exit(-1);
    } else if(atoi(argv[3]) < 0){
        printf("Error: I need at least 1 thread\n");
        printf("%s", usage_message); // erro
        exit(-1);
    }

    if(DEBUG){
        printf("Running on: [debug mode]\n");
        printf("Samples: [%s]\n", argv[1]);
        printf("Function id: [%s]\n", argv[2]);
        printf("Threads: [%s]\n", argv[3]);
        printf("Array size on memory: [%.2LFGB]\n", ((long double) atoi(argv[1]) * sizeof(long double)) / 1000000000.0);
    }

    srand(time(NULL));

    int size = atoi(argv[1]);
    struct function target_function = functions[atoi(argv[2])];
    int n_threads = atoi(argv[3]);

    samples = malloc(size * sizeof(long double));

    long double estimate;

    if(n_threads == 1){
        if(DEBUG){
            printf("Running sequential version\n");
        }

        timer.c_start = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
        gettimeofday(&timer.v_start, NULL);
        estimate = monte_carlo_integrate(target_function.f,
                                         uniform_sample(target_function.interval,
                                                        samples,
                                                        size),
                                         size);
        
        timer.c_end = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
        gettimeofday(&timer.v_end, NULL);
    } else {
        if(DEBUG){
            printf("Running parallel version\n");
        }

        timer.c_start = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
        gettimeofday(&timer.v_start, NULL);

        // Your pthreads code goes here
        thread_data_array data_array = THREAD_DATA_ARRAY_INITIALIZER(target_function.f,
                                                                     uniform_sample(target_function.interval, 
                                                                                    samples, 
                                                                                    size),
                                                                     size);
        int i;
        pthread_t thread_pool[n_threads];
        for (i = 0; i < sizeof(thread_pool) / sizeof(thread_pool[0]); ++i)
        {
            pthread_create(&thread_pool[i], NULL, monte_carlo_integrate_thread, &data_array);
        }
        for (i = 0; i < sizeof(thread_pool) / sizeof(thread_pool[0]); ++i)
        {
            pthread_join(thread_pool[i], NULL);
        }
        estimate = data_array.sum / (long double) data_array.size;
        // Your pthreads code ends here

        timer.c_end = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
        gettimeofday(&timer.v_end, NULL);

        if(DEBUG && VERBOSE){
            print_array(results, n_threads);
        }
    }

    if(DEBUG){
        if(VERBOSE){
            print_array(samples, size);
            printf("Estimate: [%.33LF]\n", estimate);
        }
        printf("%.16LF, [%f, clock], [%f, clock_gettime], [%f, gettimeofday]\n",
               estimate,
               (double) (timer.c_end - timer.c_start) / (double) CLOCKS_PER_SEC,
               (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
               (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0,
               (double) (timer.v_end.tv_sec - timer.v_start.tv_sec) +
               (double) (timer.v_end.tv_usec - timer.v_start.tv_usec) / 1000000.0);
    } else {
        printf("%.16LF, %f\n",
               estimate,
               (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
               (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
    }
    return 0;
}
