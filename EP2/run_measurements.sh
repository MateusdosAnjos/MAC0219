#! /bin/bash

set -o xtrace

MEASUREMENTS=15
ITERATIONS=8
INITIAL_SIZE=4096
THREADS=32

SIZE=$INITIAL_SIZE

NAMES=("mandelbrot_omp_noio")
# "mandelbrot_seq_noio" "mandelbrot_omp_noio" "mandelbrot_pth_noio"
make
mkdir results

for NAME in ${NAMES[@]}; do
    mkdir results/$NAME

    for ((i=1; i<=$ITERATIONS; i++)); do
        if [ "$NAME" != "mandelbrot_seq_noio" ]; then
            for ((j=1; j<=$THREADS; j*=2)); do
                perf stat -r $MEASUREMENTS ./$NAME -0.188 -0.012 0.554 0.754 $SIZE $j >> triple_spiral.log 2>&1
            done
        else
            perf stat -r $MEASUREMENTS ./$NAME -0.188 -0.012 0.554 0.754 $SIZE >> triple_spiral.log 2>&1
        fi
    done

    mv *.log results/$NAME
    rm output.ppm
done
