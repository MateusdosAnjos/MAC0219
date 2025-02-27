#! /bin/bash

set -o xtrace

MEASUREMENTS=10
ITERATIONS=10
INITIAL_SIZE=16
THREADS=8

SIZE=$INITIAL_SIZE

NAMES=("mandelbrot_seq" "mandelbrot_seq_noio" "mandelbrot_pth" "mandelbrot_pth_noio" "mandelbrot_omp" "mandelbrot_omp_noio")

make
mkdir results

for NAME in ${NAMES[@]}; do
    mkdir results/$NAME

    for ((i=1; i<=$ITERATIONS; i++)); do
        if [ "$NAME" != "mandelbrot_seq" ] && [ "$NAME" != "mandelbrot_seq_noio" ]; then
            for ((j=1; j<=$THREADS; j++)); do
                perf stat -r $MEASUREMENTS ./$NAME -2.5 1.5 -2.0 2.0 $SIZE $j >> full.log 2>&1
                perf stat -r $MEASUREMENTS ./$NAME -0.8 -0.7 0.05 0.15 $SIZE $j >> seahorse.log 2>&1
                perf stat -r $MEASUREMENTS ./$NAME 0.175 0.375 -0.1 0.1 $SIZE $j >> elephant.log 2>&1
                perf stat -r $MEASUREMENTS ./$NAME -0.188 -0.012 0.554 0.754 $SIZE $j >> triple_spiral.log 2>&1
            done
            SIZE=$(($SIZE * 2))
        else
            perf stat -r $MEASUREMENTS ./$NAME -2.5 1.5 -2.0 2.0 $SIZE >> full.log 2>&1
            perf stat -r $MEASUREMENTS ./$NAME -0.8 -0.7 0.05 0.15 $SIZE >> seahorse.log 2>&1
            perf stat -r $MEASUREMENTS ./$NAME 0.175 0.375 -0.1 0.1 $SIZE >> elephant.log 2>&1
            perf stat -r $MEASUREMENTS ./$NAME -0.188 -0.012 0.554 0.754 $SIZE >> triple_spiral.log 2>&1
            SIZE=$(($SIZE * 2))
        fi
    done

    SIZE=$INITIAL_SIZE

    mv *.log results/$NAME
    rm output.ppm
done
