OUTPUT=mandelbrot

IMAGE=.ppm

CC=gcc
CC_OPT=-std=c11
CC_OMPI=mpicc

CC_OMP=-fopenmp
CC_PTH=-pthread

CUDA=nvcc

.PHONY: all
all: $(OUTPUT)_cuda_noio $(OUTPUT)_ompi_noio $(OUTPUT)_omp_noio $(OUTPUT)_pth_noio $(OUTPUT)_seq_noio

$(OUTPUT)_cuda: $(OUTPUT)_cuda.cu
	$(CUDA) -o $(OUTPUT)_cuda $(OUTPUT)_cuda.cu

$(OUTPUT)_cuda_noio: $(OUTPUT)_cuda_noio.cu
	$(CUDA) -o $(OUTPUT)_cuda_noio $(OUTPUT)_cuda_noio.cu

$(OUTPUT)_ompi: $(OUTPUT)_ompi.c
	$(CC_OMPI) -o $(OUTPUT)_ompi $(OUTPUT)_ompi.c

$(OUTPUT)_ompi_noio: $(OUTPUT)_ompi_noio.c
	$(CC_OMPI) -o $(OUTPUT)_ompi_noio $(OUTPUT)_ompi_noio.c

$(OUTPUT)_omp: $(OUTPUT)_omp.c
	$(CC) -o $(OUTPUT)_omp $(CC_OPT) $(CC_OMP) $(OUTPUT)_omp.c

$(OUTPUT)_omp_noio: $(OUTPUT)_omp_noio.c
	$(CC) -o $(OUTPUT)_omp_noio $(CC_OPT) $(CC_OMP) $(OUTPUT)_omp_noio.c

$(OUTPUT)_pth: $(OUTPUT)_pth.c
	$(CC) -o $(OUTPUT)_pth $(CC_OPT) $(CC_PTH) $(OUTPUT)_pth.c

$(OUTPUT)_pth_noio: $(OUTPUT)_pth_noio.c
	$(CC) -o $(OUTPUT)_pth_noio $(CC_OPT) $(CC_PTH) $(OUTPUT)_pth_noio.c

$(OUTPUT)_seq: $(OUTPUT)_seq.c
	$(CC) -o $(OUTPUT)_seq $(CC_OPT) $(OUTPUT)_seq.c

$(OUTPUT)_seq_noio: $(OUTPUT)_seq_noio.c
	$(CC) -o $(OUTPUT)_seq_noio $(CC_OPT) $(OUTPUT)_seq_noio.c

.PHONY: clean
clean:
	rm $(OUTPUT)_cuda $(OUTPUT)_cuda_noio $(OUTPUT)_omp $(OUTPUT)_ompi $(OUTPUT)_ompi_noio $(OUTPUT)_omp_noio $(OUTPUT)_pth $(OUTPUT)_pth_noio $(OUTPUT)_seq $(OUTPUT)_seq_noio *$(IMAGE)
