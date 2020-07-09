/******************************************************************************
* FILE: mpi_hello.c
* DESCRIPTION:
*   MPI tutorial example code: Simple hello world program
* AUTHOR: Blaise Barney
* LAST REVISED: 03/05/10
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#define  MASTER		0

int update(int myoffset, int chunk, int myid, int *data) {
    int i;
    int mysum;
    /* Perform addition to each of my array elements and keep my sum */
    mysum = 0;
    for(i=myoffset; i < myoffset + chunk; i++) {
        mysum += data[i];
    }
    printf("Task %d mysum = %d\n", myid, mysum);
    return(mysum);
}

int main (int argc, char *argv[])
{
    int numtasks, taskid, len;
    int n = 18;
    int i, j;
    int data[n], sum[4];
    int offset;
    int mysum = 0;
    int verificador = 0;
    int chunksize;
    int leftover;
    int rc, dest, tag1, tag2, source;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    chunksize = (n / numtasks);
    leftover = (n % numtasks);

    MPI_Get_processor_name(hostname, &len);


    // printf("Hello from task %d!\n", taskid);
    if (taskid == MASTER) {
        printf("MASTER: Number of MPI tasks is: %d\n", numtasks);
        
        /* Initialize the array */
        for(i = 0; i < n; i++) {
            data[i] = i;
            verificador += i;
        }
        printf("Initialized array sum = %d\n", verificador);
        printf("numtasks= %d  chunksize= %d  leftover= %d\n", numtasks, chunksize, leftover);

        /* Send each task its portion of the array - master keeps 1st part plus leftover elements */
        offset = chunksize + leftover;
        for (dest = 1; dest < numtasks; dest++) {
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            MPI_Send(&data[offset], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
            printf("Sent %d elements to task %d offset= %d\n", chunksize, dest, offset);
            offset = offset + chunksize;
        }

        /* Master does its part of the work */
        offset = 0;
        sum[0] = update(offset, (chunksize + leftover), taskid, data);

        /* Wait to receive results from each task */
        for (i = 1; i < numtasks; i++) {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            MPI_Recv(&sum[i], chunksize, MPI_DOUBLE, source, tag2,
                     MPI_COMM_WORLD, &status);
        }

        /* Get final sum and print sample results */
        //MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

        printf("Sample results: \n");
        offset = 0;
        for (i = 0; i < numtasks; i++) {
            printf("%d\n", sum[i]);
        }

    }  /* end of master section */
    if (taskid > MASTER) {

        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&data[offset], chunksize, MPI_DOUBLE, source, tag2,
                 MPI_COMM_WORLD, &status);

        /* Do my part of the work */
        mysum = update(offset, chunksize, taskid, data);

        /* Send my results back to the master task */
        dest = MASTER;
        MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        MPI_Send(&mysum, chunksize, MPI_DOUBLE, MASTER, tag2, MPI_COMM_WORLD);

        /* Use sum reduction operation to obtain final sum */
        // MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    } /* end of non-master */

    MPI_Finalize();
}