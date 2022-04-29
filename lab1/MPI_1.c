#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

// 1800 -->  ~30 sec
#define N 1800
#define epsilon 10e-5

void MatrixMultVector(double *A, double *vector, double *result, int height, int width)
{
    for (int i = 0; i < height; i++)
    {
        // result[i] = 0;
        for (int j = 0; j < width; j++)
        {
            result[i] += A[i * width + j] * vector[j];
        }
    }
}

void parallelMatrixMultVector(double *A_part, double *vector, double *result,
                              int A_p_width, int A_p_height, int procSize, int procRank)
{
    // A_part is N * rowsNum[procRank]
    // res is vector (one column matrix) of A_part height
    double *tmp = (double *)calloc(A_p_height, sizeof(double));
    MatrixMultVector(A_part, vector, tmp, A_p_height, A_p_width);

    // array of numbers of elemets send by each process
    int *recvcounts = (int *)calloc(procSize, sizeof(int));
    for (int i = 0; i < procSize; ++i)
    {
        recvcounts[i] = A_p_width / procSize;
        recvcounts[i] += (i < A_p_width % procSize);
    }

    int *displs = (int *)calloc(procSize, sizeof(int));
    for (int i = 1; i < procSize; ++i)
    {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    MPI_Allgatherv(tmp, A_p_height, MPI_DOUBLE, result, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

    free(tmp);
    free(displs);
    free(recvcounts);
}

double sqrNormOfVector(double *vector)
{
    double result = 0;

    for (int i = 0; i < N; i++)
    {
        result += vector[i] * vector[i];
    }

    return result;
}

double *SubVectors(double *vect1, double *vect2)
{
    double *result = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        result[i] = vect1[i] - vect2[i];
    }
    return result;
}

double *SumVectors(double *vect1, double *vect2)
{
    double *result = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        result[i] = vect1[i] + vect2[i];
    }
    return result;
}

double VectorsMult(double *vect1, double *vect2)
{
    double result = 0;
    for (int i = 0; i < N; i++)
    {
        result += vect1[i] * vect2[i];
    }
    return result;
}

double *ScalarVectMult(double *vector, double scalar)
{
    double *result = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        result[i] = vector[i] * scalar;
    }

    return result;
}

int main(int argc, char **argv)
{
    int procSize;
    int procRank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);

    // Matrix A initializing and sending parts to other processes
    // Matrix A generated on 0 process
    double *A_full = NULL;

    if (procRank == 0)
    {
        A_full = (double *)malloc(N * N * sizeof(double));
        for (int i = 0; i < N; i++)
        {
            for (int j = i; j < N; j++)
            {
                if (i == j)
                {
                    A_full[i * N + j] = (double)(rand() % 100) + 200;
                    // A_full[i * N + j] = 2.0;
                }
                else
                {
                    A_full[i * N + j] = A_full[j * N + i] = (double)(rand() % 100);
                    // A_full[i * N + j] = A_full[j * N + i] = 1.0;
                }
            }
        }
    }

    double start = MPI_Wtime();

    int count1 = N / procSize;
    int count2 = N % procSize;

    int *rowsNum = (int *)calloc(procSize, sizeof(int));
    // number of elements in each process part
    int *sendcounts = (int *)calloc(procSize, sizeof(int));

    for (int i = 0; i < procSize; i++)
    {
        rowsNum[i] = count1 + (i < count2);
        sendcounts[i] = rowsNum[i] * N;
    }

    // i-ое значение в массиве displs определяет смещение относительно
    // начала sendbuf для данных, посылаемых процессу i
    int *displs = (int *)calloc(procSize, sizeof(int));
    displs[0] = 0;
    for (int i = 1; i < procSize; i++)
    {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    double *A_part = (double *)malloc(sizeof(double) * sendcounts[procRank]);
    MPI_Scatterv(A_full, sendcounts, displs, MPI_DOUBLE, A_part, sendcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // A_part dimension
    int A_p_width = N;
    int A_p_height = rowsNum[procRank];

    // Vector b generating and sending
    // b is common for all processes
    double *b = (double *)malloc(N * sizeof(double));
    if (procRank == 0)
    {
        for (int i = 0; i < N; i++)
        {
            b[i] = (double)(rand() % 100);
            // b[i] = N + 1;
        }
    }

    // sending vector b to every process
    MPI_Bcast(b, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // b dimension
    int b_width = 1;
    int b_height = N;

    double *X = (double *)calloc(N, sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *prev_r = (double *)malloc(N * sizeof(double));
    double *z = (double *)malloc(N * sizeof(double));

    double square_b_len = 0; // square of vector b length
    square_b_len = VectorsMult(b, b);
    for (int i = 0; i < N; i++)
    {
        // square_b_len += b[i] * b[i];
        r[i] = b[i]; // r0 = b - Ax0, x0 - нулеой вектор
        z[i] = r[i]; // z0 = r0
    }

    // Print matrix A, b, X
    //==========================================
    /*for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f ", A_full[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
    for (int i = 0; i < N; i++)
    {
        printf("%f ", b[i]);
    }
    printf("\n");
    printf("\n");

    for (int i = 0; i < N; i++)
    {
        printf("%f ", X[i]);
    }
    printf("\n");*/
    //==========================================

    double *temp = (double *)calloc(N, sizeof(double));
    double alpha, beta;
    int cntr = 0;
    int cnt = 0;
    double restrictor = epsilon * epsilon * square_b_len;

    do
    {
        // A*z(n)
        // MatrixMultVector(A_part, z, temp, A_p_height, A_p_width);
        parallelMatrixMultVector(A_part, z, temp, A_p_width, A_p_height, procSize, procRank);

        // alpha(n+1) = (r(n), r(n)) / (Az(n), z(n))
        double oldsqrR = VectorsMult(r, r);
        alpha = oldsqrR / VectorsMult(temp, z);

        // X(n+1) = X(n) + alpha(n+1)*z(n)
        X = SumVectors(X, ScalarVectMult(z, alpha));

        prev_r = r;
        // r(n+1) = r(n) - alpha(n+1)Az(n)
        r = SubVectors(r, ScalarVectMult(temp, alpha));

        // beta(n+1) = (r(n+1), r(n+1)) / (r(n), r(n))
        double temp_res = VectorsMult(r, r);
        beta = temp_res / oldsqrR;

        // z(n+1) = r(n+1) + beta(n+1)z(n)
        z = SumVectors(r, ScalarVectMult(z, beta));

        if (temp_res < restrictor)
        {
            cntr++;
        }

        cnt++;
    } while (cntr < 3);

    double finish = MPI_Wtime();
    double total_time = finish - start;

    if (procRank == 0)
    {
        printf("Iterations count: %d\n", cnt);
        printf("Total time: %f \n", total_time);
    }

    // Print result
    /*if (procRank == 0)
    {
        printf("\n");
        for (int i = 0; i < N; i++)
        {
            printf("%f ", X[i]);
        }
        printf("\n");
    }*/

    free(A_full);
    free(b);
    free(r);
    free(prev_r);
    free(z);
    free(X);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
