#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

// 1800 --> ~28 sec

#define N 1800
#define epsilon 10e-5

// OLD Without using matrix transpose
void MatrixMultVector(double *A, double *vector, double *result, int height, int width)
{
    for (int j = 0; j < width; j++)
    {
        // result[i] = 0;
        for (int i = 0; i < height; i++)
        {
            result[i] += A[j * height + i] * vector[j];
        }
    }
}

// with using matrix A trancpose
/*void MatrixMultVector(double* A, double* vector, double* result, int height, int width)
{
    for (int j = 0; j < width; j++)
    {
        //result[i] = 0;
        for (int i = 0; i < height; i++)
        {
            double Ap = A[i * width + j];
            double vj = vector[j];
            double mul = Ap * vj;
            double prev_res = result[i];
            result[i] += mul;
        }
    }
}*/

void parallelMatrixMultVector(double *A_part, double *vector, double *result,
                              int A_p_width, int A_p_height, int procSize, int procRank)
{
    // A_part is N * rowsNum[procRank]
    // res is vector (one column matrix) of A_part height
    double *tmp = (double *)calloc(A_p_height, sizeof(double));
    MatrixMultVector(A_part, vector, tmp, A_p_height, A_p_width);

    double *full_result = NULL;
    if (procRank == 0)
    {
        full_result = (double *)calloc(A_p_height, sizeof(double));
    }

    MPI_Reduce(tmp, full_result, A_p_height, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // array of numbers of elemets send by each process
    int *sendcounts = (int *)calloc(procSize, sizeof(int));
    for (int i = 0; i < procSize; ++i)
    {
        sendcounts[i] = A_p_height / procSize;
        sendcounts[i] += (i < A_p_height % procSize);
    }

    int *displs = (int *)calloc(procSize, sizeof(int));
    displs[0] = 0;
    for (int i = 1; i < procSize; ++i)
    {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    MPI_Scatterv(full_result, sendcounts, displs, MPI_DOUBLE, result, sendcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (procRank == 0)
    {
        free(full_result);
    }
}

double *SubVectors(double *vect1, double *vect2, int size)
{
    double *result = (double *)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
    {
        result[i] = vect1[i] - vect2[i];
    }
    return result;
}

double *SumVectors(double *vect1, double *vect2, int size)
{
    double *result = (double *)malloc(size * sizeof(double));
    for (int i = 0; i < size; i++)
    {
        result[i] = vect1[i] + vect2[i];
    }
    return result;
}

double *ScalarVectMult(double *vector, double scalar, int size)
{
    double *result = (double *)malloc(size * sizeof(double));

    for (int i = 0; i < size; i++)
    {
        result[i] = vector[i] * scalar;
    }

    return result;
}

void VectorParallelMult(double *vect1, double *vect2, double *res, int size)
{
    double tmp = 0;
    for (int i = 0; i < size; i++)
    {
        tmp += vect1[i] * vect2[i];
    }

    MPI_Allreduce(&tmp, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int procSize;
    int procRank;
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

    // number of elements in each process part
    int *sendcounts = (int *)calloc(procSize, sizeof(int));

    for (int i = 0; i < procSize; i++)
    {
        // rowsNum[i] = N / procSize + (i < (N % procSize));
        // sendcounts[i] = rowsNum[i] * N;

        sendcounts[i] = N / procSize + (i < N % procSize);
        sendcounts[i] *= N;
    }

    // i-ое значение в массиве displs определяет смещение относительно
    // начала sendbuf для данных, посылаемых процессу i
    int *displs = (int *)calloc(procSize, sizeof(int));
    displs[0] = 0;
    for (int i = 1; i < procSize; i++)
    {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    double *temp_part = (double *)malloc(sizeof(double) * sendcounts[procRank]);
    MPI_Scatterv(A_full, sendcounts, displs, MPI_DOUBLE, temp_part, sendcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // A_part dimension
    // int A_p_width = rowsNum[procRank];
    int A_p_width = sendcounts[procRank] / N;
    int A_p_height = N;

    double *A_part = (double *)malloc(sizeof(double) * sendcounts[procRank]);

    A_part = temp_part; // if not to transpose matrix A

    // matrix A transposing
    /*for (int i = 0; i < A_p_width; i++) {
        for (int j = 0; j < A_p_height; j++) {
            A_part[i + j * A_p_width] = temp_part[j + A_p_height * i];
        }
    }*/

    // generating vector b on zero process
    double *b_full = NULL;
    if (procRank == 0)
    {
        b_full = (double *)malloc(sizeof(double) * N);
        for (int i = 0; i < N; i++)
        {
            b_full[i] = (double)(rand() % 100);
            // b[i] = N + 1;
        }
    }

    for (int i = 0; i < procSize; i++)
    {
        sendcounts[i] /= N;
    }

    displs[0] = 0;
    for (int i = 1; i < procSize; i++)
    {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    double *b_part = (double *)malloc(sizeof(double) * sendcounts[procRank]);

    // vector b scattering
    MPI_Scatterv(b_full, sendcounts, displs, MPI_DOUBLE, b_part, sendcounts[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // b dimension
    int b_width = 1;
    int b_part_height = sendcounts[procRank];
    // int b_height = N;

    double *X = (double *)calloc(b_part_height, sizeof(double));
    double *r = (double *)malloc(b_part_height * sizeof(double));
    double *prev_r = (double *)malloc(b_part_height * sizeof(double));
    double *z = (double *)malloc(b_part_height * sizeof(double));

    // square of vector b length
    double square_b_len = 0;
    VectorParallelMult(b_part, b_part, &square_b_len, b_part_height);

    for (int i = 0; i < b_part_height; i++)
    {
        r[i] = b_part[i]; // r0 = b - Ax0, x0 - нулеой вектор
        z[i] = r[i];      // z0 = r0
    }

    int counter = 0;

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

    double tmp_scalar = 0;
    double r_len = 0;
    double *temp = (double *)calloc(b_part_height, sizeof(double));
    double alpha, beta;

    int cntr = 0;
    int cnt = 0;
    double restrictor = epsilon * epsilon * square_b_len;

    do
    {
        // temp = A*z(n)
        // MatrixMultVector(A_part, z, temp, A_p_height, A_p_width);
        parallelMatrixMultVector(A_part, z, temp, A_p_width, A_p_height, procSize, procRank);

        // alpha(n+1) = (r(n), r(n)) / (Az(n), z(n))
        VectorParallelMult(r, r, &r_len, b_part_height);
        VectorParallelMult(temp, z, &tmp_scalar, b_part_height);
        // double oldsqrR = VectorParallelMult(r, r, b_part_height);
        alpha = r_len / tmp_scalar;

        // X(n+1) = X(n) + alpha(n+1)*z(n)
        X = SumVectors(X, ScalarVectMult(z, alpha, b_part_height), b_part_height);

        prev_r = r;
        // r(n+1) = r(n) - alpha(n+1)Az(n)
        r = SubVectors(r, ScalarVectMult(temp, alpha, b_part_height), b_part_height);

        // beta(n+1) = (r(n+1), r(n+1)) / (r(n), r(n))
        double temp_res = 0;
        VectorParallelMult(r, r, &temp_res, b_part_height);
        beta = temp_res / r_len;

        // z(n+1) = r(n+1) + beta(n+1)z(n)
        z = SumVectors(r, ScalarVectMult(z, beta, b_part_height), b_part_height);

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

    double *result = (double *)calloc(N, sizeof(double));
    MPI_Gatherv(X, sendcounts[procRank], MPI_DOUBLE, result, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Print result
    if (procRank == 0)
    {
        printf("\n");
        for (int i = 0; i < N; i++)
        {
            printf("%f ", result[i]);
        }
        printf("\n");
    }

    free(A_full);
    free(b_full);
    free(r);
    free(prev_r);
    free(z);
    free(X);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}