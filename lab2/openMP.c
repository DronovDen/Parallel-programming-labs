#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <omp.h>

#define N 1400
#define epsilon 10e-5

int main(int argc, char **argv)
{
    // srand(time(0));

    double *A = (double *)malloc(N * N * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *prev_r = (double *)malloc(N * sizeof(double));
    double *z = (double *)malloc(N * sizeof(double));
    double *X = (double *)calloc(N, sizeof(double));
    double *Az = (double *)malloc(N * sizeof(double));

    omp_set_num_threads(atoi(argv[1]));

    // INITIALIZING OF VARIABLES
    for (int i = 0; i < N; i++)
    {
        b[i] = (double)(rand() % 100);

        for (int j = i; j < N; j++)
        {
            if (i == j)
            {
                A[i * N + j] = (double)(rand() % 100) + 200;
            }
            else
            {
                A[i * N + j] = A[j * N + i] = (double)(rand() % 100);
            }
        }
    }

    double start = omp_get_wtime();

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        X[i] = 0;
    }

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        r[i] = b[i]; // r0 = b - Ax0, x0 - нулеой вектор
        z[i] = r[i]; // z0 = r0
    }

    // Print matrix A, b, X
    /*for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            printf("%f ", A[i * N + j]);
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

    double square_b_len = 0; // square of vector b length
#pragma omp parallel for reduction(+ \
                                   : square_b_len)
    for (int i = 0; i < N; i++)
    {
        square_b_len += b[i] * b[i];
    }

    double restrictor = epsilon * epsilon * square_b_len;

    double temp = 0;
    double alpha = 0;
    double beta = 0;
    int cntr = 0;
    int iterations = 0;
    double sqrR = 0;
    double oldsqrR = 0;

    // double oldsqrR = VectorsMult(r, r);
#pragma omp parallel for reduction(+ \
                                   : oldsqrR)
    for (int i = 0; i < N; i++)
    {
        oldsqrR += r[i] * r[i];
    }

#pragma omp parallel
    while (cntr < 3)
    {
#pragma omp single
        {
            sqrR = 0;
            temp = 0;
        }

        // A*z(n)
#pragma omp for
        for (int i = 0; i < N; i++)
        {
            Az[i] = 0;
            for (int j = 0; j < N; j++)
            {
                Az[i] += A[i * N + j] * z[j];
            }
        }

#pragma omp for reduction(+ \
                          : temp)
        for (int i = 0; i < N; i++)
        {
            temp += Az[i] * z[i];
        }

        // alpha(n+1) = (r(n), r(n)) / (Az(n), z(n))
        alpha = oldsqrR / temp;

        // X(n+1) = X(n) + alpha(n+1)*z(n)
#pragma omp for
        for (int i = 0; i < N; i++)
        {
            X[i] = X[i] + alpha * z[i];
        }

        //  r(n+1) = r(n) - alpha(n+1)Az(n)
#pragma omp for
        for (int i = 0; i < N; i++)
        {
            r[i] = r[i] - alpha * Az[i];
        }

        // beta(n+1) = (r(n+1), r(n+1)) / (r(n), r(n))
#pragma omp for reduction(+ \
                          : sqrR)
        for (int i = 0; i < N; i++)
        {
            sqrR += r[i] * r[i];
        }
#pragma omp single
        {
            beta = sqrR / oldsqrR;
            oldsqrR = sqrR;
        }

        // z(n+1) = r(n+1) + beta(n+1)z(n)
#pragma omp for
        for (int i = 0; i < N; i++)
        {
            z[i] = r[i] + beta * z[i];
        }

#pragma omp single
        if (sqrR < restrictor)
        {
            cntr++;
        }
#pragma omp single
        {
            iterations++;
        }
    }

    double finish = omp_get_wtime();
    double total_time = (finish - start);
    printf("Iterations count: %d \n", iterations);
    printf("Total time: %f \n", total_time);

    // Print result
    printf("\n");
    for (int i = 0; i < N; i++)
    {
        printf("%f ", X[i]);
    }

    free(A);
    free(b);
    free(r);
    free(z);
    free(X);
    return 0;
}