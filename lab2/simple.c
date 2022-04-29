#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>

#define N 1400
#define epsilon 10e-5

double sqrNormOfVector(double *vector)
{
    double result = 0;
    for (int i = 0; i < N; i++)
    {
        result += vector[i] * vector[i];
    }

    return result;
}

double *MatrixMultVector(double *A, double *vector)
{
    double *result = (double *)calloc(N, sizeof(double));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            result[i] += A[i * N + j] * vector[j];
        }
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

int main()
{
    // srand(time(0));

    double *A = (double *)malloc(N * N * sizeof(double));
    double *b = (double *)malloc(N * sizeof(double));
    double *r = (double *)malloc(N * sizeof(double));
    double *prev_r = (double *)malloc(N * sizeof(double));
    double *z = (double *)malloc(N * sizeof(double));
    double *X = (double *)calloc(N, sizeof(double));

    // INITIALIZING OF VARIABLES
    for (int i = 0; i < N; i++)
    {
        b[i] = (double)(rand() % 100);
        r[i] = b[i]; // r0 = b - Ax0, x0 - нулеой вектор
        z[i] = r[i]; // z0 = r0
        X[i] = 0;

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
    double start = clock();

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

    double *temp = NULL;
    double alpha = 0;
    double beta = 0;
    int cntr = 0;
    int cnt = 0;

    double square_b_len = 0; // square of vector b length
    for (int i = 0; i < N; i++)
    {
        square_b_len += b[i] * b[i];
    }

    double restrictor = epsilon * epsilon * square_b_len;

    while (cntr < 3)
    {
        temp = MatrixMultVector(A, z); // A*z(n)

        double oldsqrR = VectorsMult(r, r);
        alpha = oldsqrR / VectorsMult(temp, z); // alpha(n+1) = (r(n), r(n)) / (Az(n), z(n))

        X = SumVectors(X, ScalarVectMult(z, alpha)); // X(n+1) = X(n) + alpha(n+1)*z(n)

        prev_r = r;
        r = SubVectors(r, ScalarVectMult(temp, alpha)); // r(n+1) = r(n) - alpha(n+1)Az(n)

        beta = VectorsMult(r, r) / oldsqrR; // beta(n+1) = (r(n+1), r(n+1)) / (r(n), r(n))

        z = SumVectors(r, ScalarVectMult(z, beta)); // z(n+1) = r(n+1) + beta(n+1)z(n)

        if (sqrNormOfVector(r) < restrictor)
        {
            cntr++;
        }
        cnt++;
    }

    double total_time = (clock() - start) / CLOCKS_PER_SEC;
    printf("Iterations count: %d \n", cnt);
    printf("Total time: %f \n", total_time);

    // Print result
    /*printf("\n");
    for (int i = 0; i < N; i++)
    {
        printf("%f ", X[i]);
    }*/

    free(A);
    free(b);
    free(r);
    free(prev_r);
    free(z);
    free(X);
    return 0;
}