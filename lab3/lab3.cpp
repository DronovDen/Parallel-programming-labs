#include <iostream>
#include <ctime>
#include <cstdlib>

#define HEIGHT1 4      // 24     // 1500
#define WIDTH_HEIGHT 4 // 8      // 1500
#define WIDTH2 4       // 9      // 1500
#define P1 2           // 6
#define P2 2           // 3

using namespace std;

void showMatrix(float *matrix, size_t height, size_t width)
{
    size_t i, j;
    for (i = 0; i < height; ++i)
    {
        for (j = 0; j < width; ++j)
        {
            cout << matrix[i * width + j] << " ";
        }
        cout << endl;
    }
    cout << "\n"
         << "\n";
}

void multiplicateMatrix(float *M1, float *M2, float *result, size_t n1, size_t n2, size_t n3)
{
    size_t i, k, j;
    for (i = 0; i < n1; ++i)
    {
        for (k = 0; k < n2; ++k)
        {
            for (j = 0; j < n3; ++j)
            {
                result[i * n3 + j] += M1[i * n2 + k] * M2[k * n3 + j];
                // tmp++;
            }
        }
    }
}

/* void multiplicateMatrix(float *A, float *B, float *C, size_t n1, size_t n2, size_t n3)
{
    for (int i = 0; i < n1; ++i)
    {
        float *c = C + i * n2;
        for (int j = 0; j < n2; ++j)
            c[j] = 0;
        for (int k = 0; k < n3; ++k)
        {
            const float *b = B + k * n2;
            float a = A[i * n3 + k];
            for (int j = 0; j < n2; ++j)
                c[j] += a * b[j];
        }
    }
} */

float *generateMatrix(size_t height, size_t width)
{
    float *tmp = (float *)calloc(height * width, sizeof(float));
    for (size_t i = 0; i < height * width; i++)
    {
        tmp[i] = rand() % 10;
    }
    return tmp;
}

int main()
{
    // srand(time(0));
    //  float *A = new float[N * N];
    //  generateMatrix(A, N);
    /*float A[N * N] = {2, 4,
                      2, 3};*/

    float A[16] = {1, 2, 4, 0,
                   4, 4, 3, 3,
                   2, 4, 0, 0,
                   1, 2, 1, 1};

    float B[16] = {0, 2, 2, 1,
                   1, 4, 2, 3,
                   2, 2, 1, 1,
                   3, 0, 2, 1};

    /*float A[N * N] = {2, 4, 2, 3, 2, 1, 5, 0,
                      8, 9, 2, 1, 1, 4, 6, 7,
                      1, 3, 5, 7, 8, 9, 2, 1,
                      7, 9, 2, 5, 1, 4, 6, 7,
                      1, 4, 6, 7, 1, 5, 5, 7,
                      2, 4, 2, 3, 2, 4, 2, 3,
                      1, 4, 6, 7, 1, 3, 5, 7,
                      0, 4, 2, 3, 2, 4, 2, 6};*/

    /* float *A = generateMatrix(HEIGHT1, WIDTH_HEIGHT);
    showMatrix(A, HEIGHT1, WIDTH_HEIGHT);
    float *B = generateMatrix(WIDTH_HEIGHT, WIDTH2);
    showMatrix(B, WIDTH_HEIGHT, WIDTH2); */

    float *C = (float *)calloc(HEIGHT1 * WIDTH2, sizeof(float));
    multiplicateMatrix(A, B, C, HEIGHT1, WIDTH_HEIGHT, WIDTH2);
    showMatrix(C, HEIGHT1, WIDTH2);

    // delete[] A;
    float total_time = clock() / CLOCKS_PER_SEC;
    cout << "Total time: " << total_time << " sec." << endl;

    // free(A);
    // free(B);
    free(C);
    return 0;
}