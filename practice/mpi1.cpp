#include <iostream>
#include <ctime>
#include "mpi.h"

void printArray(int *arr)
{
    for (int i = 0; i < 3; ++i)
    {
        std::cout << arr[i] << " ";
    }
    std::cout << "\n\n";
}

void printMatrix(int *arr, int height, int width)
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            std::cout << arr[i * width + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void fillPoints(int *lapTime, int *lapPoints)
{
    if (lapTime[0] < lapTime[1] && lapTime[0] < lapTime[2])
    {
        lapPoints[0] = 10;
        if (lapTime[1] < lapTime[2])
        {
            lapPoints[1] = 8;
            lapPoints[2] = 5;
        }
        else
        {
            lapPoints[1] = 5;
            lapPoints[2] = 8;
        }
    }
    if (lapTime[1] < lapTime[0] && lapTime[1] < lapTime[2])
    {
        lapPoints[1] = 10;
        if (lapTime[0] < lapTime[2])
        {
            lapPoints[0] = 8;
            lapPoints[2] = 5;
        }
        else
        {
            lapPoints[0] = 5;
            lapPoints[2] = 8;
        }
    }
    if (lapTime[2] < lapTime[0] && lapTime[2] < lapTime[1])
    {
        lapPoints[1] = 10;
        if (lapTime[0] < lapTime[1])
        {
            lapPoints[0] = 8;
            lapPoints[1] = 5;
        }
        else
        {
            lapPoints[1] = 5;
            lapPoints[0] = 8;
        }
    }
    if (lapTime[0] < lapTime[2] && lapTime[0] == lapTime[1])
    {
        lapPoints[0] = 10;
        lapPoints[1] = 10;
        lapPoints[2] = 8;
    }
    if (lapTime[1] < lapTime[0] && lapTime[2] == lapTime[1])
    {
        lapPoints[0] = 8;
        lapPoints[1] = 10;
        lapPoints[2] = 10;
    }
    if (lapTime[0] < lapTime[1] && lapTime[0] == lapTime[2])
    {
        lapPoints[0] = 10;
        lapPoints[1] = 8;
        lapPoints[2] = 10;
    }
    if (lapTime[0] == lapTime[1] && lapTime[0] == lapTime[2])
    {
        lapPoints[0] = 10;
        lapPoints[1] = 10;
        lapPoints[2] = 10;
    }
}

int main(int argc, char **argv)
{
    // 5 processes!!!

    srand(time(0));

    int procSize;
    int procRank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);

    int bloknot[3 * 30];

    if (procRank == 0)
    {
        for (int i = 0; i < 30; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                bloknot[i * 3 + j] = (rand() % 7 != 0) ? ((3 + rand()) % 10) : 0;
            }
        }
    }

    if (procRank == 0)
        printMatrix(bloknot, 30, 3);

    int lapTime[3] = {0, 0, 0};
    int lapPoints[3] = {0, 0, 0};
    int resultPoints[3] = {0, 0, 0};

    int *rowsNum = (int *)calloc(procSize, sizeof(int));
    // number of elements in each process part
    int *sendcounts = (int *)calloc(procSize, sizeof(int));
    for (int i = 0; i < procSize; i++)
    {
        rowsNum[i] = 6;
        sendcounts[i] = 18;
    }

    int *bloknot_part = (int *)malloc(sizeof(double) * 18);
    MPI_Scatter(bloknot, 18, MPI_INT, bloknot_part, 18, MPI_INT, 0, MPI_COMM_WORLD);

    for (int j = 0; j < 3; ++j)
    {
        lapTime[j] = 0;
        for (int k = 0; k < 6; ++k)
        {
            lapTime[j] += bloknot_part[k * 3 + j];
        }
    }

    fillPoints(lapTime, lapPoints);

    for (int j = 0; j < 3; ++j)
    {
        resultPoints[j] += lapPoints[j];
    }

    int result[3] = {0, 0, 0};
    // MPI_Allreduce(resultPoints, result, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int j = 0; j < 3; ++j)
    {
        MPI_Allreduce(&resultPoints[j], &result[j], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    int given_places[3] = {0, 0, 0};
    int first_place = 10000;
    int max_points = 0;

    if (procRank == 0)
    {
        printArray(result);
    }

    if (procRank == 0)
    {
        for (int k = 0; k < 3; ++k)
        {
            max_points = 0;
            for (int i = 0; i < 3; i++)
            {
                if (result[i] >= max_points && given_places[i] == 0)
                {
                    max_points = result[i];
                    first_place = i;
                }
            }
            given_places[first_place] = k + 1;
        }
    }

    if (procRank == 0)
    {
        printArray(given_places);
    }

    MPI_Finalize();

    return 0;
}