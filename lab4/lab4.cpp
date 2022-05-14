#include <iostream>
#include "mpi.h"

#define WIDTH 200
#define HEIGHT 200
//#define DEAD false
//#define ALIVE true

enum State
{
    DEAD,
    ALIVE
}

void
initField(bool *field, int height, int width)
{
    for (int i = 0; i < height * width; ++i)
    {
        field[i] = DEAD;
    }

    field[2 * WIDTH + 0] = ALIVE;
    field[2 * WIDTH + 1] = ALIVE;
    field[2 * WIDTH + 2] = ALIVE;
    field[1 * WIDTH + 2] = ALIVE;
    field[0 * WIDTH + 1] = ALIVE;
}

void calculateStopFlags(bool *fieldPart, bool **history, int iterCounter, int width, int height, bool *local_similarity)
{
    for (int i = 0; i < iterCounter; i++)
    {
        local_similarity[i] = false;
    }

    int simiarityCount = 0;
    bool exitFlag = false;
    for (int i = 0; i < iterCounter; i++)
    {
        bool *part = history[i]; // itearating through field parts history
        for (int j = 0; j < height; j++)
        {
            for (int k = 0; k < width; k++)
            {
                if (part[j * width + k] == fieldPart[j * width + k])
                {
                    simiarityCount++;
                }
                else
                {
                    exitFlag = true;
                    break;
                }
            }
            if (exitFlag)
            {
                exitFlag = false;
                break;
            }
        }

        if (simiarityCount == width * height) // cuurent field state repeats one of the previous ones
        {
            local_similarity[i] = true;
        }
        else
        {
            local_similarity[i] = false;
            simiarityCount = 0;
        }
    }
}

void simualteEra(bool *part, int width, int height, bool *firstRow, bool *penultRow)
{
    int *aliveAround = new int[width * height];
    for (int i = 0; i < width * height; ++i)
    {
        aliveAround[i] = 0;
    }

    for (int i = 0; i < width; i++)
    {
        firstRow[i] = part[width + i];                 // first row (zero row is before)
        penultRow[i] = part[width * (height - 2) + i]; // penultimate row
    }

    for (int i = 1; i < height - 1; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int idx = i * width + j;
            if (j == (width - 1))
            {
                aliveAround[idx] += int(part[idx - 1]) + int(part[idx + 1]) + int(part[idx + width]) + int(part[idx - width]) + int(part[idx + width - 1]) + int(part[idx - width - 1]) + int(part[idx + 1 - width]) + int(part[idx + 1 - 2 * width]);
            }
            else if (j == 0)
            {
                aliveAround[idx] += int(part[idx - 1]) + int(part[idx + 1]) + int(part[idx + width]) + int(part[idx - width]) + int(part[idx + 1 - width]) + int(part[idx + 1 + width]) + int(part[idx + 2 * width - 1]) + int(part[idx - 1 + width]);
            }
            else
            {
                aliveAround[idx] += int(part[idx - 1]) + int(part[idx + 1]) + int(part[idx + 1 + width]) + int(part[idx - 1 + width]) + int(part[idx + width]) + int(part[idx - width]) + int(part[idx + 1 - width]) + int(part[idx - 1 - width]);
            }
        }
    }
    for (int i = 1; i < height - 1; i++)
    {
        for (int j = 0; j < width; j++)
        {
            int idx = i * width + j;
            if (part[idx] == ALIVE)
            {
                if (aliveAround[idx] > 3 || aliveAround[idx] < 2)
                {
                    part[idx] = DEAD;
                }
                else
                {
                    part[idx] = ALIVE;
                }
            }
            else
            {
                if (aliveAround[idx] == 3)
                {
                    part[idx] = ALIVE;
                }
            }
        }
    }
    delete[] aliveAround;
}

void simulateRemainingCells(bool *part, bool *topLine, bool *bottomLine, int width, int displ)
{
    int *aliveAround = new int[width];
    for (int i = 0; i < width; ++i)
    {
        aliveAround[i] = 0;
    }

    for (int i = 0; i < width; i++)
    {
        if (i == 0)
        {
            aliveAround[i] += int(topLine[i]) + int(topLine[i + 1]) + int(topLine[width - 1]) + int(bottomLine[i]) + int(bottomLine[i + 1]) + int(bottomLine[width - 1]) + int(part[i + 1 + displ]) + int(part[width - 1 + displ]);
        }
        else if (i == width - 1)
        {
            aliveAround[i] += int(topLine[i]) + int(topLine[i - 1]) + int(topLine[0]) + int(bottomLine[i]) + int(bottomLine[i - 1]) + int(bottomLine[0]) + int(part[i - 1 + displ]) + int(part[0 + displ]);
        }
        else
        {
            aliveAround[i] += int(topLine[i]) + int(topLine[i - 1]) + int(topLine[i + 1]) + int(bottomLine[i]) + int(bottomLine[i - 1]) + int(bottomLine[i + 1]) + int(part[i - 1 + displ]) + int(part[i + 1 + displ]);
        }
    }

    for (int i = 0; i < width; i++)
    {
        if (part[i + displ] == ALIVE)
        {
            if (aliveAround[i] > 3 || aliveAround[i] < 2)
            {
                part[i + displ] = DEAD;
            }
            else
            {
                part[i + displ] = ALIVE;
            }
        }
        else
        {
            if (aliveAround[i] == 3)
            {
                part[i + displ] = ALIVE;
            }
        }
    }
    delete[] aliveAround;
}

bool allPartsRepeat(bool *general_similarity, int procSize, int iterCounter, int iterations)
{
    for (int i = 0; i < iterCounter; i++)
    {
        int counter = 0;
        for (int j = 0; j < procSize; j++)
        {
            if (general_similarity[j * iterations + i] == true)
            {
                counter++;
            }
            else
            {
                break;
            }
        }
        if (counter == procSize) // every process have duplicate part on a certain iteration
        {
            return true;
        }
        else
        {
            counter = 0;
        }
    }
    return false;
}

int startSimulation(int procSize, int procRank)
{
    int *rowsNum = new int[procSize];
    int *sendcounts = new int[procSize];
    int *displs = new int[procSize];

    int rowsDefault = (HEIGHT / procSize);
    int dif = HEIGHT - (rowsDefault * procSize);
    int displ = 0;

    for (int i = 0; i < procSize; i++)
    {
        rowsNum[i] = rowsDefault;
        if (i < dif)
        {
            rowsNum[i] += 1;
        }
        sendcounts[i] = rowsNum[i] * WIDTH;
        displs[i] = displ;

        displ += sendcounts[i];
    }

    bool *field = NULL;
    bool *fieldPart = new bool[rowsNum[procRank] * WIDTH];

    bool *topLine = new bool[WIDTH];
    bool *bottomLine = new bool[WIDTH];

    bool *firstRow = new bool[WIDTH];
    bool *penultRow = new bool[WIDTH];

    // array that stores history of field parts states for each iteration
    int itersRestrictor = 2000;
    bool **history = new bool *[itersRestrictor];
    for (int i = 0; i < itersRestrictor; i++)
    {
        history[i] = new bool[rowsNum[procRank] * WIDTH];
    }

    bool *general_similarity = new bool[procSize * itersRestrictor];
    bool *local_similarity = new bool[itersRestrictor];
    for (int i = 0; i < itersRestrictor; i++)
    {
        local_similarity[i] = false;
    }

    if (procRank == 0)
    {
        field = new bool[HEIGHT * WIDTH];
        initField(field, HEIGHT, WIDTH);
    }

    MPI_Scatterv(field, sendcounts, displs, MPI_C_BOOL, fieldPart, sendcounts[procRank], MPI_C_BOOL, 0, MPI_COMM_WORLD);

    int next_rank = procRank + 1;
    int prev_rank = procRank - 1;
    if (next_rank == procSize)
    {
        next_rank = 0;
    }
    if (prev_rank == -1)
    {
        prev_rank = procSize - 1;
    }

    int iterCounter = 0;
    while (iterCounter < itersRestrictor)
    {
        MPI_Request sendRequestUp;
        MPI_Request sendRequestDown;
        MPI_Request recvRequestUp;
        MPI_Request recvRequestDown;
        MPI_Request recvRequsetVectors;

        MPI_Isend(fieldPart, WIDTH, MPI_C_BOOL, prev_rank, 1, MPI_COMM_WORLD, &sendRequestUp);
        MPI_Isend(fieldPart + WIDTH * (rowsNum[procRank] - 1), WIDTH, MPI_C_BOOL, next_rank, 0, MPI_COMM_WORLD, &sendRequestDown);

        MPI_Irecv(topLine, WIDTH, MPI_C_BOOL, prev_rank, 0, MPI_COMM_WORLD, &recvRequestUp);      // receiving upper line from previous process
        MPI_Irecv(bottomLine, WIDTH, MPI_C_BOOL, next_rank, 1, MPI_COMM_WORLD, &recvRequestDown); // receiving bottom line from next process
        memcpy(history[iterCounter], fieldPart, sendcounts[procRank]);

        calculateStopFlags(fieldPart, history, iterCounter, WIDTH, rowsNum[procRank], local_similarity);

        MPI_Iallgather(local_similarity, itersRestrictor, MPI_C_BOOL, general_similarity, itersRestrictor, MPI_C_BOOL, MPI_COMM_WORLD, &recvRequsetVectors);

        simualteEra(fieldPart, WIDTH, rowsNum[procRank], firstRow, penultRow);

        // Simulating remainig cells of field part using lines recieved from top and bottom processes
        MPI_Wait(&sendRequestUp, MPI_STATUS_IGNORE);
        MPI_Wait(&recvRequestUp, MPI_STATUS_IGNORE);
        simulateRemainingCells(fieldPart, topLine, firstRow, WIDTH, 0);

        MPI_Wait(&sendRequestDown, MPI_STATUS_IGNORE);
        MPI_Wait(&recvRequestDown, MPI_STATUS_IGNORE);
        simulateRemainingCells(fieldPart, penultRow, bottomLine, WIDTH, WIDTH * (rowsNum[procRank] - 1));

        MPI_Wait(&recvRequsetVectors, MPI_STATUS_IGNORE);
        if (allPartsRepeat(general_similarity, procSize, iterCounter, itersRestrictor))
        {
            break;
        }

        iterCounter++;
    }

    delete[] bottomLine;
    delete[] fieldPart;
    delete[] topLine;
    delete[] local_similarity;
    delete[] general_similarity;
    delete[] firstRow;
    delete[] penultRow;
    for (int i = 0; i < itersRestrictor; i++)
    {
        delete[] history[i];
    }
    delete[] history;

    return iterCounter;
}

int main(int argc, char **argv)
{
    int procSize;
    int procRank;

    double start, finish;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    start = MPI_Wtime();
    int iterationsNum = startSimulation(procSize, procRank);
    finish = MPI_Wtime();

    if (procRank == 0)
    {
        std::cout << "Total time: " << finish - start << " sec." << std::endl;
        std::cout << "Iterations count: " << iterationsNum << std::endl;
    }

    MPI_Finalize();
    return 0;
}
