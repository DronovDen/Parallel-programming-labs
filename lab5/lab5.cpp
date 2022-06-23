#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

#define WEIGHT_FACTOR 10000
#define TASKS_PER_PROCESS 1000
#define MIN_TASKS_TO_SHARE 10
#define ITERATIONS_COUNT 5
#define REQUEST 0
#define REPLY 1

double global_result = 0.0;
double global_result_sum = 0.0;
int procRank, procSize;

int *tasks;
int tasks_remaining;

pthread_mutex_t mutex_tasks;
pthread_mutex_t mutex_tasks_remaining;

pthread_t thread_receiver;

void setTasksWeight(int *tasks, int count, int iterCounter)
{
    pthread_mutex_lock(&mutex_tasks);
    for (int i = 0; i < count; ++i)
    {
        tasks[i] = abs(50 - i % 100) * abs(procRank - (iterCounter % procSize)) * WEIGHT_FACTOR;
    }
    pthread_mutex_unlock(&mutex_tasks);
}

void completeTasks(int *tasks, int *tasks_executed)
{
    pthread_mutex_lock(&mutex_tasks_remaining);

    for (int i = 0; tasks_remaining; ++i, --tasks_remaining)
    {
        pthread_mutex_unlock(&mutex_tasks_remaining);

        pthread_mutex_lock(&mutex_tasks);
        int task_weight = tasks[i];
        pthread_mutex_unlock(&mutex_tasks);

        for (int j = 0; j < task_weight; ++j)
        {
            global_result += sin(j);
        }

        ++(*tasks_executed);
        pthread_mutex_lock(&mutex_tasks_remaining);
    }
    pthread_mutex_unlock(&mutex_tasks_remaining);
}

void *receiverThreadStart(void *args)
{
    int tasks_to_share;
    int rankRequestedTasks;

    while (true)
    {
        // receiving process rank that requests tasks
        MPI_Recv(&rankRequestedTasks, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (rankRequestedTasks == procRank)
            break;

        pthread_mutex_lock(&mutex_tasks_remaining);
        if (tasks_remaining >= MIN_TASKS_TO_SHARE)
        {
            tasks_to_share = tasks_remaining / 2;
            tasks_remaining -= tasks_to_share;

            // sending number of tasks that shares
            MPI_Send(&tasks_to_share, 1, MPI_INT, rankRequestedTasks, REPLY, MPI_COMM_WORLD);

            pthread_mutex_lock(&mutex_tasks);
            // sending tasks
            MPI_Send(tasks + tasks_remaining - 1, tasks_to_share, MPI_INT, rankRequestedTasks, REPLY, MPI_COMM_WORLD);
            pthread_mutex_unlock(&mutex_tasks);
        }
        else
        {
            tasks_to_share = 0;

            MPI_Send(&tasks_to_share, 1, MPI_INT, rankRequestedTasks, REPLY, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex_tasks_remaining);
    }
    return NULL;
}

void *workerThreadStart(void *args)
{
    tasks = (int *)malloc(sizeof(int) * TASKS_PER_PROCESS);

    double start_time, total_time;
    double minTime, maxTime;

    for (int iterCounter = 0; iterCounter < ITERATIONS_COUNT; ++iterCounter)
    {
        setTasksWeight(tasks, TASKS_PER_PROCESS, iterCounter);

        pthread_mutex_lock(&mutex_tasks_remaining);
        tasks_remaining = TASKS_PER_PROCESS;
        pthread_mutex_unlock(&mutex_tasks_remaining);
        int tasks_executed = 0;
        int tasks_received;

        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();

        completeTasks(tasks, &tasks_executed);

        // Start requesting tasks from other processes
        //(sequental iterating)
        for (int rank_iter = 0; rank_iter < procSize; ++rank_iter)
        {

            if (rank_iter == procRank)
                continue;

            // Notificating other process that current process is free
            // and want to execute extra tasks (sending other process current procRank)
            MPI_Send(&procRank, 1, MPI_INT, rank_iter, REQUEST, MPI_COMM_WORLD);

            // Receive number of extra tasks
            MPI_Recv(&tasks_received, 1, MPI_INT, rank_iter, REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (tasks_received > 0)
            {
                // Receive tasks
                MPI_Recv(tasks, tasks_received, MPI_INT, rank_iter, REPLY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                pthread_mutex_lock(&mutex_tasks_remaining);
                tasks_remaining = tasks_received;
                pthread_mutex_unlock(&mutex_tasks_remaining);

                completeTasks(tasks, &tasks_executed);
            }
        }

        total_time = MPI_Wtime() - start_time;

        MPI_Allreduce(&total_time, &minTime, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

        MPI_Allreduce(&total_time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        // Print info about this iteration
        if (procRank == 0)
        {
            printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
            printf("Iteration %d\n", iterCounter + 1);
            printf("Disbalance time: %lf sec.\n", maxTime - minTime);
            printf("Disbalance percentage: %lf %%\n", (maxTime - minTime) / maxTime * 100);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        for (int rank_iter = 0; rank_iter < procSize; rank_iter++)
        {
            if (procRank == rank_iter)
            {
                printf("\tProcess #%d\n", procRank + 1);
                printf("\t\texecuted tasks:\t%d\n", tasks_executed);
                printf("\t\tglobal result:\t%lf\n", global_result);
                printf("\t\titeration time:\t%lf sec.\n", total_time);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    // Terminate Thread 'Receiver'
    MPI_Send(&procRank, 1, MPI_INT, procRank, 0, MPI_COMM_WORLD);

    MPI_Allreduce(&global_result, &global_result_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    free(tasks);
    return NULL;
}

void createThreads()
{
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // creating joinable thread for receiving requests from worker thread
    pthread_create(&thread_receiver, &attr, receiverThreadStart, NULL);
    pthread_attr_destroy(&attr);

    // starting work on main thread
    workerThreadStart(NULL);

    // waitng for recieiver thread to finish
    pthread_join(thread_receiver, NULL);
}

int main(int argc, char **argv)
{
    int required = MPI_THREAD_MULTIPLE;
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != required)
    {
        fprintf(stderr, "Provided level doesn't match the required!\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);

    pthread_mutex_init(&mutex_tasks, NULL);
    pthread_mutex_init(&mutex_tasks_remaining, NULL);

    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    createThreads();
    double total_time = MPI_Wtime() - start_time;

    pthread_mutex_destroy(&mutex_tasks);
    pthread_mutex_destroy(&mutex_tasks_remaining);

    if (procRank == 0)
    {
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        printf("Total time:\t%lf sec.\n", total_time);
        printf("Global result sum:\t%lf\n", global_result_sum);
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}