#include <iostream>
#include <ctime>

int main()
{
    srand(time(0));
    int sum = 0;
    int iterations_count[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
    {
        while (sum < 1000000)
        {
            int number = (rand() % 25 + 10) * 10000; // from 100000 to 350000
            std::cout << number << " ";
            sum += number;
            iterations_count[i]++;
        }
        sum = 0;
        std::cout << "\n";
    }

    for (int i = 0; i < 3; i++)
    {
        std::cout << iterations_count[i] << " ";
    }

    return 0;
}