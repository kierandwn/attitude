#include <cstdlib>
#include <iostream>
#include <ctime>

#include "include/helper.h"

int random_not(int i)
{
    int j = (std::rand() / RAND_MAX);
    int HI, LO;

    if (i == 0)
    {
        HI = 2;
        LO = 1;
    }
    else if (i == 1)
    {
        HI = 2;
        LO = 0;
    }
    else if (i == 2)
    {
        HI = 1;
        LO = 0;
    }
    return (j == 1) ? HI : LO;
}

int randomise_order()
{
    int i = 2 * (std::rand() / RAND_MAX);
    int j = random_not(i);
    int k = random_not(j);

    return i * 100 + j * 10 + k;
}